"""Ported from cna_maste `python/cnamaste/spatio_genomic_counts.py` at 42b96ff (2026-11-01):
the definitions this repository's analysis reaches, verbatim less the methods it never calls, with their imports
pointed at `scripts.cnamaste_port`.

Per-spot arrays `(n_obs, n_spots)`, dense or CSR: `SpatioGenomicCounts`, and `single_X` as `SpatioGenomicPairCounts`.

`SpatioGenomicCounts` is one per-spot quantity, `(n_obs, n_spots)`: a read depth, an allele count, a baseline.
`SpatioGenomicPairCounts` is `single_X`, `(n_obs, 2, n_spots)`: two `SpatioGenomicCounts`, read depth (`rdr`, channel 0)
and B allele count (`baf`, channel 1), in one backing. `.rdr` and `.baf` return the `SpatioGenomicCounts`.

A run builds them as CSR; a capture records them dense, and its replays hand them back dense.
Every operation computes the same bits on both backings: sums in numpy's order (segments in sequence on `axis=0`,
spots pairwise on `axis=1`).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import scipy.sparse as sp
from numba import njit, prange


@njit(cache=True)
def _pairwise_leaf(idx, val, lo, hi, a, n):
    """Sum positions `[a, a + n)`, at most 128, bitwise as numpy's `pairwise_sum` leaf.

    Below 8 in sequence, else eight interleaved partial sums combined as numpy does, then the remainder; absent is 0.0.

    Args: `idx`, `val` a row's sorted stored indices and values; `lo`, `hi` its slice; `a`, `n` the dense range.
    Returns: the float sum.
    """
    if n < 8:
        res = 0.0
        for p in range(lo, hi):
            res += val[p]
        return res
    r = np.zeros(8)
    m = n - (n % 8)
    for p in range(lo, hi):
        j = idx[p] - a
        if j < m:
            r[j % 8] += val[p]
    res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]))
    for p in range(lo, hi):
        if idx[p] - a >= m:
            res += val[p]
    return res


@njit(cache=True)
def _pairwise(idx, val, lo, hi, a, n):
    """Sum dense row positions `[a, a + n)` from stored entries, bitwise as numpy's `pairwise_sum`.

    Blocks of at most 128 by `_pairwise_leaf`; larger ranges split at half rounded down to a multiple of 8.
    Iterative with an explicit stack: numba's cache does not reload a self-recursive function intact.

    Args: `idx`, `val` sorted stored indices and values; `lo`, `hi` the row's slice; `a`, `n` the dense range.
    Returns: the float sum.
    """
    depth = 64
    fa, fn, flo, fhi = (
        np.zeros(depth, np.int64),
        np.zeros(depth, np.int64),
        np.zeros(depth, np.int64),
        np.zeros(depth, np.int64),
    )
    fmid, stage = np.zeros(depth, np.int64), np.zeros(depth, np.int64)
    res = np.zeros(depth)
    top, rtop = 1, 0
    fa[0], fn[0], flo[0], fhi[0], stage[0] = a, n, lo, hi, 0
    while top > 0:
        t = top - 1
        if fn[t] <= 128:
            res[rtop] = _pairwise_leaf(idx, val, flo[t], fhi[t], fa[t], fn[t])
            rtop += 1
            top -= 1
        elif stage[t] == 0:
            n2 = fn[t] // 2
            n2 -= n2 % 8
            mid = flo[t]
            while mid < fhi[t] and idx[mid] < fa[t] + n2:
                mid += 1
            fmid[t], stage[t] = mid, 1
            fa[top], fn[top], flo[top], fhi[top], stage[top] = fa[t], n2, flo[t], mid, 0
            top += 1
        elif stage[t] == 1:
            n2 = fn[t] // 2
            n2 -= n2 % 8
            stage[t] = 2
            fa[top], fn[top], flo[top], fhi[top], stage[top] = fa[t] + n2, fn[t] - n2, fmid[t], fhi[t], 0
            top += 1
        else:
            res[rtop - 2] = res[rtop - 2] + res[rtop - 1]
            rtop -= 1
            top -= 1
    return res[0]


@njit(parallel=True, cache=True)
def _row_sums(indptr, indices, data, n_cols):
    out = np.empty(len(indptr) - 1)
    for i in prange(len(indptr) - 1):
        out[i] = _pairwise(indices, data, indptr[i], indptr[i + 1], 0, n_cols)
    return out


@njit(parallel=True, cache=True)
def _row_sums_in_sequence(indptr, indices, data):
    """Sum each CSR row's stored entries in column order from 0.0, as numpy reduces a non-contiguous axis.

    Args: `indptr`, `indices`, `data` a CSR.
    Returns: the per-row sums.
    """
    out = np.zeros(len(indptr) - 1)
    for i in prange(len(indptr) - 1):
        acc = 0.0
        for p in range(indptr[i], indptr[i + 1]):
            acc += data[p]
        out[i] = acc
    return out


@njit(cache=True)
def _column_sums(indptr, indices, data, n_cols):
    """Compute `np.sum(axis=0)` of a C-order float array from its CSR, each column over the rows in sequence.

    Args: `indptr`, `indices`, `data` a CSR; `n_cols` its width.
    Returns: the per-column sums.
    """
    out = np.zeros(n_cols)
    for i in range(len(indptr) - 1):
        for p in range(indptr[i], indptr[i + 1]):
            out[indices[p]] += data[p]
    return out


def canonical(m: Any, dtype: Any = None) -> sp.csr_matrix:
    """Return `m` as canonical CSR, with no explicit zeros, cast to `dtype` where given.

    Args: `m` a sparse matrix; `dtype` optional.
    Returns: the canonical CSR.
    """
    out = sp.csr_matrix(m, dtype=dtype if dtype is not None else m.dtype)
    out.eliminate_zeros()
    out.sort_indices()
    out.indptr, out.indices = out.indptr.astype(np.int64), out.indices.astype(np.int64)
    return out


class SpatioGenomicCounts:
    """One per-spot quantity `(n_obs, n_spots)`, dense or canonical CSR; every operation gives the same bits from either.

    A `SpatioGenomicPairCounts` holds two of them, its read-depth and allele channels.

    Args: `values` a dense array or sparse matrix (made canonical CSR).
    """

    ndim = 2

    def __init__(self, values: Any) -> None:
        # NB a dense array kept as given: `sum` keeps numpy's order.
        self.values = canonical(values) if sp.issparse(values) else np.asarray(values)

    @classmethod
    def of(cls, a: Any) -> SpatioGenomicCounts:
        return a if isinstance(a, cls) else cls(a)

    @property
    def sparse(self) -> bool:
        return sp.issparse(self.values)

    @property
    def shape(self) -> tuple[int, int]:
        return self.values.shape

    @property
    def dtype(self) -> np.dtype:
        return self.values.dtype

    def sum(self, axis: int | None = None) -> Any:
        """Compute `np.sum(values, axis)` bitwise on either backing.

        Dense summed by numpy in its layout's order; CSR as numpy sums the C-order array: integers in any order, floats
        over segments in sequence (`axis=0`), spots pairwise (`axis=1`). A run's arrays are C-order, so both agree.

        Args: `axis` 0, 1 or None.
        Returns: the sum.
        """
        if not self.sparse:
            return np.sum(self.values, axis=axis)
        m = self.values
        if self.dtype.kind in "iub":
            return (
                m.sum(axis=axis, dtype=self.dtype)
                if axis is None
                else np.asarray(m.sum(axis=axis, dtype=self.dtype)).ravel()
            )
        if axis == 0:
            return _column_sums(m.indptr, m.indices, m.data.astype(np.float64), m.shape[1])
        if axis == 1:
            return _row_sums(m.indptr, m.indices, m.data.astype(np.float64), m.shape[1])
        raise ValueError("sum of float CSR over both axes: give an axis.")

    def spot_sums(self, idx: Any = None) -> np.ndarray:
        """Compute `np.sum(values[:, idx], axis=1)` bitwise: per segment, spots `idx` summed in sequence, `idx`'s order.

        numpy's advanced indexing returns the selection F-ordered, so it reduces in sequence, not pairwise.
        Every spot, summed as `sum(axis=1)`, where `idx` is None.

        Args: `idx` spot indices or None.
        Returns: the (n_obs,) sums.
        """
        if idx is None:
            return self.sum(axis=1)
        if not self.sparse:
            return np.sum(np.asfortranarray(self.values[:, idx]), axis=1)
        sub = canonical(self.values[:, idx])
        if self.dtype.kind in "iub":
            return np.asarray(sub.sum(axis=1, dtype=self.dtype)).ravel()
        return _row_sums_in_sequence(sub.indptr, sub.indices, sub.data.astype(np.float64))

    def count_positive(self, axis: int | None = None) -> Any:
        """Compute `(values > 0).sum(axis)`.

        Args: `axis` 0, 1 or None.
        Returns: the counts of positive entries.
        """
        if not self.sparse:
            return (self.values > 0).sum(axis=axis)
        m = self.values
        positive = m.data > 0
        if axis is None:
            return int(np.count_nonzero(positive))
        if axis == 0:
            return np.bincount(m.indices[positive], minlength=m.shape[1])
        return (
            np.add.reduceat(positive.astype(np.int64), m.indptr[:-1]) * (np.diff(m.indptr) > 0)
            if m.nnz
            else np.zeros(m.shape[0], np.int64)
        )


class SpatioGenomicPairCounts:
    """`single_X` `(n_obs, 2, n_spots)`: channel 0 read depth (`rdr`), 1 BAF numerator (`baf`), one backing.

    Each channel is a `SpatioGenomicCounts`, one per-spot quantity; `SpatioGenomicPairCounts` pairs them as the HMM's
    two observations.

    Args: `rdr`, `baf` each `(n_obs, n_spots)`, a `SpatioGenomicCounts`, dense array or CSR.
    """

    ndim = 3

    def __init__(self, rdr: Any, baf: Any) -> None:
        self.rdr, self.baf = SpatioGenomicCounts.of(rdr), SpatioGenomicCounts.of(baf)

    @classmethod
    def of(cls, x: Any) -> SpatioGenomicPairCounts:
        """Return `x` as `SpatioGenomicPairCounts`: itself, or the dense `(n_obs, 2, n_spots)` a capture records, kept
        dense.

        Args: `x` a `SpatioGenomicPairCounts` or `(n_obs, 2, n_spots)` array.
        Returns: a `SpatioGenomicPairCounts`.
        """
        if isinstance(x, cls):
            return x
        x = np.asarray(x)
        return cls(x[:, 0, :], x[:, 1, :])

    @property
    def shape(self) -> tuple[int, int, int]:
        return (self.rdr.shape[0], 2, self.rdr.shape[1])

    def spot_sums(self, idx: Any = None) -> np.ndarray:
        """Compute `np.sum(single_X[:, :, idx], axis=-1)`.

        Args: `idx` spot indices or None.
        Returns: the `(n_obs, 2)` sums.
        """
        return np.stack([self.rdr.spot_sums(idx), self.baf.spot_sums(idx)], axis=1)


NPZ_CHANNELS = ("rdr", "baf")
"""`SpatioGenomicPairCounts`' channels, each written under `<name>__<channel>` in turn: read depth, then the BAF
numerator."""


def from_npz(loaded: Any, name: str) -> Any:
    """The entry `name` of a `binned_data.npz`, as the container the pages take.

    Reads this repository's sparse entries (`npz_entries`) and the dense entry CalicoST writes, and a file written
    before the layout, alike: a `(n_obs, n_spots)` array comes back a `SpatioGenomicCounts` and `single_X` a
    `SpatioGenomicPairCounts`, sparse where the file holds it sparse.

    Args: `loaded` an open `.npz` (or any mapping of its entries); `name` the entry.
    Returns: the `SpatioGenomicPairCounts` or `SpatioGenomicCounts`.
    Raises: `KeyError` where the file carries the entry under neither layout.
    """
    keys = set(getattr(loaded, "files", None) or loaded.keys())
    if f"{name}__{NPZ_CHANNELS[0]}__data" in keys or f"{name}__{NPZ_CHANNELS[0]}" in keys:
        return SpatioGenomicPairCounts(*(from_npz(loaded, f"{name}__{c}") for c in NPZ_CHANNELS))
    if f"{name}__data" in keys:
        shape = tuple(int(n) for n in np.asarray(loaded[f"{name}__shape"]))
        csr = sp.csr_matrix(
            (np.asarray(loaded[f"{name}__data"]), np.asarray(loaded[f"{name}__indices"]),
             np.asarray(loaded[f"{name}__indptr"])), shape=shape)  # fmt: skip
        return SpatioGenomicCounts(csr)
    if name not in keys:
        raise KeyError(f"{name} is in neither layout of the file: {sorted(keys)}")
    values = np.asarray(loaded[name])
    return SpatioGenomicPairCounts.of(values) if values.ndim == 3 else SpatioGenomicCounts(values)
