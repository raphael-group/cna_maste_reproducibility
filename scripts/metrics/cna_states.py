"""A run's clone profile as CNA states and fields, per segment or gene, for every method alike.

- `get_cna_states`, `get_uncalled` -- the six CNA states, and the mark of no call;
- `get_profile_fields` -- the fields each method's profile gives, `ascn`, `Uest`, `tcn`, `state` and `tcn_state`;
- `get_ab_states` -- phase-free AB pairs as CNA states;
- `get_called_states` -- a run's CNA state per segment and clone;
- `get_field_table` -- a run's fields per segment and clone;
- `get_comparison_field`, `get_phase_free_pairs`, `get_field_order` -- the field two runs are compared at, AB pairs
  phase-free, and a field's values in table order;
- `get_profile_copies` -- a run's A, B and total copies per segment and clone, where its profile gives copies.

A CNA state is one of Numbat's six, Numbat's `loh` written `cnloh`, by an AB pair's balance, its alleles' copies
equal, and its total copies:

    |            | total 2     | above 2        | below 2    |
    | balanced   | neu (1:1)   | bamp (2:2)     | bdel (0:0) |
    | unbalanced | cnloh (2:0) | amp (2:1, 3:1) | del (1:0)  |

each pair as Numbat's file descriptions name it (https://kharchenkolab.github.io/numbat/articles/descriptions.html):
"`l11_x`: ... CNV state 1:1 (neutral)", "`l10_x`: ... CNV state 1:0 (deletion)", "`l20_x`: ... CNV state 2:0
(CNLoH)", "`l00_x`: ... CNV state 0:0 (biallelic deletion)", "`l22_x`: ... CNV state 2:2 (biallelic amplification)",
"`p_amp`: ... amplification states (2:1, 3:1)".

Each method's fields per clone:

- an AB-pair profile (cnamaste, CalicoST, a simulation's truth): `ascn` `A|B` as the table phases it, `tcn` `A + B`,
  `state` of the phase-free pair;
- Numbat, states: `state`; `ascn` and `tcn` the state's one pair where it has one (`NUMBAT_PAIRS`), `unknown` for
  `amp` and `bamp`;
- Clonalscope, total copies: `Uest`, its cluster's mean coverage fold change per region (`AssignCluster`'s "Uest":
  "Average fold changes of all regions for each identify subclone/cluster", Clonalscope 07dd0c5 `R/AssignCluster.R`),
  `tcn`, its total copies, and `state` (`read_method_input.read_clonalscope`);
- every method: `tcn_state`, its state of total copies alone, `del` (< 2), `neu` (2) or `amp` (> 2), as Clonalscope's
  `state` and copytyping-benchmark's `analysis/scripts/cn_transform.py` (1204da9) read total copies: "DEL (-1) = total
  CN < 2, NEU (0) = total CN 2, AMP (1) = total CN > 2". It follows from `state`, each state's total copies below,
  at or above 2 (`TCN_STATE_OF`).

A field is `uncalled` where the run's profile has no value for the clone.

Two runs are compared at the finest field both give (`get_comparison_field`): their AB pairs, phase-free, where both
give pairs; `tcn_state` where either gives total copies alone; else `state`.
"""

from __future__ import annotations

import re
from typing import NamedTuple

import numpy as np
import pandas as pd

CNA_STATES = ("neu", "del", "amp", "cnloh", "bamp", "bdel")
"""The CNA states: Numbat's `cnv_state` values, its `loh` written `cnloh`."""

UNCALLED = "uncalled"
"""A field's value where the run calls nothing."""

UNKNOWN = "unknown"
"""`ascn` and `tcn` of a Numbat state that is not one pair."""

TCN_STATES = ("neu", "del", "amp")
"""The CNA states of total copies alone: 2, below 2 and above 2."""

TCN_STATE_OF = {"neu": "neu", "cnloh": "neu", "del": "del", "bdel": "del", "amp": "amp", "bamp": "amp"}
"""Each CNA state's state of total copies alone: its total copies at 2 (`neu`, `cnloh`), below (`del`, `bdel`) or
above (`amp`, `bamp`)."""

NUMBAT_PAIRS = {"neu": (1, 1), "del": (1, 0), "cnloh": (2, 0), "bdel": (0, 0)}
"""The Numbat states of one AB pair, its `l11_x`, `l10_x`, `l20_x` and `l00_x` states; `amp` is 2:1 or 3:1, and
`bamp` is left unknown."""


class ProfileCopies(NamedTuple):
    """A run's copies per segment and clone (`get_profile_copies`): `clones` the clone numbers, `a`, `b` and `total`
    each `(segments, clones)`, NaN where the run calls nothing or gives no allele copies."""

    clones: list[int]
    a: np.ndarray
    b: np.ndarray
    total: np.ndarray


def get_cna_states() -> tuple[str, ...]:
    """Get the CNA states, `neu`, `del`, `amp`, `cnloh`, `bamp` and `bdel`, in table order."""
    return CNA_STATES


def get_uncalled() -> str:
    """Get the mark of a field the run calls nothing on, `uncalled`."""
    return UNCALLED


def get_clone_count(profile: pd.DataFrame) -> int:
    """Get the clones of a `MethodInput.profile`: its `clone<k> A`, else its `clone<k>` state columns.

    Args:
        profile: the profile.

    Returns:
        The number of clones.
    """
    pattern = r"clone\d+ A" if "clone0 A" in profile else r"clone\d+"
    return sum(1 for column in profile.columns if re.fullmatch(pattern, str(column)))


def get_profile_fields(profile: pd.DataFrame) -> tuple[str, ...]:
    """Get the fields a profile gives, in the gene-level table's order.

    Args:
        profile: a `MethodInput.profile`: AB pairs, Clonalscope's states with `clone<k> Uest`, or Numbat's states.

    Returns:
        `ascn`, `tcn`, `state`, `tcn_state` for AB pairs and Numbat; `Uest`, `tcn`, `state`, `tcn_state` for
        Clonalscope.
    """
    if "clone0 Uest" in profile:
        return ("Uest", "tcn", "state", "tcn_state")
    return ("ascn", "tcn", "state", "tcn_state")


def get_ab_states(codes: np.ndarray) -> np.ndarray:
    """Get phase-free AB pairs' CNA states, by their balance and total copies.

    Args:
        codes: phase-free codes, `minor * 1_000 + major` (`supervised_metrics.phase_free`), any shape.

    Returns:
        The state of each, one of `CNA_STATES`, in `codes`' shape.
    """
    codes = np.asarray(codes)
    minor, major = codes // 1_000, codes % 1_000
    total = minor + major
    balanced = np.where(total > 2, "bamp", np.where(total < 2, "bdel", "neu"))
    unbalanced = np.where(total > 2, "amp", np.where(total < 2, "del", "cnloh"))
    return np.where(minor == major, balanced, unbalanced)


def get_called_states(profile: pd.DataFrame) -> np.ndarray:
    """Get a run's CNA state per segment and clone, from its profile's AB pairs or its states.

    Args:
        profile: a `MethodInput.profile`, per clone `k` its `clone<k> A` and `clone<k> B`, or its state `clone<k>`,
            Numbat's `loh` read as `cnloh`.

    Returns:
        `(segments, clones)` states, "" where the run calls none.

    Raises:
        ValueError: where a state is not one of `CNA_STATES`.
    """
    n = get_clone_count(profile)
    if "clone0 A" in profile:
        a = profile[[f"clone{k} A" for k in range(n)]].to_numpy(dtype=float)
        b = profile[[f"clone{k} B" for k in range(n)]].to_numpy(dtype=float)
        called = ~(np.isnan(a) | np.isnan(b))
        codes = np.where(called, np.minimum(a, b) * 1_000 + np.maximum(a, b), 0).astype(np.int64)
        return np.where(called, get_ab_states(codes), "")
    states = profile[[f"clone{k}" for k in range(n)]].fillna("").to_numpy(dtype=str)
    # NB Numbat's `loh` is its CNLoH state (R/vis.R `cnv_labels`).
    states = np.where(states == "loh", "cnloh", states)
    unknown = set(states.ravel()) - {*CNA_STATES, ""}
    if unknown:
        raise ValueError(f"states {sorted(unknown)} not in CNA_STATES")
    return states


def get_field_table(profile: pd.DataFrame) -> pd.DataFrame:
    """Get a run's fields per segment and clone, `uncalled` where it calls nothing.

    Args:
        profile: a `MethodInput.profile`.

    Returns:
        `clone<k> <field>` per clone and `get_profile_fields(profile)` field, a row per profile row, its index.
    """
    states = get_called_states(profile)
    fields = get_profile_fields(profile)
    columns: dict[str, np.ndarray] = {}
    for k in range(states.shape[1]):
        found = {
            "state": states[:, k].astype(object),
            "tcn_state": np.array([TCN_STATE_OF.get(state, "") for state in states[:, k]], dtype=object),
        }
        if "clone0 A" in profile:
            a = profile[f"clone{k} A"].to_numpy(dtype=float)
            b = profile[f"clone{k} B"].to_numpy(dtype=float)
            found["ascn"] = np.array([f"{x:.0f}|{y:.0f}" for x, y in zip(a, b, strict=True)], dtype=object)
            found["tcn"] = np.array([f"{x + y:.0f}" for x, y in zip(a, b, strict=True)], dtype=object)
        elif "Uest" in fields:
            found["Uest"] = profile[f"clone{k} Uest"].to_numpy(dtype=object)
            found["tcn"] = np.array([f"{x:.0f}" for x in profile[f"clone{k} tcn"].to_numpy(dtype=float)], dtype=object)
        else:
            pairs = [NUMBAT_PAIRS.get(state) for state in states[:, k]]
            found["ascn"] = np.array([f"{p[0]}|{p[1]}" if p else UNKNOWN for p in pairs], dtype=object)
            found["tcn"] = np.array([f"{p[0] + p[1]}" if p else UNKNOWN for p in pairs], dtype=object)
        for field in fields:
            columns[f"clone{k} {field}"] = np.where(states[:, k] != "", found[field], UNCALLED)
    return pd.DataFrame(columns, index=profile.index)


def get_comparison_field(profile_a: pd.DataFrame, profile_b: pd.DataFrame) -> str:
    """Get the field two runs' calls are compared at, the finest both give.

    Args:
        profile_a: a `MethodInput.profile`.
        profile_b: another.

    Returns:
        `ascn` where both give AB pairs (cnamaste, CalicoST), `tcn_state` where either gives total copies alone
        (Clonalscope), else `state` (Numbat).
    """
    profiles = (profile_a, profile_b)
    if all("clone0 A" in profile for profile in profiles):
        return "ascn"
    if any("clone0 Uest" in profile for profile in profiles):
        return "tcn_state"
    return "state"


def get_phase_free_pairs(values: pd.DataFrame) -> pd.DataFrame:
    """Get `ascn` values phase-free, `A|B` written `(major,minor)`, as `palette.state_colour` keys a pair; other
    values, e.g. `uncalled`, as given.

    Args:
        values: `ascn` values, `A|B`.

    Returns:
        The values, `values`' shape.
    """

    def get_pair(value: str) -> str:
        a, sep, b = str(value).partition("|")
        if not sep:
            return str(value)
        return f"({max(int(a), int(b))},{min(int(a), int(b))})"

    return values.apply(lambda column: column.map(get_pair))


def get_field_order(field: str, values: list[pd.DataFrame]) -> list[str]:
    """Get a field's values in table order, `uncalled` left out.

    Args:
        field: `state`, `tcn_state` or `ascn` (`get_comparison_field`).
        values: the runs' values of the field, `ascn` phase-free (`get_phase_free_pairs`).

    Returns:
        `CNA_STATES` for `state`, `TCN_STATES` for `tcn_state`; for `ascn`, the pairs found, (1,1) first, then by
        total copies, then by the major allele's.

    Raises:
        ValueError: where a field is not one of the three.
    """
    if field == "state":
        return list(CNA_STATES)
    if field == "tcn_state":
        return list(TCN_STATES)
    if field != "ascn":
        raise ValueError(f"no order of field {field!r}")
    found = {str(v) for table in values for v in np.ravel(table.to_numpy())} - {UNCALLED}
    pairs = [tuple(int(x) for x in pair.strip("()").split(",")) for pair in found]
    return [f"({a},{b})" for a, b in sorted(pairs, key=lambda p: (p != (1, 1), p[0] + p[1], p[0]))]


def get_profile_copies(profile: pd.DataFrame) -> ProfileCopies | None:
    """Get a run's copies per segment and clone: an AB-pair profile's `A`, `B` and `A + B`; Clonalscope's total copies,
    `tcn`, its `A` and `B` NaN; None for a profile of states alone (Numbat).

    Args:
        profile: a `MethodInput.profile`.

    Returns:
        The clones, as numbered by its `clone<k> A`, else `clone<k> tcn`, columns, and their copies; or None.
    """
    field = "A" if any(re.fullmatch(r"clone\d+ A", str(c)) for c in profile.columns) else "tcn"
    clones = sorted(int(m.group(1)) for c in profile.columns if (m := re.fullmatch(rf"clone(\d+) {field}", str(c))))
    if not clones:
        return None
    if field == "A":
        a = profile[[f"clone{k} A" for k in clones]].to_numpy(dtype=float)
        b = profile[[f"clone{k} B" for k in clones]].to_numpy(dtype=float)
        return ProfileCopies(clones, a, b, a + b)
    total = profile[[f"clone{k} tcn" for k in clones]].to_numpy(dtype=float)
    return ProfileCopies(clones, np.full_like(total, np.nan), np.full_like(total, np.nan), total)
