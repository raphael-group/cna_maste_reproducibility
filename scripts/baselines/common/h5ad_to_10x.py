#!/usr/bin/env python
"""Write a CalicoST simulation's filtered_feature_bc_matrix.h5ad as a 10x matrix directory.

Runpeng Luo (2026-10-30)

Dependencies: anndata, numpy, pandas, scipy (calicost_env).

Usage:
    python -I common/h5ad_to_10x.py <h5ad> <ordered_genes_tsv> <out_dir>
            h5ad               # simulated_data_related/<sim>/filtered_feature_bc_matrix.h5ad
            ordered_genes_tsv  # Cell Ranger GRCh38-2020-A genes/ordered_genes.tsv (Gene, Ensembl_id, ...)
            out_dir            # output filtered_feature_bc_matrix/

Inputs:
    h5ad: spot x gene integer UMI counts; var_names are gene symbols only.
    ordered_genes_tsv: symbol -> Ensembl ID table of the reference the simulation was built on.
Outputs:
    out_dir/matrix.mtx.gz: gene x spot counts (MatrixMarket integer).
    out_dir/barcodes.tsv.gz: spot barcodes (obs_names).
    out_dir/features.tsv.gz: Ensembl ID, gene symbol, "Gene Expression" (Space Ranger layout).
Notes/References:
    Clonalscope's EstRegionCov joins features column 1 to the gene BED's geneID (Ensembl), and
    numbat/build_count_mat.R reads gene symbols from column 2, so both columns are needed.
    Symbols that are absent from ordered_genes_tsv or occur more than once there keep the symbol
    as column 1, which matches no Ensembl-keyed BED row.
    Simulation deposit: https://doi.org/10.5281/zenodo.14175627
"""

import gzip
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def main(h5ad: str, ordered_genes: str, out_dir: str) -> None:
    """Convert one simulation.

    Args:
        h5ad: Path to the simulated AnnData file.
        ordered_genes: Path to ordered_genes.tsv.
        out_dir: Output 10x matrix directory.
    """
    adata = ad.read_h5ad(h5ad)
    symbols = pd.Series(adata.var_names.astype(str))
    ref = pd.read_csv(ordered_genes, sep="\t")
    ref = ref[~ref["Gene"].duplicated(keep=False)]
    ids = symbols.map(dict(zip(ref["Gene"], ref["Ensembl_id"])))
    print(f"{h5ad}: {adata.n_obs} spots x {adata.n_vars} genes; {int(ids.notna().sum())} mapped to Ensembl")
    ids = ids.fillna(symbols)

    x = adata.X if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    x = sp.csc_matrix(x.T).astype(np.int64)
    os.makedirs(out_dir, exist_ok=True)
    with gzip.open(os.path.join(out_dir, "matrix.mtx.gz"), "wb") as f:
        scipy.io.mmwrite(f, x, field="integer")
    pd.Series(adata.obs_names).to_csv(os.path.join(out_dir, "barcodes.tsv.gz"), sep="\t", header=False, index=False)
    pd.DataFrame({"id": ids, "name": symbols, "type": "Gene Expression"}).to_csv(
        os.path.join(out_dir, "features.tsv.gz"), sep="\t", header=False, index=False
    )


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("usage: common/h5ad_to_10x.py <h5ad> <ordered_genes_tsv> <out_dir>")
    main(*sys.argv[1:])
