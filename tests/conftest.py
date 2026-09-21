# tests/conftest.py

import numpy as np
import pandas as pd
import pytest
import scanpy as sc


np.random.seed(42) 

@pytest.fixture
def single_cell_adata():
    """Small synthetic single-cell AnnData object for testing."""


    n_cells = 20
    n_genes = 10

    X = np.random.poisson(
        lam=5,
        size=(n_cells, n_genes),
    )

    adata = sc.AnnData(X)

    adata.obs_names = [
        f"cell_{i:03d}"
        for i in range(n_cells)
    ]

    adata.var_names = [
        f"gene_{i:03d}"
        for i in range(n_genes)
    ]

    adata.obs["n_genes_by_counts"] = (
        np.count_nonzero(X, axis=1)
    )

    adata.obs["total_counts"] = X.sum(axis=1)

    adata.obs["pct_counts_mt"] = np.random.uniform(
        0,
        10,
        n_cells,
    )

    return adata