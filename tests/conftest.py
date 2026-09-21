# tests/conftest.py

import numpy as np
import pandas as pd
import pytest
import scanpy as sc


np.random.seed(42) 

import numpy as np
import anndata as ad
import pytest


@pytest.fixture
def small_adata():
    X = np.array([
        [10, 0, 0, 5, 0],
        [8,  0, 1, 4, 0],
        [0,  9, 8, 0, 1],
        [0, 10, 7, 0, 0],
        [1,  0, 0, 1, 10],
    ])

    adata = ad.AnnData(X)

    adata.var_names = [
        "GeneA",
        "GeneB",
        "GeneC",
        "GeneD",
        "MT-GeneE",
    ]

    adata.obs_names = [
        "Cell1",
        "Cell2",
        "Cell3",
        "Cell4",
        "Cell5",
    ]
    
    adata.obs["pct_counts_mt"] = [5.0, 10.0, 25.0, 15.0, 30.0]

    return adata


@pytest.fixture
def doublet_adata():
    np.random.seed(42)

    X = np.random.poisson(
        lam=5,
        size=(100, 50)
    )

    adata = sc.AnnData(X)

    adata.obs_names = [
        f"Cell{i}" for i in range(100)
    ]

    adata.var_names = [
        f"Gene{i}" for i in range(50)
    ]

    return adata
