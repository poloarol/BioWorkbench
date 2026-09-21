
import pytest

from unittest.mock import patch
from src.filtering import (
    _filter_cells, _filter_genes, 
    _filter_mitochondrial_genes, _identify_doublets,
    data_filtering
)



## filter cells

def test_filtering_removes_low_quality_cells(small_adata):
    result = _filter_cells(small_adata, min_genes=3)

    assert result.n_obs < small_adata.n_obs

def test_filtering_keeps_good_cells(small_adata):
    result = _filter_cells(small_adata, min_genes=1)

    assert result.n_obs == small_adata.n_obs

def test_filtering_can_remove_all_cells(small_adata):
    result = _filter_cells(small_adata, min_genes=100)

    assert result.n_obs == 0

def test_filtering_does_not_modify_original(small_adata):
    original_n_obs = small_adata.n_obs

    _filter_cells(small_adata, min_genes=3)

    assert small_adata.n_obs == original_n_obs

## filter genes

# def test_filter_genes_default_threshold(small_adata):
#     result = _filter_genes(small_adata)

#     assert result.n_vars == 2
#     assert list(result.var_names) == ["GeneA", "GeneD"]

def test_filter_genes_min_cells_2(small_adata):
    result = _filter_genes(small_adata, min_cells=2)

    assert result.n_vars == 5


def test_filter_genes_strict_threshold(small_adata):
    result = _filter_genes(small_adata, min_cells=4)

    assert result.n_vars == 0


def test_filter_genes_does_not_modify_original(small_adata):
    original_n_vars = small_adata.n_vars
    original_genes = list(small_adata.var_names)

    result = _filter_genes(small_adata, min_cells=3)

    assert result is not small_adata
    assert small_adata.n_vars == original_n_vars
    assert list(small_adata.var_names) == original_genes


## Filter Mitochondrial Genes

def test_filter_mitochondrial_genes_default_threshold(small_adata):
    result = _filter_mitochondrial_genes(small_adata)

    assert list(result.obs_names) == ["Cell1", "Cell2", "Cell4"]


def test_filter_mitochondrial_genes_custom_threshold(small_adata):
    result = _filter_mitochondrial_genes(
        small_adata,
        max_mt_percentage=10
    )

    assert list(result.obs_names) == ["Cell1"]


def test_filter_mitochondrial_genes_removes_threshold_value(small_adata):
    result = _filter_mitochondrial_genes(
        small_adata,
        max_mt_percentage=10
    )

    assert "Cell2" not in result.obs_names


def test_filter_mitochondrial_genes_does_not_modify_original(small_adata):
    original_n_obs = small_adata.n_obs

    result = _filter_mitochondrial_genes(
        small_adata,
        max_mt_percentage=20
    )

    assert result is not small_adata
    assert small_adata.n_obs == original_n_obs


def test_filter_mitochondrial_genes_requires_mt_percentage(small_adata):
    del small_adata.obs["pct_counts_mt"]

    with pytest.raises(KeyError):
        _filter_mitochondrial_genes(small_adata)


## Identifying Doublets

# def test_identify_doublets_adds_results(doublet_adata):
#     result = _identify_doublets(doublet_adata)

#     assert "doublet_score" in result.obs
#     assert "predicted_doublet" in result.obs



def test_identify_doublets_passes_expected_rate(small_adata):
    with patch("src.filtering.sc.pp.scrublet") as mock_scrublet:
        _identify_doublets(
            small_adata,
            expected_doublet_rate=0.10
        )

        mock_scrublet.assert_called_once_with(
            small_adata,
            expected_doublet_rate=0.10
        )


def test_identify_doublets_uses_default_rate(small_adata):
    with patch("src.filtering.sc.pp.scrublet") as mock_scrublet:
        _identify_doublets(small_adata)

        mock_scrublet.assert_called_once_with(
            small_adata,
            expected_doublet_rate=0.05
        )


def test_identify_doublets_returns_same_adata(small_adata):
    with patch("src.filtering.sc.pp.scrublet"):
        result = _identify_doublets(small_adata)

    assert result is small_adata


# Data Filtering Aggregator


def test_data_filtering_non_spatial(small_adata):
    with patch("src.filtering._filter_cells") as mock_cells, \
         patch("src.filtering._filter_genes") as mock_genes, \
         patch("src.filtering._filter_mitochondrial_genes") as mock_mt, \
         patch("src.filtering._identify_doublets") as mock_doublets:

        mock_cells.return_value = small_adata
        mock_genes.return_value = small_adata
        mock_mt.return_value = small_adata
        mock_doublets.return_value = small_adata

        result = data_filtering(
            small_adata,
            min_genes=100,
            min_cells=5,
            max_mt_percentage=15,
            expected_doublet_rate=0.1,
            is_spatial=False,
        )

    mock_cells.assert_called_once_with(small_adata, 100)
    mock_genes.assert_called_once_with(small_adata, 5)
    mock_mt.assert_called_once_with(small_adata, 15)
    mock_doublets.assert_called_once_with(small_adata, 0.1)

    assert result is small_adata


def test_data_filtering_spatial(small_adata):
    with patch("src.filtering._filter_cells") as mock_cells, \
         patch("src.filtering._filter_genes") as mock_genes, \
         patch("src.filtering._identify_blanks") as mock_blanks, \
         patch("src.filtering._filter_mitochondrial_genes") as mock_mt, \
         patch("src.filtering._identify_doublets") as mock_doublets:

        mock_cells.return_value = small_adata
        mock_genes.return_value = small_adata
        mock_blanks.return_value = small_adata

        result = data_filtering(
            small_adata,
            min_genes=100,
            min_cells=5,
            blank_max_percentage=8,
            is_spatial=True,
        )

    mock_cells.assert_called_once_with(small_adata, 100)
    mock_genes.assert_called_once_with(small_adata, 5)

    mock_blanks.assert_called_once_with(
        small_adata,
        max_blank_percentage=8
    )

    mock_mt.assert_not_called()
    mock_doublets.assert_not_called()

    assert result is small_adata


def test_data_filtering_defaults(small_adata):
    with patch("src.filtering._filter_cells") as mock_cells, \
         patch("src.filtering._filter_genes") as mock_genes, \
         patch("src.filtering._filter_mitochondrial_genes") as mock_mt, \
         patch("src.filtering._identify_doublets") as mock_doublets:

        mock_cells.return_value = small_adata
        mock_genes.return_value = small_adata
        mock_mt.return_value = small_adata
        mock_doublets.return_value = small_adata

        data_filtering(small_adata)

    mock_cells.assert_called_once_with(small_adata, 200)
    mock_genes.assert_called_once_with(small_adata, 3)
    mock_mt.assert_called_once_with(small_adata, 20)
    mock_doublets.assert_called_once_with(small_adata, 0.05)


def test_data_filtering_passes_results_between_steps(small_adata):
    filtered_cells = small_adata.copy()
    filtered_genes = small_adata.copy()
    filtered_mt = small_adata.copy()
    filtered_doublets = small_adata.copy()

    with patch("src.filtering._filter_cells", return_value=filtered_cells) as mock_cells, \
         patch("src.filtering._filter_genes", return_value=filtered_genes) as mock_genes, \
         patch("src.filtering._filter_mitochondrial_genes", return_value=filtered_mt) as mock_mt, \
         patch("src.filtering._identify_doublets", return_value=filtered_doublets) as mock_doublets:

        result = data_filtering(small_adata)

    mock_cells.assert_called_once_with(small_adata, 200)
    mock_genes.assert_called_once_with(filtered_cells, 3)
    mock_mt.assert_called_once_with(filtered_genes, 20)
    mock_doublets.assert_called_once_with(filtered_mt, 0.05)

    assert result is filtered_doublets
