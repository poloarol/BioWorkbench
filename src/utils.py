
import scanpy as sc
import squidpy as sq
import numpy as np

def load_data(filepath: str, 
            cell_columns: list[str] = None, 
            genes: list[str] = None,
            is_spatial: bool = False) -> dict[str, sc.AnnData]:
    """
    Load a dataset from a file and return a dictionary of AnnData objects.

    Parameters:
    - filepath: str
        The path to the file to load.
    - cell_columns: list[str], optional
        A list of column names to use as cell metadata. If None, all columns will be used.
    - genes: list[str], optional
        A list of column names to use as gene metadata. If None, all columns will be used.
    - is_spatial: bool, optional
        Whether the dataset is spatial. If True, the function will attempt to load spatial data.

    Returns:
    - data: dict[str, sc.AnnData]
        A dictionary containing the loaded data.
    """
    
    adata = None
    sdata = None
    
    if genes is None and cell_columns is None:
        adata = sc.read(filepath)
        _calculate_qc_metrics(adata, is_spatial=is_spatial)
        adata.layers['raw'] = adata.X.copy()
        
        if is_spatial:
            blank_genes = adata.var_names[adata.var_names.str.startswith("Blank-")].tolist()
            blank_counts = adata[:, blank_genes].X.sum(axis=1)
            adata.obs["blank_counts"] = blank_counts
            adata.obs["pct_counts_blank"] = (
                adata.obs["blank_counts"] / adata.obs["total_counts"]
            ) * 100
            adata.obsm["spatial"] = adata.obs[["center_x", "center_y"]].to_numpy()
    else:
        adata = sc.read(filepath, backed='r')
        sdata = adata[cell_columns, genes].copy()
        _calculate_qc_metrics(sdata, is_spatial=is_spatial)
        sdata.layers['raw'] = sdata.X.copy()
        if is_spatial:
            blank_genes = sdata.var_names[sdata.var_names.str.startswith("Blank-")].tolist()
            blank_counts = sdata[:, blank_genes].X.sum(axis=1)
            sdata.obs["blank_counts"] = blank_counts
            sdata.obs["pct_counts_blank"] = (
                sdata.obs["blank_counts"] / sdata.obs["total_counts"]
            ) * 100
            sdata.obsm["spatial"] = sdata.obs[["center_x", "center_y"]].to_numpy()
    
    return {'raw': adata, 'subset': sdata}


def _calculate_qc_metrics(adata, is_spatial=False):
    """
    Calculate QC metrics for an AnnData object.

    For spatial/MERFISH data, calculate the standard cell/spot
    metrics without top-N gene percentages.

    For single-cell data, calculate additional QC metrics based
    on available gene annotations.
    """

    if is_spatial:
        sc.pp.calculate_qc_metrics(
            adata,
            percent_top=None,
            inplace=True,
        )
        return

    # --------------------------------------------------------
    # Single-cell QC
    # --------------------------------------------------------

    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=None,
        inplace=True,
        log1p=True,
    )


def write_to_disk(adata: sc.AnnData, filepath: str) -> None:
    """
    Write an AnnData object to disk.

    Parameters:
    - adata: sc.AnnData
        The AnnData object to write to disk.
    - filepath: str
        The path to the file to write.

    Returns:
    - None
        The function writes the AnnData object to disk.
    """
    adata.write(filepath)


def calculate_module_score(
    adata: sc.AnnData,
    genes: list[str],
    score_name: str = "module_score",
):
    """
    Calculate a module score for a set of genes in an AnnData object.

    Parameters:
    - adata: sc.AnnData
        The AnnData object containing the gene expression data.
    - genes: list[str]
        A list of gene names to include in the module score.
    - score_name: str, optional (default="module_score")
        The name of the column in adata.obs to store the module score.

    Returns:
    - adata: sc.AnnData
        The AnnData object with the calculated module score added to adata.obs.
    """
    genes = [gene for gene in genes if gene in adata.var_names]

    if not genes:
        raise ValueError("None of the selected genes were found in the dataset.")

    sc.tl.score_genes(
        adata,
        gene_list=genes,
        score_name=score_name,
    )

    return adata