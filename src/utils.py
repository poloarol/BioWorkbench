
import scanpy as sc


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
        calculate_qc_metrics(adata)
        adata.layers['raw'] = adata.X.copy()
    else:
        adata = sc.read(filepath, backed='r')
        sdata = adata[cell_columns, genes].copy()
        calculate_qc_metrics(sdata)
        sdata.layers['raw'] = sdata.X.copy()
    
    if is_spatial:
        ...
    
    return {'raw': adata, 'subset': sdata}


def calculate_qc_metrics(adata: sc.AnnData, metrics: list[str] = ['mt', 'ribo']) -> None:
    """
    Calculate quality control metrics for an AnnData object.

    Parameters:
    - adata: sc.AnnData
        The AnnData object to calculate QC metrics for.
    - metrics: list[str]
        A list of metrics to calculate. Default is ['mt', 'ribo'].

    Returns:
    - None
        The function modifies the AnnData object in place.
    """
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=metrics, inplace=True, log1p=True)


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