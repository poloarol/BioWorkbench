# BioWorkbench

**An interactive computational workbench for single-cell and spatial genomics.**

BioWorkbench is a lightweight Streamlit application for exploring, analyzing, and visualizing single-cell and spatial genomics datasets.

The project is designed to provide an interactive layer on top of reproducible computational workflows — allowing researchers to explore quality control metrics, dimensionality reduction, clustering, cell-type annotations, and spatial organization without having to repeatedly modify analysis code.

> **Status:** Active development

## Overview

Modern single-cell and spatial genomics analyses generate complex datasets that often require substantial computational expertise to explore.

BioWorkbench provides a structured interface for common analysis and visualization tasks while keeping the underlying data in the standard **AnnData** format.

The application currently focuses on:

* Single-cell RNA-seq analysis
* Cell-level quality control and filtering
* Dimensionality reduction and clustering
* Cell-type annotation
* Interactive exploratory visualization
* Cell-resolved spatial transcriptomics
* Platform-aware spatial quality control

The spatial component is being developed around technologies including **MERFISH, Xenium, and CosMx**, with **Visium HD** planned as a separate high-resolution spatial workflow.

## Current Features

### Single-cell analysis

* AnnData (`.h5ad`) input
* Cell and gene filtering
* Mitochondrial QC
* Doublet detection with Scrublet
* Normalization
* Highly variable gene selection
* PCA
* Neighborhood graph construction
* UMAP
* Leiden clustering
* Marker gene analysis
* CellTypist-based cell-type annotation
* Custom annotation support
* Export of processed AnnData objects

### Spatial analysis

BioWorkbench is being extended beyond conventional single-cell workflows to support cell-resolved spatial transcriptomics.

Current development includes:

* Technology selection for spatial datasets
* Platform-specific filtering parameters
* Spatial QC visualization
* Cell-level spatial coordinates
* Visualization of categorical and continuous spatial features
* MERFISH-oriented QC metrics
* Spatial visualization based on standardized `AnnData.obsm["spatial"]` coordinates

The spatial architecture is designed to accommodate differences between technologies rather than treating all spatial datasets as interchangeable.

For example, imaging-based platforms such as MERFISH, Xenium, and CosMx may provide cell morphology, segmentation, blank/negative probe information, and transcript-level QC metrics that differ substantially from conventional scRNA-seq QC.

## Analysis Workflow

The current single-cell workflow follows a modular analysis structure:

```text
AnnData
   │
   ▼
Quality Control & Filtering
   │
   ▼
Normalization
   │
   ▼
Highly Variable Genes
   │
   ▼
PCA
   │
   ▼
Neighbourhood Graph
   │
   ▼
UMAP / Leiden Clustering
   │
   ▼
Marker Genes
   │
   ▼
Cell-type Annotation
   │
   ▼
Processed AnnData
```

For spatial datasets, the workflow is being extended with platform-aware QC and spatial visualization:

```text
Spatial AnnData
      │
      ▼
Technology-specific QC
      │
      ├── Cell-level metrics
      ├── Transcript / count metrics
      ├── Blank / negative probe metrics
      └── Morphology / segmentation metrics
      │
      ▼
Spatial Visualization
      │
      ▼
Spatial Analysis
```

## Spatial Data Model

BioWorkbench uses `AnnData` as the common data structure while recognizing that different spatial technologies represent biological space differently.

### Cell-resolved spatial technologies

For technologies such as:

* MERFISH
* Xenium
* CosMx

the primary analytical unit can be represented as:

```text
cell × gene
```

with spatial coordinates stored in:

```python
adata.obsm["spatial"]
```

and cell-level metadata stored in:

```python
adata.obs
```

This allows BioWorkbench to use a common visualization layer while retaining technology-specific QC and metadata.

### High-resolution spatial transcriptomics

Visium HD is being treated as a distinct workflow because its analytical units are high-resolution spatial bins rather than inherently segmented cells.

Planned support will therefore address:

* Bin-level QC
* Tissue image visualization
* Spatial clustering
* Spatially variable genes
* Tissue domains
* Cell-type composition / deconvolution

## Configuration

Analysis parameters are being separated from the Streamlit interface using YAML configuration.

For example:

```text
config/
└── params.yaml
```

This allows filtering parameters to vary according to dataset type and spatial technology rather than applying conventional scRNA-seq thresholds indiscriminately.

Conceptually:

```yaml
single_cell:
  ...

spatial:
  MERFISH:
    ...

  Xenium:
    ...

  CosMX:
    ...

  Visium HD:
    ...
```

This design is intended to make platform-specific analysis behavior explicit, reproducible, and easier to extend.

## Project Structure

```text
BioWorkbench/
│
├── app.py
│
├── pages/
│   ├── 01_filtering.py
│   ├── 02_clustering.py
│   └── 03_annotation.py
│
├── src/
│   ├── filtering.py
│   ├── clustering.py
│   ├── plotting.py
│   └── utils.py
│
├── config/
│   └── params.yaml
│
├── notebooks/
│
└── README.md
```

The application separates the Streamlit interface from computational and visualization logic wherever practical.

## Technology Stack

BioWorkbench is built primarily with Python and the scientific Python ecosystem.

* **Python**
* **Streamlit** — interactive application interface
* **AnnData** — data structure
* **Scanpy** — single-cell analysis
* **Squidpy** — spatial analysis and visualization
* **CellTypist** — cell-type annotation
* **Scrublet** — doublet detection
* **Pandas** — data manipulation
* **NumPy** — numerical computing
* **Matplotlib / Seaborn / Plotly** — visualization

## Design Philosophy

BioWorkbench is not intended to replace specialized analysis pipelines or computational workflows.

Instead, it provides an **interactive analytical workbench** on top of those workflows.

The goal is to separate:

```text
Reproducible computational analysis
              │
              ▼
        Validated results
              │
              ▼
        BioWorkbench
              │
              ▼
 Interactive exploration & visualization
```

This makes the application useful both as a research tool and as a potential interface for delivering computational analyses to collaborators and clients.

## Roadmap

### Single-cell

* [x] AnnData input
* [x] QC and filtering
* [x] Doublet detection
* [x] PCA
* [x] UMAP
* [x] Leiden clustering
* [x] Marker gene analysis
* [x] CellTypist annotation
* [x] Custom annotations
* [ ] Expanded differential expression
* [ ] Improved analysis provenance
* [ ] Project-level configuration

### Spatial transcriptomics

#### Cell-resolved

* [x] Spatial technology selection
* [x] Platform-aware filtering configuration
* [x] Spatial coordinate handling
* [x] Spatial QC visualization
* [x] MERFISH-oriented workflows
* [ ] Xenium workflows
* [ ] CosMx workflows
* [ ] Spatial domains
* [ ] Spatially variable genes

#### High-resolution

* [ ] Visium HD
* [ ] Tissue image integration
* [ ] Bin-level QC
* [ ] Spatial clustering
* [ ] Spatially variable genes
* [ ] Tissue domains
* [ ] Cell-type composition / deconvolution

## Development

Clone the repository and create a Python environment:

```bash
git clone https://github.com/poloarol/BioWorkbench.git
cd BioWorkbench

python -m venv env
```

Activate the environment and install the project dependencies.

Then launch the application with:

```bash
streamlit run app.py
```

## Project Status

BioWorkbench is an evolving project.

The single-cell workflow provides the foundation, while the current development effort is expanding the application toward **cell-resolved and high-resolution spatial genomics**.

The architecture is intentionally being developed around modular analysis components, platform-aware configuration, and standardized biological data structures so that additional spatial technologies can be incorporated without duplicating the entire analysis framework.

## Author

**Paul A. Wambo**

Computational scientist working at the intersection of computational biology, machine learning, and scientific software.

GitHub: [@poloarol](https://github.com/poloarol)

---

*BioWorkbench is developed as an open computational tool for exploratory analysis and visualization of single-cell and spatial genomics data.*
