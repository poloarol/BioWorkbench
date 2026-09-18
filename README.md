# BioWorkbench

**An interactive computational workbench for single-cell and spatial genomics.**

BioWorkbench is a lightweight Streamlit application for exploring, analyzing, and visualizing single-cell and spatial genomics datasets.

The project provides an interactive layer on top of reproducible computational workflows, allowing researchers to explore quality control metrics, dimensionality reduction, clustering, cell-type annotations, and spatial organization without repeatedly modifying analysis code.

> **Status:** Active development

## Overview

Modern single-cell and spatial genomics analyses generate complex datasets that often require substantial computational expertise to explore.

BioWorkbench provides a structured interface for common analysis and visualization tasks while keeping the underlying data in the standard **AnnData** format.

The project currently focuses on:

* Single-cell RNA-seq analysis
* Cell and gene quality control
* Filtering and preprocessing
* Dimensionality reduction
* Clustering
* Marker gene analysis
* Cell-type annotation
* Cell-resolved spatial transcriptomics
* Platform-aware spatial quality control
* Exploratory visualization

Spatial development is designed to accommodate multiple technologies, including **MERFISH, Xenium, CosMX, and Visium HD**.

## Current Features

### Single-cell analysis

BioWorkbench currently provides a modular single-cell analysis workflow including:

* AnnData (`.h5ad`) input
* Cell and gene filtering
* Mitochondrial quality control
* Doublet detection with Scrublet
* Normalization
* Highly variable gene selection
* PCA
* Neighborhood graph construction
* UMAP
* Leiden clustering
* Marker gene analysis
* CellTypist-based cell-type annotation
* Custom cell-type annotations
* Export of processed AnnData objects

### Exploratory visualization

Reusable plotting functions are provided for exploring analytical results.

Current visualization functionality includes:

* PCA visualization
* UMAP visualization
* Categorical variables
* Continuous variables
* Automatic categorical color handling
* Continuous color scales
* Legends for categorical variables
* Colorbars for continuous variables
* Spatial visualization using `AnnData.obsm["spatial"]`
* Visualization of multiple spatial variables

Categorical variables are displayed as discrete groups, while continuous variables are represented using continuous color scales.

### Spatial analysis

BioWorkbench is being extended beyond conventional single-cell workflows to support **cell-resolved spatial transcriptomics**.

The current architecture includes:

* Spatial technology selection
* Technology-specific configuration
* Platform-aware filtering
* Cell-level spatial coordinates
* Spatial quality-control visualization
* Categorical spatial visualization
* Continuous spatial visualization
* MERFISH-oriented QC

Spatial coordinates are standardized through:

```python
adata.obsm["spatial"]
```

The spatial architecture is deliberately platform-aware rather than assuming that all spatial datasets follow the same QC model as scRNA-seq.

For example, imaging-based platforms such as MERFISH, Xenium, and CosMX can provide blank/negative probe measurements, segmentation information, morphology measurements, and other QC metrics that differ from conventional single-cell sequencing assays.

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
Neighborhood Graph
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
Visualization
   │
   ▼
Processed AnnData
```

For spatial datasets, the workflow incorporates technology-specific QC and spatial visualization:

```text
Spatial AnnData
      │
      ▼
Technology-specific QC
      │
      ├── Cell-level metrics
      ├── Count metrics
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

BioWorkbench uses **AnnData** as the common data structure while recognizing that different spatial technologies represent biological space differently.

### Cell-resolved spatial technologies

For technologies such as:

* MERFISH
* Xenium
* CosMX

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

This provides a common interface for spatial visualization while retaining technology-specific metadata and QC metrics.

### High-resolution spatial transcriptomics

**Visium HD** is treated as a distinct workflow because its analytical units are high-resolution spatial bins rather than inherently segmented cells.

Planned support includes:

* Bin-level QC
* Tissue image visualization
* Spatial clustering
* Spatially variable genes
* Tissue domains
* Cell-type composition / deconvolution

## Platform-Aware Configuration

Analysis parameters are separated from the Streamlit interface using YAML configuration.

```text
config/
└── params.yaml
```

The current configuration distinguishes between single-cell and spatial workflows:

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

This allows filtering behavior to reflect the characteristics of the underlying assay rather than applying conventional scRNA-seq thresholds indiscriminately.

For example, the current spatial configuration includes parameters for:

* Minimum genes
* Minimum cells
* Blank-probe thresholds
* Platform-specific QC behavior

The configuration architecture is intended to make analytical assumptions explicit, reproducible, and easier to extend as additional spatial technologies are implemented.

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
|   |__ 04_visualization.py
|   |__ 05_download.py
│
├── src/
│   ├── filtering.py
│   ├── clustering.py
│   ├── plotting.py
│   └── utils.py
|   |__ visualizer.py
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
* **AnnData** — biological data structure
* **Scanpy** — single-cell analysis
* **Squidpy** — spatial analysis
* **CellTypist** — cell-type annotation
* **Scrublet** — doublet detection
* **Pandas** — data manipulation
* **NumPy** — numerical computing
* **Matplotlib / Seaborn / Plotly** — visualization

## Design Philosophy

BioWorkbench is not intended to replace specialized analysis pipelines or computational workflows.

Instead, it provides an **interactive analytical workbench** on top of those workflows.

The intended architecture is:

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

This separation allows computationally intensive or reproducibility-critical analysis to remain in dedicated workflows while BioWorkbench provides an accessible interface for downstream exploration.

The same architecture can support both research use and the delivery of computational analyses to collaborators and clients.

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
* [x] Gene expression visualizer
* [x] Gene-set / module-score visualization
* [x] Categorical domain visualization
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
* [x] Categorical spatial visualization
* [x] Continuous spatial visualization
* [ ] Xenium workflows
* [ ] CosMX workflows
* [X] Spatial domains
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

1. Clone the repository and create a Python environment:

```bash
git clone https://github.com/poloarol/BioWorkbench.git
cd BioWorkbench

python -m venv env
```

Activate the environment and install the project dependencies.

Then launch the application with:

```bash
streamlit run app.py --server.maxUploadSize=500
```

2. Build docker container from scratch

```bash
docker build -t bioworkbench .
```

The launch the application with

```bash
docker run -p 8501:8501 bioworkbench --server.maxUploadSize=500
```

3. Download docker container form DockerHub

```bash
coming soon ...
```

## Project Status

BioWorkbench is an evolving project.

The current implementation provides the foundation for interactive single-cell analysis and is expanding toward a broader workbench for **single-cell, cell-resolved spatial, and high-resolution spatial genomics**.

Development is centered around:

* Modular computational components
* Platform-aware configuration
* Reusable visualization functions
* Standardized biological data structures
* Separation between analytical workflows and exploratory interfaces

The long-term goal is to provide an interactive workbench that can sit on top of validated computational analyses and make their results easier for researchers, collaborators, and clients to explore.

## License

BioWorkbench is free and open-source software released under the
GNU General Public License, version 3 or any later version
(GPL-3.0-or-later).

Copyright © 2026 Paul Wambo.

You are free to use, study, modify, and redistribute BioWorkbench
under the terms of the GNU GPL. If you distribute BioWorkbench or
a modified version of it, you must comply with the applicable
requirements of the GPL, including providing recipients with the
corresponding source code and preserving the applicable copyright
and license notices.

BioWorkbench is provided WITHOUT WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the `LICENSE` file for the complete terms and conditions.

For the full license text, see [LICENSE](LICENSE).

## Citation

If you use BioWorkbench in academic research, publications,
presentations, or other scholarly work, please cite the software
using the information provided in [`CITATION.cff`](CITATION.cff).

Citation information will be updated as BioWorkbench develops and
as associated publications and persistent software identifiers
become available.

## Author

**Paul A. Wambo**

Computational scientist working at the intersection of computational biology, machine learning, and scientific software.

GitHub: [@poloarol](https://github.com/poloarol)

---

*BioWorkbench is developed as an open computational tool for exploratory analysis and visualization of single-cell and spatial genomics data.*
