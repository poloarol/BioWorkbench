# 🧬 BioWorkbench

BioWorkbench is a modular, interactive **single-cell and spatial transcriptomics analysis workbench** built with Python and Streamlit.

It provides an interactive interface for common computational biology workflows while keeping the underlying analysis logic in reusable Python modules. The application is built around **AnnData** and the Scanpy ecosystem and is designed to support both exploratory analysis and reproducible computational workflows.

> 🚧 **Active development**
>
> BioWorkbench currently supports conventional single-cell transcriptomics and an emerging imaging-based spatial workflow, with **MERFISH currently supported**. Additional spatial technologies are under development.

---

## Overview

BioWorkbench is designed to sit between a full computational pipeline and a purely visual data explorer.

The goal is to provide a lightweight interface for:

* quality control and filtering
* dimensionality reduction
* clustering
* marker-gene analysis
* cell-type annotation
* interactive visualization
* spatial analysis and visualization

while keeping the underlying computational functions independent of the user interface.

The core workflow is organized around an AnnData object:

```text
                         AnnData (.h5ad)
                                │
                                ▼
                     ┌─────────────────────┐
                     │  Quality Control &  │
                     │      Filtering      │
                     └──────────┬──────────┘
                                │
                                ▼
                     ┌─────────────────────┐
                     │ Normalization & HVG │
                     │     Selection       │
                     └──────────┬──────────┘
                                │
                                ▼
                     ┌─────────────────────┐
                     │         PCA         │
                     └──────────┬──────────┘
                                │
                                ▼
                     ┌─────────────────────┐
                     │  UMAP + Leiden      │
                     │     Clustering      │
                     └──────────┬──────────┘
                                │
                                ▼
                     ┌─────────────────────┐
                     │ Marker Genes &      │
                     │ Cluster Analysis    │
                     └──────────┬──────────┘
                                │
                                ▼
                     ┌─────────────────────┐
                     │    CellTypist       │
                     │    Annotation       │
                     └──────────┬──────────┘
                                │
                                ▼
                         Annotated AnnData
```

For supported spatial datasets, spatial information can be incorporated into quality control, visualization, and spatial-domain analysis.

---

## Features

### 🔬 Single-cell transcriptomics

BioWorkbench currently supports a conventional Scanpy-based single-cell workflow including:

* AnnData (`.h5ad`) input
* cell and gene filtering
* quality-control metrics
* mitochondrial gene filtering
* Scrublet-based doublet detection
* highly variable gene selection
* normalization and log transformation
* PCA
* nearest-neighbor graph construction
* UMAP
* Leiden clustering
* marker-gene analysis
* CellTypist-based cell-type annotation
* custom annotations
* interactive QC and analysis visualization
* processed AnnData export

### 🧬 Spatial transcriptomics

BioWorkbench is being extended to support cell-resolved and high-resolution spatial transcriptomics workflows.

**Currently supported:**

* MERFISH

The current MERFISH workflow includes:

* cell-level spatial coordinates
* spatial quality-control metrics
* blank-transcript filtering
* cell morphology and segmentation-derived metadata
* spatial visualization
* spatial-domain identification
* visualization of annotations in tissue coordinates
* comparison of transcriptional and spatial organization

**Planned spatial technologies:**

* Xenium
* CosMX
* Visium HD

Additional spatial functionality is planned, including tissue-image overlays, segmentation visualization, transcript-level visualization, and spatial statistics.

---

## Analysis Workflow

### 1. Quality Control & Filtering

For conventional single-cell datasets, BioWorkbench provides configurable filtering parameters including:

* minimum genes detected per cell
* minimum cells expressing a gene
* maximum mitochondrial gene percentage
* expected doublet rate

Doublet detection is performed using **Scrublet**.

The application maintains intermediate AnnData objects so that users can inspect different stages of filtering, including predicted doublets and retained singlets.

For imaging-based spatial datasets, QC is handled differently because conventional scRNA-seq metrics are not necessarily applicable.

The MERFISH workflow can incorporate measurements such as:

* genes detected per cell
* transcript counts
* blank transcripts
* percentage of blank transcripts
* cell area
* cell volume
* segmentation-derived measurements

---

### 2. Dimensionality Reduction & Clustering

BioWorkbench provides a standard Scanpy-based workflow:

```text
Normalization
      ↓
Highly Variable Gene Selection
      ↓
PCA
      ↓
Nearest-Neighbor Graph
      ↓
UMAP
      ↓
Leiden Clustering
      ↓
Marker-Gene Analysis
```

Relevant parameters can be adjusted interactively, including:

* number of highly variable genes
* number of principal components
* number of neighbors
* UMAP `min_dist`
* Leiden resolution

The resulting embeddings and cluster assignments are stored in the AnnData object.

---

### 3. Marker-Gene Analysis

Marker genes can be identified for existing categorical groups, such as Leiden clusters.

The current implementation uses Scanpy's `rank_genes_groups` workflow with a **t-test-based method**.

Results include:

* marker-gene rankings
* scores
* log fold changes
* adjusted p-values

The application also provides cluster-level summaries such as cell counts and proportions.

---

### 4. Cell-Type Annotation

BioWorkbench integrates **CellTypist** for model-based cell-type annotation.

Users can:

1. Browse available CellTypist models
2. Select a reference model
3. Run automated annotation
4. Apply majority voting
5. Inspect prediction confidence
6. Visualize annotations on UMAP

Annotation results are stored in `adata.obs`.

The application also supports custom annotation workflows using JSON input.

---

### 5. Visualization

BioWorkbench supports visualization of both categorical and continuous cell-level measurements.

#### UMAP

UMAPs can be colored by variables stored in `adata.obs`, including:

* cluster assignments
* cell-type annotations
* QC metrics
* other metadata

#### Spatial Visualization

For supported spatial datasets, cell-level measurements can be visualized using tissue coordinates.

Both categorical and continuous variables are supported.

For example:

```text
                 UMAP
                  │
                  │
            transcriptional
             organization
                  │
                  ▼
             tissue space
                  │
                  │
              spatial
             organization
```

This allows users to examine whether transcriptional clusters or cell-type annotations correspond to spatial organization.

---

### 6. Export

Processed AnnData objects can be exported from the application as `.h5ad` files.

Marker-gene results can also be exported from the clustering workflow.

---

## Input Data

BioWorkbench uses **AnnData** as its primary data structure.

At minimum, a single-cell dataset should contain a gene-expression matrix in:

```python
adata.X
```

along with appropriate cell and gene metadata.

The general AnnData structure is:

```text
AnnData
├── X        expression matrix
├── obs      cell metadata
├── var      gene metadata
├── obsm     embeddings / coordinates
├── uns      unstructured metadata
└── layers   alternative expression matrices
```

### MERFISH

The current MERFISH workflow expects cell-level spatial coordinates.

Coordinates are standardized internally through:

```python
adata.obsm["spatial"]
```

MERFISH datasets may additionally contain platform-specific information such as:

```text
cell coordinates
cell morphology
segmentation measurements
blank-transcript counts
cell area
cell volume
```

Input requirements may evolve as support for additional spatial technologies is added.

---

## Configuration

Default analysis parameters are stored in:

```text
config/params.yaml
```

Configuration is organized by dataset modality and, where appropriate, spatial technology.

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

Only technologies currently exposed by the application should be considered supported.

Configuration values provide defaults for the interactive workflow and can be adjusted through the application where supported.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/poloarol/BioWorkbench.git
cd BioWorkbench
```

Create a virtual environment:

```bash
python -m venv .venv
```

Activate the environment.

### macOS / Linux

```bash
source .venv/bin/activate
```

### Windows

```powershell
.venv\Scripts\Activate.ps1
```

Install the project dependencies:

```bash
pip install -r requirements.txt
```

> **Note:** BioWorkbench is under active development and its dependency specification may change as the project matures.

---

## Running BioWorkbench

Start the Streamlit application with:

```bash
streamlit run app.py
```

The application will open in your browser.

Upload an `.h5ad` dataset through the application to begin the analysis workflow.

---

## Testing

BioWorkbench uses **pytest** for automated testing of its computational functionality.

Run the complete test suite from the repository root:

```bash
python -m pytest
```

To run a specific test file:

```bash
python -m pytest tests/<test_file>.py
```

To run a specific test or group of tests:

```bash
python -m pytest tests/<test_file>.py -k "<test_name>"
```

The test suite is intended to provide deterministic checks of the reusable analysis and utility functions while keeping tests independent of large biological datasets where possible.

When adding or modifying computational functionality, corresponding tests should be added or updated where appropriate.

---

## Project Structure

```text
BioWorkbench/
│
├── app.py                         # Main Streamlit application
│
├── pages/                         # Streamlit application pages
│   ├── 01_filtering.py            # QC and filtering
│   ├── 02_clustering.py           # Clustering and spatial domains
│   ├── 03_annotation.py           # Cell-type annotation
│   ├── 04_visualization.py        # Gene/module visualization
│   └── 05_download.py             # Data and result export
│
├── src/                           # Reusable computational code
│   ├── filtering.py               # Filtering and doublet detection
│   ├── clustering.py              # Clustering and spatial analysis
│   ├── plotting.py                # Visualization functions
│   └── utils.py                   # Data loading and shared utilities
│
├── config/
│   └── params.yaml                # Analysis configuration
│
├── tests/                         # Automated tests
│
├── notebooks/                     # Exploratory analysis and development
│
├── requirements.txt               # Python dependencies
├── Dockerfile                     # Container configuration
├── CITATION.cff                   # Software citation metadata
├── LICENSE                        # GNU GPL-3.0-or-later
└── README.md                      # Project documentation
```

The application interface is separated from the reusable computational layer.

This allows analysis functions to be reused independently of Streamlit in scripts, notebooks, pipelines, or future interfaces.

---

## Technology Stack

| Component            | Technology                  |
| -------------------- | --------------------------- |
| Language             | Python                      |
| Application          | Streamlit                   |
| Data structure       | AnnData                     |
| Single-cell analysis | Scanpy                      |
| Cell-type annotation | CellTypist                  |
| Doublet detection    | Scrublet                    |
| Data manipulation    | Pandas                      |
| Visualization        | Matplotlib, Plotly, Seaborn |
| Testing              | pytest                      |

---

## Roadmap

BioWorkbench is under active development. The roadmap is divided into implemented functionality, work in development, and longer-term planned functionality.

### Implemented

* [x] AnnData input
* [x] Single-cell QC
* [x] Cell and gene filtering
* [x] Mitochondrial filtering
* [x] Scrublet doublet detection
* [x] Highly variable gene selection
* [x] PCA
* [x] UMAP
* [x] Leiden clustering
* [x] Marker-gene analysis
* [x] CellTypist annotation
* [x] Custom annotation workflow
* [x] Processed AnnData export
* [x] MERFISH spatial input
* [x] MERFISH-specific QC
* [x] Spatial coordinate visualization
* [x] Spatial annotation visualization
* [x] Spatial-domain identification
* [x] Automated tests

### In Development

* [ ] Improved spatial-domain analysis
* [ ] Improved spatial QC visualization
* [ ] More robust annotation workflows
* [ ] Expanded spatial data validation
* [ ] Additional spatial visualization capabilities
* [ ] Improved reproducibility and configuration management

### Planned

#### Additional Spatial Technologies

* [ ] Xenium support
* [ ] CosMX support
* [ ] Visium HD support

#### Spatial Analysis

* [ ] Spatially variable gene analysis
* [ ] Spatial marker analysis
* [ ] Tissue-image overlays
* [ ] Segmentation overlays
* [ ] Transcript-level spatial visualization
* [ ] Cell neighborhood analysis
* [ ] Cell-cell interaction analysis

#### Annotation

* [ ] Spatially aware cell-type annotation
* [ ] Multiple annotation models
* [ ] Marker-based annotation
* [ ] Interactive annotation editing
* [ ] Annotation confidence and uncertainty visualization

#### Intelligent Analysis

* [ ] Context-aware gene and gene-set information retrieval
* [ ] Agentic AI-assisted biological interpretation

#### Software Engineering

* [ ] Containerized deployment
* [ ] Improved session-state management
* [ ] Modular dataset abstractions
* [ ] Reproducible analysis configurations
* [ ] Continuous integration
* [ ] Expanded documentation

---

## Design Philosophy

### Interactive Analysis

BioWorkbench allows users to explore datasets and adjust common analysis parameters without repeatedly modifying analysis scripts.

### Modularity

Analysis functions are separated from the Streamlit interface so that computational components can be reused in scripts, notebooks, pipelines, or future applications.

### AnnData-First Workflows

AnnData serves as the central representation for expression data, metadata, embeddings, and analysis results.

### Modality-Aware Analysis

Different biological technologies have different data structures and quality-control considerations.

BioWorkbench therefore aims to share common computational components where appropriate while exposing modality-specific functionality where necessary.

### Reproducibility

Configuration, reusable analysis functions, and automated tests are intended to make the computational behaviour of BioWorkbench easier to inspect, validate, and reproduce.

---

## Future Direction

The longer-term goal is to extend BioWorkbench from a single-cell analysis interface into a lightweight computational workbench for **single-cell and spatial biology**.

The architecture is intended to allow different modalities to share common analytical components while maintaining modality-specific workflows.

```text
                           BioWorkbench
                                │
                ┌───────────────┴───────────────┐
                │                               │
          Single-cell                     Spatial Biology
                │                               │
          ┌─────┴─────┐                ┌────────┼────────┐
          │           │                │        │        │
       scRNA-seq   Other            MERFISH   Xenium   Other
                   assays                              spatial
                │                               │
                └───────────────┬───────────────┘
                                │
                                ▼
                     Shared computational
                         components
                                │
                  ┌─────────────┼─────────────┐
                  │             │             │
                 QC       Visualization   Annotation
                  │             │             │
                  └─────────────┼─────────────┘
                                │
                                ▼
                       Biological insight
```

The objective is not to replace established single-cell or spatial analysis frameworks. Instead, BioWorkbench aims to provide a lightweight interactive layer for exploration, visualization, annotation, and downstream computational analysis.

---

## Project Status

> 🚧 **Active development**

BioWorkbench is currently functional for conventional single-cell transcriptomics workflows and has an emerging MERFISH spatial workflow.

The codebase, interfaces, APIs, and supported technologies may change as development continues.

---

## License

BioWorkbench is free software released under the **GNU General Public License, version 3 or any later version (GPL-3.0-or-later)**.

Copyright © 2026 Paul Wambo.

You are free to use, study, modify, and redistribute BioWorkbench under the terms of the GNU GPL. If you distribute BioWorkbench or a modified version of it, you must comply with the applicable requirements of the GPL.

See the [`LICENSE`](LICENSE) file for the complete license terms.

---

## Citation

If you use BioWorkbench in academic research, publications, presentations, or other scholarly work, please cite the software using the metadata provided in [`CITATION.cff`](CITATION.cff).

Citation information will be updated as the project develops and persistent software identifiers or associated publications become available.

---

## Disclaimer

BioWorkbench is research software intended to support exploratory and computational analysis of biological datasets.

It is **not intended to replace validated clinical, diagnostic, regulatory, or other professionally validated analytical workflows**.

Users are responsible for evaluating the suitability, accuracy, and reproducibility of results generated using the software for their intended application.

---

## Author

**Paul Wambo**

Computational scientist working across bioinformatics, machine learning, and computational biology.

GitHub: [@poloarol](https://github.com/poloarol)
