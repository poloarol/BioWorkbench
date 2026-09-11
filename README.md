# 🧬 BioWorkbench

**BioWorkbench** is a lightweight, interactive Streamlit application for exploratory analysis of single-cell transcriptomics data.

The application takes an **AnnData (`.h5ad`) object** as input and provides an interactive workflow for quality control, filtering, dimensionality reduction, clustering, visualization, marker-gene exploration, and model-based cell-type annotation using [CellTypist](https://www.celltypist.org/).

The project is being developed as a modular computational biology workbench, with planned support for **spatial transcriptomics and spatial imaging datasets**, including **10x Visium, Visium HD, MERFISH, and related spatial modalities**.

## Project History

BioWorkbench is a Python/Streamlit rewrite and extension of an earlier
R Shiny application for interactive single-cell transcriptomics analysis.

**Previous implementation:**  
[Single-Cell Transcriptomics App](https://github.com/poloarol/single-cell-transcriptomics-app)

The original application provided interactive quality control, dimensionality
reduction, clustering, visualization, and cell-type annotation. BioWorkbench
reimplements these capabilities using Python and the AnnData/Scanpy ecosystem
while introducing a modular architecture designed for future single-cell and
spatial transcriptomics workflows.


---

## Overview

BioWorkbench is designed to make common single-cell analysis tasks accessible through a lightweight interactive interface while keeping the underlying analysis logic in reusable Python modules.

The current workflow is:

```text
                    AnnData (.h5ad)
                           │
                           ▼
                  ┌─────────────────┐
                  │  Quality Control │
                  │   & Filtering    │
                  └────────┬────────┘
                           │
                           ▼
                  ┌─────────────────┐
                  │  Normalization   │
                  │  & HVG selection │
                  └────────┬────────┘
                           │
                           ▼
                  ┌─────────────────┐
                  │       PCA        │
                  └────────┬────────┘
                           │
                           ▼
                  ┌─────────────────┐
                  │ UMAP + Leiden   │
                  │   Clustering    │
                  └────────┬────────┘
                           │
                           ▼
                  ┌─────────────────┐
                  │ Marker Genes &  │
                  │ Cluster Summary │
                  └────────┬────────┘
                           │
                           ▼
                  ┌─────────────────┐
                  │    CellTypist   │
                  │   Annotation    │
                  └────────┬────────┘
                           │
                           ▼
                 Annotated AnnData
```

---

## Current Features

### 🔬 Quality Control & Filtering

Upload an `.h5ad` file and interactively configure basic single-cell QC and filtering parameters.

Current filtering options include:

* Minimum genes detected per cell
* Minimum cells expressing a gene
* Maximum mitochondrial gene percentage
* Expected doublet rate
* Scrublet-based doublet detection

The application retains intermediate AnnData objects so that users can compare:

* Input/subset data
* Filtered data including predicted doublets
* Filtered singlet data

QC summaries include:

* Number of cells
* Genes detected per cell
* Counts per cell
* Highly expressed genes
* Doublet/singlet statistics
* QC visualizations

---

### 📊 Dimensionality Reduction & Clustering

BioWorkbench provides an interactive implementation of a standard Scanpy-based analysis workflow.

#### PCA

The PCA workflow performs:

1. Total-count normalization
2. Log transformation
3. Highly variable gene selection
4. PCA

Parameters can be adjusted interactively, including:

* Number of highly variable genes
* Number of principal components

#### UMAP & Leiden clustering

Following PCA, users can run:

* Nearest-neighbor graph construction
* UMAP
* Leiden clustering

Configurable parameters include:

* Number of neighbors
* UMAP `min_dist`
* Leiden resolution

The resulting embeddings and cluster labels are stored directly in the AnnData object.

---

### 🧬 Marker Gene Analysis

BioWorkbench automatically identifies marker genes for detected clusters.

The current implementation uses:

* Cluster-level differential expression
* t-test-based ranking
* Adjusted p-values
* Log fold changes
* Marker scores

The interface displays the top marker genes for each cluster and provides a cluster-level summary including cell counts and proportions.

---

### 🏷️ Cell-Type Annotation

BioWorkbench integrates **CellTypist** for automated cell-type annotation.

Users can:

1. Browse available CellTypist reference models
2. Select a reference model
3. Run model-based annotation
4. Apply CellTypist majority voting
5. Inspect annotation confidence scores
6. Visualize predictions on UMAP embeddings

The resulting annotations are stored in `adata.obs`, including:

```text
celltypist_predicted_labels
celltypist_majority_voting
celltypist_conf_score
```

The application also supports uploading custom JSON-based annotations for additional annotation workflows.

---

### 💾 Export

Processed AnnData objects can be exported directly from the application as `.h5ad` files.

Marker-gene results can also be exported from the clustering workflow.

---

## Technology Stack

BioWorkbench is built primarily with the Python scientific-computing ecosystem.

| Component            | Technology                  |
| -------------------- | --------------------------- |
| Application          | Streamlit                   |
| Data structure       | AnnData                     |
| Single-cell analysis | Scanpy                      |
| Cell-type annotation | CellTypist                  |
| Doublet detection    | Scrublet                    |
| Data manipulation    | Pandas                      |
| Visualization        | Matplotlib, Plotly, Seaborn |
| Language             | Python                      |

---

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
└── notebooks/
```

### `app.py`

Defines the Streamlit application and navigation between analysis modules.

### `pages/`

Contains the user-facing Streamlit workflows:

* `01_filtering.py` — QC, filtering, and doublet detection
* `02_clustering.py` — PCA, UMAP, Leiden clustering, and marker analysis
* `03_annotation.py` — CellTypist and custom annotation

### `src/`

Contains reusable computational functions separated from the Streamlit interface:

* `filtering.py` — cell/gene filtering and doublet detection
* `clustering.py` — normalization, HVG selection, PCA, UMAP, clustering, and marker identification
* `plotting.py` — visualization functions
* `utils.py` — data loading, QC metric calculation, and AnnData utilities

This separation is intended to allow the computational components to evolve independently of the user interface.

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

Activate it on macOS/Linux:

```bash
source .venv/bin/activate
```

On Windows:

```powershell
.venv\Scripts\activate
```

Install the required Python packages:

```bash
pip install streamlit scanpy celltypist scrublet pandas matplotlib seaborn plotly
```

> A dedicated dependency file will be added as the project matures.

---

## Running BioWorkbench

Start the Streamlit application with:

```bash
streamlit run app.py
```

The application will open in your browser.

Upload an `.h5ad` file from the **Cell and Gene Filtering** page to begin the analysis.

---

## Input Data

The current version expects an **AnnData `.h5ad` file** containing a gene-expression matrix.

The basic data model is:

```text
AnnData
├── X       expression matrix
├── obs     cell metadata
├── var     gene metadata
├── obsm    embeddings
├── uns     unstructured analysis metadata
└── layers  alternative expression matrices
```

BioWorkbench calculates common QC metrics and maintains the AnnData object throughout the analysis workflow.

---

## Example Workflow

A typical analysis might look like:

```text
1. Upload .h5ad
       ↓
2. Inspect QC metrics
       ↓
3. Filter low-quality cells and genes
       ↓
4. Identify potential doublets
       ↓
5. Remove predicted doublets
       ↓
6. Normalize expression data
       ↓
7. Identify highly variable genes
       ↓
8. Run PCA
       ↓
9. Construct neighborhood graph
       ↓
10. Run UMAP
       ↓
11. Run Leiden clustering
       ↓
12. Identify cluster marker genes
       ↓
13. Annotate cell types with CellTypist
       ↓
14. Inspect annotations on UMAP
       ↓
15. Export processed AnnData
```

---

## Roadmap

BioWorkbench is currently focused on **single-cell transcriptomics**, but the project is intended to expand into a broader computational workbench for single-cell and spatial biology.

### Current

* [x] AnnData input
* [x] QC metric calculation
* [x] Cell filtering
* [x] Gene filtering
* [x] Mitochondrial filtering
* [x] Doublet detection
* [x] PCA
* [x] Highly variable gene selection
* [x] UMAP
* [x] Leiden clustering
* [x] Marker-gene identification
* [x] CellTypist annotation
* [x] Custom annotations
* [x] Processed AnnData export

### Planned

#### Spatial transcriptomics

* [ ] Spatial AnnData loading
* [ ] Visium support
* [ ] Visium HD support
* [ ] Spatial quality control
* [ ] Tissue-coordinate visualization
* [ ] Spatial clustering/domain identification
* [ ] Spatially variable gene analysis
* [ ] Spatial marker analysis

#### Imaging-based spatial transcriptomics

* [ ] MERFISH support
* [ ] Xenium support
* [ ] Cell-coordinate visualization
* [ ] Segmentation-aware visualization
* [ ] Transcript-level spatial visualization
* [ ] Cell neighborhood analysis
* [ ] Cell-cell interaction analysis

#### Expanded annotation

* [ ] Spatially aware cell-type annotation
* [ ] Multiple annotation models
* [ ] Marker-based annotation
* [ ] Interactive annotation editing
* [ ] Annotation confidence and uncertainty visualization

#### Software engineering

* [ ] Formal dependency management
* [ ] Automated testing
* [ ] Containerized deployment
* [ ] Improved session-state management
* [ ] Modular dataset abstraction
* [ ] Reproducible analysis configurations
* [ ] Expanded documentation

---

## Design Philosophy

BioWorkbench is intended to sit between a full computational pipeline and a purely visual data explorer.

The goal is to provide:

**Interactive analysis**

Allow users to explore datasets and adjust common analysis parameters without repeatedly modifying scripts.

**Modularity**

Keep analysis functions separate from the user interface so that individual components can be reused in scripts, notebooks, pipelines, or future interfaces.

**AnnData-first workflows**

Use AnnData as the central representation for expression data and analysis results.

**Extensibility**

Design the application so that the same interface can eventually support multiple biological modalities rather than building independent applications for every data type.

---

## Future Direction

The longer-term goal is to extend BioWorkbench from a single-cell analysis interface into a lightweight computational workbench for **single-cell and spatial biology**.

The planned architecture will allow different modalities to share common analytical components while exposing modality-specific functionality where appropriate.

```text
                         BioWorkbench
                              │
              ┌───────────────┴───────────────┐
              │                               │
        Single-cell                      Spatial Biology
              │                               │
        ┌─────┴─────┐              ┌──────────┼──────────┐
        │           │              │          │          │
      scRNA-seq   Other          Visium    MERFISH    Xenium
                   assays         /HD
        │           │              │          │          │
        └───────────┴──────────────┴──────────┴──────────┘
                              │
                              ▼
                    Shared computational
                         components
                              │
                  ┌───────────┼───────────┐
                  │           │           │
                 QC      Visualization  Annotation
                  │           │           │
                  └───────────┼───────────┘
                              │
                              ▼
                       Biological insight
```

The objective is not to replace established single-cell or spatial analysis frameworks, but to provide a **lightweight interactive layer** for exploration, visualization, annotation, and downstream analysis.

---

## Status

🚧 **Active development**

The current release is primarily a **single-cell transcriptomics workbench**. Spatial functionality is planned and the codebase is being developed with future spatial data support in mind.

---

## Author

**Paul Wambo**

Computational scientist working across bioinformatics, machine learning, and computational biology.

GitHub: [@poloarol](https://github.com/poloarol)

---

## License

License information will be added as the project matures.
