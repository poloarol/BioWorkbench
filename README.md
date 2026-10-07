![BioWorkbench logo](docs/images/bioworkbench-logo.png)

**An interactive computational workbench for single-cell and spatial transcriptomics.**

BioWorkbench is a lightweight, open-source [Streamlit](https://streamlit.io/) application for exploring, analyzing, and visualizing single-cell and spatial omics datasets through an interactive interface.

It is built around **scverse** ecosystem, while keeping the underlying computational functions separate from the user interface. The goal is to make common exploratory analyses easier to perform, inspect, and reproduce without requiring users to repeatedly modify analysis scripts.

> 🚧 **BioWorkbench is currently in beta.**
>
> The current beta focuses on conventional single-cell RNA-seq and MERFISH spatial transcriptomics. The project is actively evolving, and feedback from researchers using real datasets is welcome.

---

## Why BioWorkbench?

Single-cell and spatial transcriptomics analyses often involve a combination of notebooks, scripts, visualization tools, and modality-specific workflows.

BioWorkbench is intended to sit between a full computational pipeline and a purely visual data explorer:

```text
                 AnnData (.h5ad)
                       │
                       ▼
              ┌─────────────────┐
              │ Quality Control  │
              │   & Filtering    │
              └────────┬────────┘
                       │
                       ▼
              ┌─────────────────┐
              │ Normalization & │
              │  HVG Selection  │
              └────────┬────────┘
                       │
                       ▼
              ┌─────────────────┐
              │       PCA       │
              └────────┬────────┘
                       │
                       ▼
              ┌─────────────────┐
              │ UMAP + Leiden   │
              │    Clustering   │
              └────────┬────────┘
                       │
                       ▼
              ┌─────────────────┐
              │ Marker Genes &  │
              │ Cluster Analysis│
              └────────┬────────┘
                       │
                       ▼
              ┌─────────────────┐
              │   CellTypist    │
              │   Annotation    │
              └────────┬────────┘
                       │
                       ▼
                Annotated AnnData
```

For spatial datasets, the workflow extends this representation with spatial quality control, coordinates, spatial domains, and tissue-level visualization.

The long-term goal is not to replace established analysis frameworks such as Scanpy, Squidpy, or Seurat. Instead, BioWorkbench aims to provide a **lightweight interactive layer for exploration, visualization, annotation, and downstream analysis**.

---

# ✨ Current capabilities

## 🔬 Single-cell RNA-seq

BioWorkbench currently provides a conventional Scanpy-based single-cell workflow including:

* AnnData (`.h5ad`) input
* Cell and gene filtering
* Quality-control metric calculation
* Mitochondrial gene filtering
* Scrublet-based doublet detection
* Highly variable gene selection
* Total-count normalization
* Log transformation
* PCA
* Nearest-neighbor graph construction
* UMAP
* Leiden clustering
* Marker-gene analysis
* Cluster-level summaries
* CellTypist-based cell-type annotation
* Custom annotation workflows
* Interactive visualization
* Processed AnnData export

### Configurable analysis parameters

Users can interactively adjust parameters including:

* Minimum genes per cell
* Minimum cells per gene
* Maximum mitochondrial percentage
* Expected doublet rate
* Number of highly variable genes
* Number of principal components
* Number of neighbors
* UMAP `min_dist`
* Leiden resolution

---

## 🧬 MERFISH spatial transcriptomics

BioWorkbench currently includes an imaging-based spatial workflow for **MERFISH** datasets.

The MERFISH workflow incorporates modality-specific information that is not captured by conventional scRNA-seq QC.

Current functionality includes:

* Cell-level spatial coordinates
* Transcript counts per cell
* Genes detected per cell
* Blank-transcript measurements
* Percentage of blank transcripts
* Cell morphology and segmentation-derived metadata
* Cell area and volume
* Spatial quality-control visualization
* Spatial filtering
* Spatial-domain identification
* Spatial visualization
* Cell-type annotation in tissue coordinates
* Gene-expression visualization in tissue space
* Module-score visualization in tissue space
* Comparison of transcriptional and spatial organization

The spatial workflow uses the AnnData representation:

```text
adata.obsm["spatial"]
```

for standardized cell coordinates.

This allows the same analyzed cells to be explored from two complementary perspectives:

```text
          UMAP                         Tissue
           │                             │
           ▼                             ▼
   transcriptional                 physical/spatial
    organization                    organization
           │                             │
           └──────────────┬──────────────┘
                          ▼
                 biological context
```

This is particularly useful for asking whether transcriptional clusters, annotations, or gene-expression patterns correspond to spatial organization within the tissue.

---

# 🔍 Analysis workflow

A typical single-cell analysis can follow:

```text
1. Load AnnData
       ↓
2. Inspect QC metrics
       ↓
3. Filter cells and genes
       ↓
4. Detect potential doublets
       ↓
5. Retain singlets
       ↓
6. Normalize expression
       ↓
7. Select highly variable genes
       ↓
8. Run PCA
       ↓
9. Construct neighborhood graph
       ↓
10. Run UMAP
       ↓
11. Run Leiden clustering
       ↓
12. Identify marker genes
       ↓
13. Annotate cell types
       ↓
14. Explore annotations and expression
       ↓
15. Export processed AnnData
```

For MERFISH, spatial QC and spatial analysis can be incorporated into the workflow.

---

# 🏷️ Cell-type annotation

BioWorkbench integrates [CellTypist](https://www.celltypist.org/) for model-based cell-type annotation.

Users can:

1. Browse available CellTypist models
2. Select a reference model
3. Run automated annotation
4. Apply majority voting
5. Inspect prediction confidence
6. Visualize predictions on UMAP
7. Explore annotations in spatial coordinates for supported datasets

Annotation results are stored in `adata.obs`.

The current workflow records CellTypist outputs such as:

```text
celltypist_predicted_labels
celltypist_majority_voting
celltypist_conf_score
```

BioWorkbench also supports custom annotation workflows using JSON input.

---

# 📊 Visualization

BioWorkbench supports interactive visualization of categorical and continuous cell-level measurements.

### UMAP

UMAPs can be colored by variables stored in `adata.obs`, including:

* Cluster assignments
* Cell-type annotations
* QC metrics
* Gene expression
* Other cell-level metadata

### Spatial visualization

For supported spatial datasets, measurements can be displayed using tissue coordinates.

Both categorical and continuous variables are supported, allowing users to visualize:

* Cell types
* Clusters
* QC metrics
* Gene expression
* Module scores
* Other cell-level metadata

The same biological variable can therefore be examined in both transcriptional and spatial contexts.

---

# 💾 Export

BioWorkbench can package selected analysis parameters and AnnData objects from the current Streamlit session.

### AnnData files

Selected AnnData objects are included as standard:

```text
.h5ad
```

files. Analysis results stored in each selected AnnData object are retained in that file.

### BioWorkbench export bundle

The Download page creates a ZIP archive named `bioworkbench.wkb` containing a manifest, selected parameters in `params.json`, and selected AnnData files under `adatas/`.

This bundle is an export format, not a restorable application session. The application currently does not import `.wkb` bundles. To continue an analysis later, retain the exported `.h5ad` files and upload the desired dataset when you return.

The export bundle can be used to:

* Save an intermediate analysis state
* Transfer an analysis between environments
* Preserve multiple selected AnnData objects and associated parameters

The ZIP-based format is intentionally used rather than Python object serialization so that the session data remains inspectable and portable.

> **Note:** For interoperability with other tools, `.h5ad` remains the preferred data-exchange format.

---

# 📥 Input data

BioWorkbench uses **AnnData** as its primary data structure.

A typical dataset is represented as:

```text
AnnData
├── X        expression matrix
├── obs      cell metadata
├── var      gene metadata
├── obsm     embeddings / spatial coordinates
├── uns      unstructured analysis metadata
└── layers   alternative expression matrices
```

At minimum, a conventional single-cell dataset should contain a gene-expression matrix in:

```python
adata.X
```

along with appropriate cell and gene metadata.

### MERFISH input

The MERFISH workflow expects cell-level spatial information and may use platform-specific measurements such as:

* Cell coordinates
* Segmentation measurements
* Cell morphology
* Blank-transcript counts
* Cell area
* Cell volume

Spatial coordinates are standardized internally through:

```python
adata.obsm["spatial"]
```

The current spatial workflow is configured for MERFISH. The input must include `center_x` and `center_y` in `adata.obs`; blank-probe genes are identified by names beginning with `Blank-`. Spatial QC plots also use the cell `volume` metadata column. Input requirements may evolve as additional spatial technologies are introduced.

---

# ⚙️ Configuration

Default analysis parameters are maintained in:

```text
config/params.yaml
```

Configuration is organized by modality and, where appropriate, spatial technology.

Conceptually:

```yaml
single_cell:
  ...

spatial:
  MERFISH:
    ...
```

Configuration values provide defaults for the interactive workflows and can be adjusted through the application where supported.

---

# 🧪 Beta testing

BioWorkbench is currently in **beta** and is being tested with real single-cell and spatial transcriptomics datasets.

The current beta focuses on:

* Single-cell RNA-seq
* MERFISH spatial transcriptomics
* Quality control and filtering
* Dimensionality reduction and clustering
* Marker-gene analysis
* Cell-type annotation
* Gene and module-score visualization
* Spatial visualization
* Spatial-domain exploration
* Exporting selected parameters and AnnData objects

## Who is this for?

BioWorkbench is primarily intended for:

* Computational biologists
* Bioinformaticians
* Single-cell researchers
* Spatial transcriptomics researchers
* Researchers working with AnnData-based datasets

You do not need to be a BioWorkbench developer to participate.

If you regularly work with `.h5ad` datasets and single-cell or spatial omics data, your feedback is particularly valuable.

## What are we looking for?

Beta testers are encouraged to use BioWorkbench with real datasets and report:

* Bugs or unexpected behavior
* Installation problems
* Environment/dependency issues
* Confusing workflows
* UI/UX problems
* Visualization issues
* Unexpected analysis results
* Performance problems
* Missing functionality
* Features that would make the workflow more useful

You do **not** need to determine whether something is technically a "bug" before reporting it.

If something is confusing, surprising, or difficult to use, that is useful feedback.

## Reporting issues

Please open a GitHub Issue and include, where possible:

1. What you were trying to do
2. What you expected to happen
3. What actually happened
4. Steps to reproduce the problem
5. Operating system
6. Installation method
7. Dataset characteristics
8. Relevant error messages
9. Screenshots, where useful

---

# 🚀 Installation

Clone the repository:

```bash
git clone https://github.com/poloarol/BioWorkbench.git
cd BioWorkbench
```

Create a virtual environment:

```bash
python -m venv .venv
```

### macOS / Linux

```bash
source .venv/bin/activate
```

### Windows

```powershell
.venv\Scripts\Activate.ps1
```

Install the dependencies:

```bash
pip install -r requirements.txt
```

---

# ▶️ Running BioWorkbench

Start the Streamlit application:

```bash
streamlit run app.py --server.maxUploadSize=500
```

The application will be available at `http://localhost:8501`. The upload widget currently limits individual `.h5ad` uploads to 500 MB. When changing the Streamlit server upload limit, note that the widget also has its own limit.

Upload an `.h5ad` dataset through the application to begin the analysis workflow.

---

# 🐳 Running with Docker

BioWorkbench includes a Dockerfile for running the application in a containerized environment.

Build the image:

```bash
docker build -t bioworkbench .
```

Run the application (the container listens on port 8501):

```bash
docker run --rm -p 8501:8501 bioworkbench
```

The image runs as an unprivileged user and includes a health check. Then open:

```text
http://localhost:8501
```

Docker provides an alternative to installing the Python dependencies directly on the host system and is particularly useful for reproducible beta testing.

---

# 🧪 Testing

BioWorkbench uses `pytest` for automated testing of its computational and utility functions, plotting behavior, and Streamlit application startup.

Run the complete test suite:

```bash
python -m pytest
```

Run a specific test file:

```bash
python -m pytest tests/<test_file>.py
```

Run a specific test:

```bash
python -m pytest tests/<test_file>.py -k "<test_name>"
```

The test suite validates reusable analysis workflows, dataset preparation and export, plotting behavior, and Streamlit application startup, including checks that workflow pages handle missing prerequisites. Tests use small synthetic datasets rather than large biological samples.

When modifying computational functionality, corresponding tests should be added or updated where appropriate.

---

# 🏗️ Project structure

```text
BioWorkbench/
│
├── app.py                         # Main Streamlit application
│
├── pages/
│   ├── 01_filtering.py            # QC and filtering
│   ├── 02_clustering.py           # PCA, UMAP, clustering, spatial domains
│   ├── 03_annotation.py            # Cell-type annotation
│   ├── 04_visualization.py        # Gene/module visualization
│   └── 05_download.py             # Parameter and AnnData export
│
├── src/
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
├── docs/
│   └── images/                    # Documentation assets
│
├── beta                          # Beta-testing guidance
├── requirements.txt               # Python dependencies
├── Dockerfile                     # Container configuration
├── CITATION.cff                   # Software citation metadata
├── LICENSE                        # GNU GPL-3.0-or-later
└── README.md
```

The application interface is separated from the reusable computational layer.

This allows the analysis functions to evolve independently of Streamlit and makes them potentially reusable in:

* Scripts
* Notebooks
* Pipelines
* Automated workflows
* Future interfaces

---

# 🧰 Technology stack

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
| Containerization     | Docker                      |

---

# 🗺️ Roadmap

BioWorkbench is being developed incrementally, with emphasis on making the existing workflows robust before expanding to additional technologies.

## Current

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
* [x] Custom annotation
* [x] Processed AnnData export
* [x] MERFISH input
* [x] MERFISH-specific QC
* [x] Spatial coordinate visualization
* [x] Spatial annotation visualization
* [x] Spatial-domain identification
* [x] Spatially variable gene analysis
* [x] Transcript-level visualization
* [x] ZIP-based export of selected parameters and AnnData objects
* [x] Automated tests
* [x] Containerized execution

## In development

* [ ] Improved spatial-domain analysis
* [ ] Improved spatial QC visualization
* [ ] Expanded spatial data validation
* [ ] More robust annotation workflows
* [ ] Additional spatial visualization capabilities
* [ ] Reload BioWorkbench `.wkb` export bundles
* [ ] Improved reproducibility and configuration management
* [ ] Expanded documentation
* [ ] Continuous integration

## Planned

### Additional spatial technologies

* [ ] Xenium
* [ ] CosMX
* [ ] Visium HD

### Spatial analysis

* [ ] Spatial marker analysis
* [ ] Cell-neighborhood analysis
* [ ] Cell-cell interaction analysis

### Annotation

* [ ] Spatially aware cell-type annotation
* [ ] Multiple annotation models
* [ ] Marker-based annotation
* [ ] Interactive annotation editing
* [ ] Annotation confidence and uncertainty visualization

### Intelligent analysis

* [ ] Context-aware gene and gene-set information retrieval
* [ ] AI-assisted biological interpretation

---

# 🧠 Design philosophy

## Interactive analysis

BioWorkbench allows researchers to explore datasets and adjust common analysis parameters without repeatedly modifying analysis scripts.

## Modularity

Analysis functions are separated from the Streamlit interface so that computational components can be reused independently of the user interface.

## AnnData-first workflows

AnnData serves as the central representation for expression data, metadata, embeddings, spatial coordinates, and analysis results.

## Modality-aware analysis

Different biological technologies have different data structures and quality-control considerations.

BioWorkbench therefore aims to share common computational components where appropriate while exposing modality-specific workflows where necessary.

## Reproducibility

Configuration, reusable analysis functions, automated tests, containerization, and portable session bundles are intended to make computational behavior easier to inspect, validate, and reproduce.

---

# 🔭 Future direction

The longer-term goal is to extend BioWorkbench from a single-cell analysis application into a lightweight computational workbench for **single-cell and spatial biology**.

```text
                           BioWorkbench
                                 │
                 ┌───────────────┴───────────────┐
                 │                               │
           Single-cell                     Spatial Biology
                 │                               │
          ┌──────┴──────┐              ┌─────────┼─────────┐
          │             │              │         │         │
       scRNA-seq     Other          MERFISH    Xenium   Other
                     assays                       spatial
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

The objective is not to replace established single-cell and spatial analysis frameworks.

Instead, BioWorkbench aims to provide a practical interactive layer that helps researchers **inspect data, explore patterns, test hypotheses, and prepare analyses for more specialized computational workflows**.

---

# 📖 Citation

If you use BioWorkbench in research, please cite the repository using the information provided in [`CITATION.cff`](CITATION.cff).

---

# 📄 License

BioWorkbench is free software released under the **GNU General Public License v3.0 or later (GPL-3.0-or-later)**.

You are free to use, study, modify, and redistribute BioWorkbench under the terms of the GNU GPL.

See [`LICENSE`](LICENSE) for the complete license terms.

---

# ⚠️ Disclaimer

BioWorkbench is research software under active development.

It is intended to support exploratory and computational analysis and should not be considered a validated clinical diagnostic tool.

Analysis results should be independently reviewed and interpreted in the context of the underlying dataset, experimental design, and appropriate biological knowledge.
