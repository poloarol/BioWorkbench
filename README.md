# BioWorkbench

**BioWorkbench** is a lightweight interactive workbench for single-cell and spatial transcriptomics analysis.

It provides a Streamlit-based interface for common preprocessing, clustering, visualization, and cell-type annotation workflows while keeping the underlying analysis logic modular and reusable.

The project currently supports conventional single-cell transcriptomics workflows and an emerging imaging-based spatial transcriptomics workflow, with **MERFISH** currently supported.

> **🚧 Active development**
>
> BioWorkbench is under active development. The current spatial workflow is focused on MERFISH, with support for additional spatial technologies planned.

---

## Features

### Single-cell transcriptomics

* AnnData (`.h5ad`) input
* Cell and gene filtering
* Quality-control metrics
* Mitochondrial gene filtering
* Doublet detection using Scrublet
* Highly variable gene selection
* Normalization and log transformation
* PCA
* Neighborhood graph construction
* UMAP
* Leiden clustering
* Marker-gene identification
* CellTypist-based cell-type annotation
* Custom annotation support
* Interactive visualization of annotations and QC metrics

### Spatial transcriptomics

BioWorkbench is being extended to support cell-resolved and high-resolution spatial transcriptomics workflows.

**Currently supported:**

* MERFISH

The current MERFISH workflow includes:

* Cell-level spatial coordinates
* Spatial quality-control metrics
* Blank-transcript filtering
* Cell morphology and segmentation-derived metadata
* Spatial visualization
* Spatial domain identification
* Visualization of annotations in both UMAP and tissue coordinates

**Planned spatial technologies:**

* Xenium
* CosMX
* Visium HD

Additional spatial functionality, including image overlays, segmentation visualization, transcript-level visualization, and spatial statistics, is planned.

---

## Analysis workflow

BioWorkbench separates the user interface from the underlying analysis functions.

```text
                    BioWorkbench
                         │
             ┌───────────┴───────────┐
             │                       │
        Single-cell               Spatial
        transcriptomics          transcriptomics
             │                       │
             │                    MERFISH
             │                       │
             └───────────┬───────────┘
                         │
              Shared analysis layer
                         │
          ┌──────────────┼──────────────┐
          │              │              │
         QC        Visualization   Annotation
          │              │              │
          └──────────────┴──────────────┘
                         │
                       AnnData
```

The application is organized around three primary analysis stages:

```text
Filtering → Clustering → Annotation
```

---

## Input data

BioWorkbench uses **AnnData** objects as its primary data structure.

Single-cell datasets should contain a gene-expression matrix in `adata.X` together with the associated cell and gene metadata.

Spatial datasets additionally require cell-level spatial information appropriate to the selected technology.

For the current MERFISH workflow, spatial coordinates are represented using:

```text
adata.obs["center_x"]
adata.obs["center_y"]
```

and are standardized internally through:

```text
adata.obsm["spatial"]
```

MERFISH datasets may also contain platform-specific metadata such as blank-transcript counts, segmentation measurements, cell area, and cell volume.

---

## Quality control

BioWorkbench does not apply identical QC procedures to every data modality.

### Single-cell datasets

The current single-cell workflow supports parameters including:

* minimum genes per cell
* minimum cells per gene
* maximum mitochondrial gene percentage
* expected doublet rate

Doublet detection is performed using Scrublet.

### Imaging-based spatial datasets

Imaging-based spatial technologies have different QC characteristics from conventional scRNA-seq.

The current MERFISH workflow therefore emphasizes platform-specific measurements such as:

* genes detected per cell
* transcript counts
* blank transcripts
* percentage of blank transcripts
* cell area
* cell volume
* segmentation-derived measurements

Conventional mitochondrial and ribosomal QC metrics are not assumed to be available or appropriate for every spatial technology.

---

## Configuration

Default filtering parameters are stored in:

```text
config/params.yaml
```

Parameters are organized by dataset modality and, where appropriate, spatial technology.

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

Only technologies currently exposed by the application are considered supported.

Configuration values provide default parameters for the interactive workflow and can be adjusted through the application where supported.

---

## Clustering and spatial domains

For single-cell datasets, BioWorkbench provides a conventional analysis workflow based on:

```text
Normalization
      ↓
Highly variable genes
      ↓
PCA
      ↓
Nearest-neighbor graph
      ↓
UMAP
      ↓
Leiden clustering
      ↓
Marker genes
```

For supported spatial datasets, BioWorkbench can additionally incorporate spatial proximity when identifying spatial domains.

The current spatial-domain workflow combines expression-based and spatial neighborhood information into a joint graph before community detection.

This allows spatial organization to be considered alongside transcriptional similarity.

---

## Visualization

BioWorkbench provides visualization of both continuous and categorical measurements.

### UMAP

UMAP visualizations can be colored by cell-level metadata, including:

* clustering assignments
* cell-type annotations
* QC metrics
* other `adata.obs` variables

Categorical variables are displayed using discrete colors and legends, while continuous variables use continuous color scales.

### Spatial visualization

For supported spatial datasets, cell-level measurements can also be visualized using tissue coordinates.

Spatial visualizations support both:

* categorical annotations
* continuous measurements

For spatial datasets, annotation results can therefore be examined in both:

```text
UMAP space
    +
tissue space
```

This makes it possible to compare transcriptional organization with spatial organization.

---

## Cell-type annotation

BioWorkbench currently supports cell-type annotation using **CellTypist**.

Available workflows include:

* CellTypist model selection
* automated prediction
* majority voting
* prediction confidence
* custom annotations

Annotation results are stored in the AnnData object and can be visualized on UMAP embeddings and, for supported spatial datasets, in tissue coordinates.

---

## Marker genes

Marker genes can be identified using Scanpy's `rank_genes_groups` workflow.

The current implementation uses a t-test-based method and provides statistics including:

* marker-gene rankings
* scores
* log fold changes
* adjusted p-values

Marker analysis is currently performed with respect to an existing categorical observation, such as a clustering assignment.

---

## Project structure

```
BioWorkbench/
├── app.py                         # Main Streamlit application
│
├── pages/                         # Streamlit application pages
│   ├── 01_filtering.py           # QC and filtering
│   ├── 02_clustering.py          # Clustering and spatial domains
│   └── 03_annotation.py          # Cell-type annotation
|   |__ 04_visualization.py
|   |__ 05_download.py
│
│
├── src/                           # Reusable analysis and visualization code
│   ├── filtering.py               # Filtering and QC functions
│   ├── clustering.py              # Clustering and spatial-domain analysis
│   ├── plotting.py                # UMAP, spatial, and QC plots
│   └── utils.py                   # Data loading and shared utilities
│
├── config/                        # Application configuration
│   └── params.yaml                # Modality/technology-specific parameters
│
├── notebooks/                     # Exploratory analysis and development
│
├── requirements.txt               # Python dependencies
├── LICENSE                        # GNU GPLv3-or-later
├── CITATION.cff                   # Software citation metadata
├── README.md                      # Project documentation
└── .gitignore
```

The `pages/` directory contains the Streamlit application pages, while reusable analysis and visualization functions are maintained in `src/`.

Configuration defaults are maintained separately in `config/`.

---

## Roadmap

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
* [x] MERFISH spatial input
* [x] MERFISH-specific QC
* [x] Spatial coordinate visualization
* [x] Spatial annotation visualization
* [x] Spatial-domain identification
* [x] Spatially variable gene analysis

### In development

* [ ] Improved spatial-domain analysis
* [ ] Improved spatial QC visualization
* [ ] More robust annotation workflows
* [ ] Expanded spatial data validation
* [ ] Additional spatial visualization capabilities

### Planned

* [ ] Xenium support
* [ ] CosMX support
* [ ] Visium HD support
* [ ] Tissue image overlays
* [ ] Segmentation overlays
* [ ] Transcript-level visualization
* [ ] Additional export functionality

---

## License

BioWorkbench is free software released under the **GNU General Public License, version 3 or any later version (GPL-3.0-or-later)**.

Copyright © 2026 Paul Wambo.

You are free to use, study, modify, and redistribute BioWorkbench under the terms of the GNU GPL. If you distribute BioWorkbench or a modified version of it, you must comply with the applicable requirements of the GPL, including the requirements concerning source code and preservation of applicable copyright and license notices.

BioWorkbench is provided **WITHOUT WARRANTY**, without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, to the extent permitted by applicable law.

See the [`LICENSE`](LICENSE) file for the complete license terms.

---

## Citation

If you use BioWorkbench in academic research, publications, presentations, or other scholarly work, please cite the software using the information provided in [`CITATION.cff`](CITATION.cff).

Citation information will be updated as BioWorkbench develops and associated publications and persistent software identifiers become available.

---

## Disclaimer

BioWorkbench is research software intended to support exploratory and computational analysis of biological datasets.

It is not intended to replace validated clinical, diagnostic, regulatory, or other professionally validated analytical workflows.

Users are responsible for evaluating the suitability, accuracy, and reproducibility of results generated using the software for their intended application.
