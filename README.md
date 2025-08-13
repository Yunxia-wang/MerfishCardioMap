# MerfishCardioMap

**3D spatial organization of heterogeneous progenitors in the zebrafish heart field pre-patterns cardiovascular development**

---

## Overview
This repository contains analysis pipelines, image-processing scripts, and data-handling utilities for reconstructing and analyzing 3D spatial transcriptomic maps from MERFISH experiments of the zebrafish heart field. The project aims to uncover how heterogeneous progenitor populations are spatially organized during early cardiovascular development.

---

## Features
- **Image Processing**
  - Mosaic splitting and replicate cropping for large `.tif` images
  - Alignment of DAPI and additional fluorescent channels
  - Support for Open-ST’s STIM and PASTE2-based spatial alignment
- **3D Reconstruction**
  - Reconstruction of embryo/tissue sections into volumetric spatial maps
  - Handling of multiple replicates per region
- **Cellular Analysis**
  - Extraction of cell coordinates and transcript counts from aligned images
  - Downstream clustering, annotation, and visualization
- **Visualization Tools**
  - Overlays for replicate boundary identification
  - 3D rendering of spatially aligned datasets

---

## Requirements
- Python ≥ 3.9
- Key Python packages:
  - `numpy`
  - `pandas`
  - `matplotlib`
  - `scikit-image`
  - `tifffile`
  - `napari`
  - [Open-ST](https://github.com/junjuew/Open-ST) (for STIM alignment)
  - [PASTE2](https://github.com/raphael-group/paste2) (for coordinate/image alignment)

Install dependencies:
```bash
pip install -r requirements.txt
