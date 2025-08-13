# MerfishCardioMap — Multipolygon Visualization Script

**3D spatial organization of heterogeneous progenitors in the zebrafish heart field pre-patterns cardiovascular development**

This repository contains Python code for visualizing **2D spatial polygon geometries** (e.g., segmented cells, tissue regions) using Plotly, and exporting the results to PDF.  
The script reads polygon geometries from a dataset (e.g., `.parquet`), assigns colors based on associated metadata, and plots them with equal aspect ratio for accurate spatial representation.

---

## Features
- **Reads MultiPolygon and Polygon geometries** from geospatial datasets.
- Handles both **single polygons** and **multiple polygons per entity**.
- **Closes polygon loops** automatically to ensure correct rendering.
- **Custom colors** for each polygon (e.g., based on cluster or annotation).
- **High-quality PDF export** using Plotly + Kaleido.
- Maintains **equal scaling** of X and Y axes for accurate spatial visualization.

---

## Dependencies
Make sure the following Python packages are installed:

```bash
pip install pandas geopandas shapely plotly kaleido

