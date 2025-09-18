# ---------------------------------------------------------------------
# Import necessary libraries
# ---------------------------------------------------------------------
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq

import os

# ---------------------------------------------------------------------
# Define base paths and filenames using Path objects
# ---------------------------------------------------------------------

# Define the base directory for the entire project. This assumes the script
# is run from the directory containing the 'Hakan_SP' folder.
# Adjust this path if your folder structure is different.
base_path = Path('./Hakan_SP')

# The primary path for input and output data
pri_path = Path('../data')

# The path for pre-clustered cell data
cluster_datapath = pri_path / '04_combined_cells_iter1'

# The output path for plots
plot_result = '08_combined_cells_cluster_plot_iter3_umap2/'
output_path_plot_base = pri_path / plot_result

# Filenames for data files
cell_meta_with_cluster_file = cluster_datapath / 'cell_metadata_combine_with_cluster_iter3.csv'
cell_by_gene_file = 'cell_by_gene.csv'
cell_metadata_pooled_file = 'cell_metadata_pooled_iter3.csv'
transform_file = 'micron_to_mosaic_pixel_transform.csv'

# ---------------------------------------------------------------------
# Load and prepare data and colors
# ---------------------------------------------------------------------

# Get the cell metadata with cluster information
annotated_cell = pd.read_csv(cell_meta_with_cluster_file)
dict_cell_clu = annotated_cell.set_index('EntityID')['clusters'].to_dict()

num_clustered = len(annotated_cell['clusters'].unique())
cell_list = [int(i) for i in annotated_cell['EntityID']]

def generate_distinct_colors(num_colors):
    """Generates a list of distinct colors from the 'Set1' colormap."""
    colors = plt.cm.get_cmap('Set1', num_colors)
    return [colors(i) for i in range(num_colors)]

def rgba_to_hex(rgba):
    """Converts an RGBA color tuple to a hex string."""
    return f'#{int(rgba[0]*255):02x}{int(rgba[1]*255):02x}{int(rgba[2]*255):02x}'

# Generate colors for the clustered cells
distinct_colors = generate_distinct_colors(num_clustered)
hex_colors = [rgba_to_hex(color) for color in distinct_colors]

custom_colors = {0: '#53b64b', 1: '#75157c', 2: '#f2a93b', 3: '#6e9af8', 4: '#e87e70'}
# Add a gray color for unclustered cells
custom_colors[num_clustered] = '#cccccc'

# Plot the colors for visualization
plt.figure(figsize=(12, 2))
for i in range(len(custom_colors)):
    plt.plot(i, 0, marker='o', markersize=20, color=custom_colors[i])
plt.xlim(-1, len(custom_colors))
plt.yticks([])
plt.title(f'{num_clustered} Distinct Colors from tab20 Colormap')
plt.show()

# ---------------------------------------------------------------------
# Main analysis loop over each region
# ---------------------------------------------------------------------
filenames = ['slide1region5','slide1region7','slide1region9','slide2region2','slide2region5','slide1region4','slide1region6','slide1region8','slide2region0','slide2region4']

for filename in filenames:
    print('=' * 40)
    print(f'Processing {filename}')
    
    # Define paths for the current region
    vizgen_dir = pri_path / filename
    resultpath = output_path_plot_base / filename
    
    # Create the result directory
    os.makedirs(resultpath, exist_ok=True)

    # Generate a new cell metadata file with cluster info
    datapath_meta = vizgen_dir / 'cell_metadata.csv'
    datameta = pd.read_csv(datapath_meta)
    
    cell_clusters = []
    for entity_id in datameta['EntityID'].tolist():
        # Map existing clusters or assign to 'unclustered' group
        cluster_id = dict_cell_clu.get(entity_id, num_clustered)
        cell_clusters.append(cluster_id)
    
    datameta['cluster'] = cell_clusters
    datameta.to_csv(vizgen_dir / cell_metadata_pooled_file, index=False)
    
    # Load data into an AnnData object using Squidpy
    adata = sq.read.vizgen(
        path=vizgen_dir,
        counts_file=cell_by_gene_file,
        meta_file=cell_metadata_pooled_file,
        transformation_file=transform_file,
    )

    # QC and Filtering
    sc.pp.calculate_qc_metrics(adata, percent_top=(10, 50, 70, 90), inplace=True)
    sc.pp.filter_cells(adata, min_counts=3, inplace=True)
    sc.pp.filter_genes(adata, min_cells=1, inplace=True)
    
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    
    # Prepare data for plotting
    adata.obs['cell_type'] = pd.Categorical([str(i) for i in adata.obs['cluster']])
    
    custom_colors_region = {
        str(c): custom_colors[c] for c in sorted(adata.obs['cluster'].unique())
    }
    
    adata.uns['cell_type_colors'] = list(custom_colors_region.values())

    # -------------------------------
    # Plot all clusters in one figure
    # -------------------------------
    fig, axs = plt.subplots(1, 1, figsize=(16, 20))
    sc.pl.embedding(
        adata,
        "spatial",
        color='cell_type',
        groups=list(adata.obs['cell_type'].unique()),
        wspace=0.4,
        size=7,
        ax=axs,
        palette=custom_colors_region,
        show=False
    )
    axs.yaxis.set_inverted(True)
    plt.savefig(resultpath / f'{filename}_all_cluster_spatial.png', dpi=600)
    plt.close(fig)

    # -------------------------------    
    # Plot each cluster individually
    # -------------------------------
    plot_cluster_ids = [str(i) for i in custom_colors if i != num_clustered]
    for cluster_id in plot_cluster_ids:
        fig, axs = plt.subplots(1, 1, figsize=(16, 20))
        sc.pl.embedding(
            adata,
            "spatial",
            color='cell_type',
            groups=cluster_id,
            wspace=0.4,
            ax=axs,
            size=7,
            palette=custom_colors_region,
            show=False,
            na_color='#cccccc'
        )
        axs.yaxis.set_inverted(True)
        plt.savefig(resultpath / f'{filename}_cluster{cluster_id}_spatial.png', dpi=600)
        plt.close(fig)
