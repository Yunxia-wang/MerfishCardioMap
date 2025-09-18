# %%
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq

import os
import decoupler as dc
import monkeybread as mb

# ===================================================================
# Configuration with Relative Paths
# ===================================================================
# Define the base directory for all project-related files.
# This makes the script portable. All other paths are relative to this one.
BASE_DIR = Path('../data')

# Define the list of filenames (regions) to process.
filenames = [
    'slide1region1', 'slide1region5', 'slide1region7', 'slide1region9',
    'slide2region2', 'slide2region5', 'slide1region4', 'slide1region6',
    'slide1region8', 'slide2region0', 'slide2region4'
]

# Define the output directories relative to the base directory.
COMBINED_OUTPUT_DIR = BASE_DIR / '01_combined_cells_001'
CLUSTER_PLOTS_DIR = BASE_DIR / '02_combined_cells_cluster_plot_001'
SPATIAL_PLOTS_DIR = BASE_DIR / '03_plot_spatial_for_cluster'

# Create the output directories if they don't exist.
os.makedirs(COMBINED_OUTPUT_DIR, exist_ok=True)
os.makedirs(CLUSTER_PLOTS_DIR, exist_ok=True)
os.makedirs(SPATIAL_PLOTS_DIR, exist_ok=True)

# %%
print("Starting combined cell analysis...")

# Helper function to add a tab character for formatting
def deal_str(data):
    return str(data) + '\t'

# ===================================================================
# Part 1: Combine Data from Selected Regions
# ===================================================================
print("\nProcessing and combining data from selected regions...")

for index, filename in enumerate(filenames):
    print('===========================================')
    print(f"Processing region: {filename}")
    
    vizgen_dir = BASE_DIR / filename

    datapath_signal = vizgen_dir / 'sum_signals.csv'
    data_signal = pd.read_csv(datapath_signal)

    # Sort and filter cells based on 'Anti-Rabbit_high_pass' signal
    data_signal_sort = data_signal.sort_values(by='Anti-Rabbit_high_pass', ascending=False)
    print(f"Before cell selection: {data_signal.shape}")
    
    data_signal_sort_pos = data_signal_sort[data_signal_sort['Anti-Rabbit_high_pass'] != 0]
    data_signal_sel = data_signal_sort_pos.iloc[0:int(data_signal_sort_pos.shape[0]*0.9), :]
    print(f"After cell selection: {data_signal_sel.shape}")

    # Process and save cell metadata
    datapath_meta = vizgen_dir / 'cell_metadata.csv'
    datameta = pd.read_csv(datapath_meta)
    datameta_sel = datameta.loc[data_signal_sel.index, :]
    datameta_sel.to_csv(vizgen_dir / 'cell_metadata01.csv', index=False)
    print("Sample metadata after selection:")
    print(datameta_sel.head())

    # Process and save cell-by-gene data
    datapath = vizgen_dir / 'cell_by_gene.csv'
    datacell = pd.read_csv(datapath)
    datacell_sel = datacell.loc[data_signal_sel.index, :]
    datacell_sel.to_csv(vizgen_dir / 'cell_by_gene01.csv', index=False)
    print("Sample cell-by-gene data after selection:")
    print(datacell_sel.head())
    
    # Concatenate dataframes
    if index == 0:
        cell_combine = datacell_sel
        meta_combine = datameta_sel
    else:
        cell_combine = pd.concat([cell_combine, datacell_sel], axis=0)
        meta_combine = pd.concat([meta_combine, datameta_sel], axis=0)        

# Apply formatting to specific columns before saving combined files
meta_combine['EntityID'] = meta_combine['EntityID'].map(deal_str)
cell_combine['cell'] = cell_combine['cell'].map(deal_str)
        
cell_combine.to_csv(COMBINED_OUTPUT_DIR / 'cell_by_gene_combine.csv', index=False)
meta_combine.to_csv(COMBINED_OUTPUT_DIR / 'cell_metadata_combine.csv', index=False)
print(f"Combined cell-by-gene shape: {cell_combine.shape}")
print(f"Combined metadata shape: {meta_combine.shape}")

# Copy the required transformation file to the output directory
os.makedirs(COMBINED_OUTPUT_DIR / 'images', exist_ok=True)
source_path = BASE_DIR / 'slide1region1' / 'images' / 'micron_to_mosaic_pixel_transform.csv'
destination_path = COMBINED_OUTPUT_DIR / 'images' / 'micron_to_mosaic_pixel_transform.csv'
if source_path.exists():
    Path(source_path).rename(destination_path)
    print(f"Copied {source_path.name} to {destination_path.parent}")
else:
    print(f"Warning: Source file not found: {source_path}")

# ===================================================================
# Part 2: Read Combined Data and Perform Clustering
# ===================================================================
print("\nReading combined data into an AnnData object...")
adata = sq.read.vizgen(
    path=COMBINED_OUTPUT_DIR,
    counts_file="cell_by_gene_combine.csv",
    meta_file="cell_metadata_combine.csv",
    transformation_file="micron_to_mosaic_pixel_transform.csv",
)

sc.pp.calculate_qc_metrics(adata, percent_top=(10, 50, 100, 130), inplace=True)
print(f"Before cell filter: {adata.X.shape}")

# Plotting quality control metrics
print("Generating QC plots...")
fig, axs = plt.subplots(1, 3, figsize=(19, 5))
data_meta_df = pd.read_csv(COMBINED_OUTPUT_DIR / 'cell_metadata_combine.csv')

axs[0].set_title("Total transcripts per cell")
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])

axs[1].set_title("Unique transcripts per cell")
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, ax=axs[1])

axs[2].set_title("Volume of segmented cells")
sns.histplot(data_meta_df["volume"], kde=False, ax=axs[2])

plt.savefig(CLUSTER_PLOTS_DIR / 'the_distribution_of_total_transcripts.png', dpi=600)
print(f"Saved QC plot to: {CLUSTER_PLOTS_DIR / 'the_distribution_of_total_transcripts.png'}")

# Filter and preprocess the AnnData object
print("Filtering and preprocessing data...")
sc.pp.filter_cells(adata, min_counts=3, inplace=True)
sc.pp.filter_genes(adata, min_cells=1, inplace=True)
print(f"After cell filter: {adata.X.shape}")

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(
    adata,
    key_added="clusters",
    resolution=0.5,
    n_iterations=2,
    directed=False,
)

# Generate and visualize distinct colors for clusters
print("Generating cluster colors...")
len_cluster = len(list(set(adata.obs['clusters'].tolist())))

def generate_distinct_colors(num_colors):
    colors = plt.cm.get_cmap('tab20', num_colors)
    return [colors(i) for i in range(num_colors)]

def rgba_to_hex(rgba):
    return '#{:02x}{:02x}{:02x}'.format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))

num_colors = len_cluster
distinct_colors = generate_distinct_colors(num_colors)
hex_colors = [rgba_to_hex(color) for color in distinct_colors]
print(f"Generated {len_cluster} hex colors.")

plt.figure(figsize=(10, 2))
for i in range(num_colors):
    plt.plot(i, 0, marker='o', markersize=20, color=hex_colors[i])
plt.xlim(-1, num_colors)
plt.yticks([])
plt.title(f"{len_cluster} Distinct Colors from tab20 Colormap")
plt.show()

custom_colors = dict(zip(list(range(len_cluster)), hex_colors))

# Plot UMAP and spatial scatter plots for combined data
print("Generating UMAP and spatial plots...")
adata.uns['clusters_colors'] = list(custom_colors.values())

fig,ax=plt.subplots(figsize=(10,8))
sc.pl.umap(adata, color=['clusters'],
            legend_fontsize=10,
            legend_fontoutline=2,
            title="Clustering of pooled cells",
            show=False, ax=ax)
plt.savefig(CLUSTER_PLOTS_DIR / 'cluster_umap_pooled_cells.png', dpi=600)
print(f"Saved UMAP plot to: {CLUSTER_PLOTS_DIR / 'cluster_umap_pooled_cells.png'}")

sq.pl.spatial_scatter(
    adata,
    shape=None,
    color=["clusters"],
    wspace=0.4,
    size=8,
    return_ax=True,
)
plt.savefig(CLUSTER_PLOTS_DIR / 'Visualize_annotation_spatial_coordinates.png', dpi=600)
print(f"Saved spatial plot to: {CLUSTER_PLOTS_DIR / 'Visualize_annotation_spatial_coordinates.png'}")

# Save the final metadata with cluster information
print("Saving final metadata with cluster information...")
cell_meta_with_cluster = pd.DataFrame(adata.obs)
cell_meta_with_cluster['EntityID'] = cell_meta_with_cluster.index
cell_meta_with_cluster['EntityID'] = cell_meta_with_cluster['EntityID'].map(deal_str)
cell_meta_with_cluster.to_csv(COMBINED_OUTPUT_DIR / 'cell_metadata_combine_with_cluster.csv', index=False)
print(f"Saved final metadata to: {COMBINED_OUTPUT_DIR / 'cell_metadata_combine_with_cluster.csv'}")


# ===================================================================
# Part 3: Visualize Spatial Clusters for Individual Regions
# ===================================================================
print("\nPlotting spatial clusters for individual regions...")
annotated_cell = pd.read_csv(COMBINED_OUTPUT_DIR / 'cell_metadata_combine_with_cluster.csv')
dict_cell_clu = annotated_cell.set_index('EntityID')['clusters'].to_dict()
cell_list = [int(i.strip()) for i in annotated_cell['EntityID'] ]
custom_colors[len_cluster] = '#f1f2f1' # Add color for unselected cells

for index, filename in enumerate(filenames):
    print('===========================================')
    print(f"Plotting spatial clusters for region: {filename}")

    vizgen_dir = BASE_DIR / filename
    resultpath = SPATIAL_PLOTS_DIR / filename
    os.makedirs(resultpath, exist_ok=True)

    # Generate new cell_metadata file with cluster information for the current region
    datapath_meta = vizgen_dir / 'cell_metadata.csv'
    datameta = pd.read_csv(datapath_meta)
    
    cell_clusters = []
    count_cell = 0
    for i in datameta['EntityID'].tolist():
        if i in cell_list:
            count_cell += 1
            cell_clusters.append(dict_cell_clu[str(i) + '\t'])
        else:
            cell_clusters.append(str(len_cluster))
    
    datameta.loc[:, ('cluster')] = cell_clusters
    datameta.to_csv(vizgen_dir / 'cell_metadata_pooled.csv', index=False)
    print(f"Number of cells in combined data found in this region: {count_cell}")
    
    # Read data for the specific region
    adata = sq.read.vizgen(
        path=vizgen_dir,
        counts_file="cell_by_gene.csv",
        meta_file="cell_metadata_pooled.csv",
        transformation_file="micron_to_mosaic_pixel_transform.csv",
    )
    
    # Preprocessing
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pp.filter_cells(adata, min_counts=3, inplace=True)
    sc.pp.filter_genes(adata, min_cells=1, inplace=True)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    
    cell_types = [str(i) for i in adata.obs['cluster'].tolist()]
    adata.obs['cell_type'] = pd.Categorical(cell_types)
    
    custom_colors_region = {}
    for cell_type in set(cell_types):
        custom_colors_region[cell_type] = custom_colors[int(cell_type)]

    adata.uns['cell_type_colors'] = list(custom_colors_region.values())

    # Plot spatial figure for all clusters
    fig, axs = plt.subplots(1, 1, figsize=(16, 20))
    sc.pl.embedding(
        adata,
        "spatial",
        color='cell_type',
        groups=cell_types,
        size=60,
        ax=axs,
        palette=custom_colors_region,
        show=False
    )
    axs.yaxis.set_inverted(True)
    plt.savefig(resultpath / f'{filename}_all_cluster_spatial.png', dpi=600)
    
    # Plot spatial figure for each cluster
    plot_cluster = set(cell_types)
    if str(len_cluster) in plot_cluster:
        plot_cluster.remove(str(len_cluster))
    
    for i in list(plot_cluster):
        fig, axs = plt.subplots(1, 1, figsize=(16, 20))
        sc.pl.embedding(
            adata,
            "spatial",
            color='cell_type',
            groups=str(i),
            size=60,
            ax=axs,
            palette=custom_colors_region,
            show=False,
            na_color='#f1f2f1'
        )
        axs.yaxis.set_inverted(True)
        plt.savefig(resultpath / f'{filename}_cluster{i}_spatial.png', dpi=600)

print("Script execution complete.")
