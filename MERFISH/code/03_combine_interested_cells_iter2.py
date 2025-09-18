# %%
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq
import os
import shutil
import decoupler as dc
import monkeybread as mb

# ===================================================================
# Configuration with Relative Paths
# ===================================================================
# Define the base directory for all project-related files.
# Change this path to the root of your project.
BASE_DIR = Path('../data')


# Define the list of regions to process.
filenames = [
    'slide1region5', 'slide1region7', 'slide1region9', 'slide2region2',
    'slide2region5', 'slide1region4', 'slide1region6', 'slide1region8',
    'slide2region0', 'slide2region4'
]

# Define the output directories relative to the base directory.
OUTPUT_DIR = BASE_DIR / '04_combined_cells_iter1'
PLOT_OUTPUT_DIR = BASE_DIR / '07_combined_cells_cluster_plot_iter2_umap2'
COMBINED_CELLS_DIR = BASE_DIR / '01_combined_cells'

# Create the output directories if they don't exist.
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(PLOT_OUTPUT_DIR, exist_ok=True)

# Helper function to add a tab character for formatting
def deal_str(data):
    return str(data) + '\t'

# %%
# ===================================================================
# Part 1: Filter and Select Cells for Iteration 2
# ===================================================================
print("Starting iterative cell filtering...")

# Load initial combined metadata from a previous step
cellmeta_iter1 = pd.read_csv(OUTPUT_DIR / 'cell_metadata_combine_with_cluster_iter1.csv')
print(f"Loaded metadata shape (iteration 1): {cellmeta_iter1.shape}")

# Filter out clusters 5 and 6 as per the original script's logic
cellmeta_iter2 = cellmeta_iter1[~cellmeta_iter1['clusters'].isin([5, 6])]
print(f"Filtered metadata shape (iteration 2): {cellmeta_iter2.shape}")

# Load the original combined metadata and filter it to match the selected cells
cellmeta_combine_iter1 = pd.read_csv(OUTPUT_DIR / 'cell_metadata_combine.csv')
cellmeta_combine_iter2 = cellmeta_combine_iter1[cellmeta_combine_iter1['EntityID'].isin(cellmeta_iter2['EntityID'].tolist())]
cellmeta_combine_iter2.to_csv(OUTPUT_DIR / 'cell_metadata_combine_iter2.csv', index=False)
print(f"Saved filtered metadata for iteration 2 to: {OUTPUT_DIR / 'cell_metadata_combine_iter2.csv'}")

# %%
# Load the original cell-by-gene data and filter for the selected cells
cell_by_gene_combine_path = COMBINED_CELLS_DIR / 'cell_by_gene_combine.csv'
cell_by_gene_combine_data = pd.read_csv(cell_by_gene_combine_path)

selected_cell_lists = cellmeta_combine_iter2['EntityID'].tolist()
print(f"Filtering cell-by-gene data for {len(selected_cell_lists)} selected cells...")

cell_by_gene_com_new = cell_by_gene_combine_data[cell_by_gene_combine_data['cell'].isin(selected_cell_lists)]
cell_by_gene_com_new.to_csv(OUTPUT_DIR / 'cell_by_gene_combine_iter2.csv', index=False)
print(f"Saved filtered cell-by-gene data to: {OUTPUT_DIR / 'cell_by_gene_combine_iter2.csv'}")

# %%
# ===================================================================
# Part 2: Prepare Gene Panel and Copy Transformation File
# ===================================================================
print("\nPreparing gene panel and file system...")

# Load the gene panel for cell type annotation
gene_panel = BASE_DIR / 'marker_data' / 'makers_annotated_merfish.csv'
df_ref_panel_ini = pd.read_csv(gene_panel)
df_ref_panel_all = df_ref_panel_ini.iloc[1:, :2]
df_ref_panel_all.columns = ['Markers', 'cell_type']
print("Loaded reference gene panel:")
print(df_ref_panel_all.head())

marker_cell = dict(zip(df_ref_panel_all['Markers'], df_ref_panel_all['cell_type']))

# Filter the cell-by-gene data to include only the genes from the marker panel
list_127 = ['cell'] + df_ref_panel_all['Markers'].tolist()
cell_by_gene_com_new_127 = cell_by_gene_com_new[list_127]
cell_by_gene_com_new_127.to_csv(OUTPUT_DIR / 'cell_by_gene_combine_iter2.csv', index=False)
print(f"Saved cell-by-gene data with 127 marker genes to: {OUTPUT_DIR / 'cell_by_gene_combine_iter2.csv'}")

# Copy the required transformation file
source_path = BASE_DIR / 'cell_seg_all' / 'images' / 'micron_to_mosaic_pixel_transform.csv'
destination_path = OUTPUT_DIR / 'images' / 'micron_to_mosaic_pixel_transform.csv'
os.makedirs(destination_path.parent, exist_ok=True)
shutil.copyfile(source_path, destination_path)
print(f"Copied transformation file to: {destination_path}")

# %%
# ===================================================================
# Part 3: Read, Preprocess, and Cluster the Filtered Data
# ===================================================================
print("\nReading filtered data into AnnData object...")
adata = sq.read.vizgen(
    path=OUTPUT_DIR,
    counts_file="cell_by_gene_combine_iter2.csv",
    meta_file="cell_metadata_combine_iter2.csv",
    transformation_file="micron_to_mosaic_pixel_transform.csv",
)

sc.pp.calculate_qc_metrics(adata, percent_top=(10, 50, 100, 120), inplace=True)
print(f"Before cell filter: {adata.X.shape}")

sc.pp.filter_cells(adata, min_counts=3, inplace=True)
sc.pp.filter_genes(adata, min_cells=1, inplace=True)
print(f"After cell filter: {adata.X.shape}")

# Preprocess the adata
adata.layers["counts"] = adata.X.copy()
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=120)
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

# %%
# Generate distinct colors for clusters
print("Generating and plotting cluster colors...")
len_cluster = len(list(set(adata.obs['clusters'].tolist())))

def generate_distinct_colors(num_colors):
    colors = plt.cm.get_cmap('tab20', num_colors)
    return [colors(i) for i in range(num_colors)]

def rgba_to_hex(rgba):
    return '#{:02x}{:02x}{:02x}'.format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))

num_colors = len_cluster
distinct_colors = generate_distinct_colors(num_colors)
hex_colors = [rgba_to_hex(color) for color in distinct_colors]
custom_colors = dict(zip(list(range(len_cluster)), hex_colors))

plt.figure(figsize=(10, 2))
for i in range(num_colors):
    plt.plot(i, 0, marker='o', markersize=20, color=hex_colors[i])
plt.xlim(-1, num_colors)
plt.yticks([])
plt.title(f"{len_cluster} Distinct Colors from tab20 Colormap")
plt.show()

# %%
# Plot UMAP and spatial scatter plots for the new clustering
print("Plotting UMAP and spatial coordinates...")
adata.uns['clusters_colors'] = list(custom_colors.values())

fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(adata, color=['clusters'],
            legend_fontsize=10,
            legend_fontoutline=2,
            title="Clustering of pooled cells",
            show=False,
            ax=ax)
plt.savefig(PLOT_OUTPUT_DIR / 'cluster_umap_pooled_cells.png', dpi=600)
plt.close(fig)

fig, ax = plt.subplots(figsize=(10, 8))
sq.pl.spatial_scatter(
    adata,
    shape=None,
    color=["clusters"],
    wspace=0.4,
    size=8,
    return_ax=True,
    ax=ax
)
plt.savefig(PLOT_OUTPUT_DIR / 'Visualize_annotation_spatial_coordinates.png', dpi=600)
plt.close(fig)

# %%
# ===================================================================
# Part 4: Cell Type Annotation
# ===================================================================
print("\nPerforming cell type annotation...")

# Save the updated metadata with the new clustering results
cell_meta_with_cluster = pd.DataFrame(adata.obs)
cell_meta_with_cluster['EntityID'] = cell_meta_with_cluster.index.map(deal_str)
cell_meta_with_cluster.to_csv(OUTPUT_DIR / 'cell_metadata_combine_with_cluster_iter2.csv')
print(f"Saved updated metadata to: {OUTPUT_DIR / 'cell_metadata_combine_with_cluster_iter2.csv'}")

# Create marker gene dictionary for plotting
marker_genes_dict = {}
for cell_type in set(df_ref_panel_all['cell_type'].tolist()):
    sele_markers = [marker for marker in df_ref_panel_all['Markers'].tolist() if marker_cell.get(marker) == cell_type]
    marker_genes_dict[cell_type] = sele_markers

# Plot dotplot of marker gene expression
fig, ax = plt.subplots(figsize=(50, 8))
sc.pl.dotplot(adata, marker_genes_dict, "clusters", dendrogram=True, show=False, ax=ax)
plt.savefig(PLOT_OUTPUT_DIR / 'marker_expression_top_gene.png', dpi=300)
plt.close(fig)

# Perform Over-Representation Analysis (ORA) with decoupler
dc.run_ora(
    mat=adata,
    net=df_ref_panel_all,
    source='cell_type',
    target='Markers',
    min_n=3,
    verbose=True,
    use_raw=False
)
acts = dc.get_acts(adata, obsm_key='ora_estimate')

# Handle infinite values from ORA
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e

# Rank sources and map clusters to cell types
df = dc.rank_sources_groups(acts, groupby='clusters', reference='rest', method='t-test_overestim_var')
ctypes_dict = df.groupby('group').head(3).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
ctypes_dict_new = {key: value[0] for key, value in ctypes_dict.items()}

print("Mapped clusters to cell types:")
print(ctypes_dict_new)

# Map colors to cell types based on the original cluster colors
dict_cell_color = {}
for cell_ty in sorted(list(ctypes_dict_new.values())):
    for key, value in ctypes_dict_new.items():
        if cell_ty == value:
            dict_cell_color[value] = custom_colors[int(key)]

# Add the new cell type annotation to the AnnData object
acts.obs["cell type"] = acts.obs["clusters"].map(ctypes_dict_new).astype("category")
acts.uns['cell type_colors'] = list(dict_cell_color.values())

# Plot the UMAP with cell type annotation
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(
    acts,
    color="cell type",
    legend_fontsize=10,
    legend_fontoutline=2,
    title="Annotation of cells",
    show=False,
    ax=ax,
    palette=dict_cell_color
)
plt.savefig(PLOT_OUTPUT_DIR / 'Annotation_umap_top_gene.png', dpi=600)
plt.show()
