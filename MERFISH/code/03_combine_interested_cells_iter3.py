# ---------------------------------------------------------------------
# Import necessary libraries
# ---------------------------------------------------------------------
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq
import decoupler as dc
import monkeybread as mb

# ---------------------------------------------------------------------
# Global parameters and file paths
# ---------------------------------------------------------------------
# Base path for input and output data
PRI_PATH = Path('../data')

# Output directory for plots and figures
OUTPUT_PATH_PLOT = PRI_PATH / '08_combined_cells_cluster_plot_iter3_umap2'
# Create the plot directory if it doesn't exist
OUTPUT_PATH_PLOT.mkdir(exist_ok=True, parents=True)

# Path for combined cell data (will be created by the script)
OUTPUT_PATH_DATA = PRI_PATH / '04_combined_cells_iter1'
OUTPUT_PATH_DATA.mkdir(exist_ok=True, parents=True)

# File names for input and output data
LAST_CELL_META_WITH_CLUSTER = OUTPUT_PATH_DATA / 'cell_metadata_combine_with_cluster_iter2.csv'
LAST_CELL_META_NO_CLUSTER = OUTPUT_PATH_DATA / 'cell_metadata_combine_iter2.csv'
CELL_BY_GENE_COMBINE_PATH = Path('/lab-share/Cardio-Chen-e2/Public/Yunxia/Hakan_SP/result/cell_seg_DAPI_mouse_rabbit/01_combined_cells/cell_by_gene_combine.csv')

# Output filenames for the new iteration
NEW_CELL_META = OUTPUT_PATH_DATA / 'cell_metadata_combine_iter3.csv'
NEW_CELL_BY_GENE = OUTPUT_PATH_DATA / 'cell_by_gene_combine_iter3.csv'
NEW_CELL_META_WITH_CLUSTER = OUTPUT_PATH_DATA / 'cell_metadata_combine_with_cluster_iter3.csv'

# Cell type annotation panel
GENE_PANEL_PATH = Path('/lab-share/Cardio-Chen-e2/Public/Yunxia/Hakan_SP/Data/marker_data/makers_annotated_merfish.csv')
GENE_PANEL_18_PATH = Path('/lab-share/Cardio-Chen-e2/Public/Yunxia/Hakan_SP/Data/marker_data/makers_annotated_merfish_18.csv')

# Helper function to generate distinct colors for plotting
def generate_distinct_colors(num_colors):
    """
    Generates a list of distinct colors from the 'Set1' colormap.
    """
    colors = plt.cm.get_cmap('Set1', num_colors)
    return [plt.cm.colors.to_hex(colors(i)) for i in range(num_colors)]

def filter_and_prepare_data(cluster_to_remove):
    """
    Loads, filters, and prepares the cell metadata and gene expression data.

    Args:
        cluster_to_remove (int or list): The cluster(s) to be removed from the analysis.

    Returns:
        tuple: A tuple containing the filtered cell metadata and cell-by-gene dataframes.
    """
    print("Step 1: Filtering and preparing data.")
    
    # Load cell metadata with clusters and filter
    cellmeta_with_cluster = pd.read_csv(LAST_CELL_META_WITH_CLUSTER)
    print(f"Initial cell count: {len(cellmeta_with_cluster)}")
    
    cellmeta_filtered = cellmeta_with_cluster[~cellmeta_with_cluster['clusters'].isin([cluster_to_remove])]
    print(f"Cell count after removing cluster {cluster_to_remove}: {len(cellmeta_filtered)}")

    # Load the combined cell metadata and filter based on the remaining cells
    cellmeta_combine_all = pd.read_csv(LAST_CELL_META_NO_CLUSTER)
    selected_entity_ids = cellmeta_filtered['EntityID'].tolist()
    cellmeta_combine_new = cellmeta_combine_all[cellmeta_combine_all['EntityID'].isin(selected_entity_ids)]
    cellmeta_combine_new.to_csv(NEW_CELL_META, index=False)
    print(f"Filtered metadata saved to {NEW_CELL_META}")

    # Load cell-by-gene data and filter based on selected cells
    cell_by_gene_combine = pd.read_csv(CELL_BY_GENE_COMBINE_PATH)
    cell_by_gene_new = cell_by_gene_combine[cell_by_gene_combine['cell'].isin(selected_entity_ids)]
    cell_by_gene_new.to_csv(NEW_CELL_BY_GENE, index=False)
    print(f"Filtered cell-by-gene data saved to {NEW_CELL_BY_GENE}")
    
    # Copy necessary image files for squidpy visualization
    src_img_path = Path('/lab-share/Cardio-Chen-e2/Public/Yunxia/Hakan_SP/result/cell_seg_all/images/micron_to_mosaic_pixel_transform.csv')
    dest_img_path = OUTPUT_PATH_DATA / 'images'
    dest_img_path.mkdir(exist_ok=True, parents=True)
    shutil.copy(src_img_path, dest_img_path)
    print(f"Copied image transformation file to {dest_img_path}")
    
    return cellmeta_combine_new, cell_by_gene_new

def create_and_process_adata(counts_file, meta_file):
    """
    Reads data into an AnnData object and performs standard preprocessing.
    """
    print(f"\nStep 2: Creating and preprocessing AnnData object from {counts_file}.")
    # Read data using squidpy's vizgen reader
    adata = sq.read.vizgen(
        path=OUTPUT_PATH_DATA,
        counts_file=counts_file.name,
        meta_file=meta_file.name,
        transformation_file="micron_to_mosaic_pixel_transform.csv",
    )
    
    # Set TF_ENABLE_ONEDNN_OPTS to 0 to suppress a warning
    os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'

    # Perform QC and filtering
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    print(f"Before cell filtering, adata shape: {adata.X.shape}")
    sc.pp.filter_cells(adata, min_counts=3, inplace=True)
    sc.pp.filter_genes(adata, min_cells=1, inplace=True)
    print(f"After cell filtering, adata shape: {adata.X.shape}")

    # Preprocess the data
    adata.layers["counts"] = adata.X.copy()
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=110)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    
    return adata

def perform_clustering_and_visualization(adata, resolution, output_prefix):
    """
    Performs Leiden clustering and visualizes the results.
    """
    print("\nStep 3: Performing Leiden clustering and visualization.")
    
    # Leiden clustering
    sc.tl.leiden(
        adata,
        key_added="clusters",
        resolution=resolution,
        n_iterations=2,
        directed=False,
    )
    
    num_clusters = len(adata.obs['clusters'].unique())
    print(f"Found {num_clusters} clusters.")
    
    # Generate and assign custom colors for plotting
    custom_colors_list = generate_distinct_colors(num_clusters)
    adata.uns['clusters_colors'] = custom_colors_list

    # UMAP plot
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color='clusters', legend_fontsize=10, legend_fontoutline=2,
               title=f"Clustering of pooled cells ({output_prefix})", show=False, ax=ax)
    plt.savefig(OUTPUT_PATH_PLOT / f'cluster_umap_pooled_cells_{output_prefix}.png', dpi=600)
    
    # Spatial plot
    sq.pl.spatial_scatter(
        adata,
        shape=None,
        color="clusters",
        wspace=0.4,
        size=8,
        return_ax=True,
    )
    plt.savefig(OUTPUT_PATH_PLOT / f'spatial_scatter_pooled_cells_{output_prefix}.png', dpi=600)
    plt.show()

    return adata

def find_and_plot_marker_genes(adata, gene_panel_path, output_prefix):
    """
    Identifies marker genes for each cluster and generates a dot plot.
    """
    print("\nStep 4: Finding and plotting marker genes.")
    
    # Find marker genes for clusters
    sc.tl.rank_genes_groups(adata, groupby="clusters", method="wilcoxon")

    # Load gene panel for cell type annotation
    df_ref_panel = pd.read_csv(gene_panel_path).iloc[1:, :2]
    df_ref_panel.columns = ['Markers', 'cell_type']
    marker_cell = dict(zip(df_ref_panel['Markers'], df_ref_panel['cell_type']))

    # Prepare dictionary for dotplot
    marker_genes_dict = {}
    for cell_type in set(df_ref_panel['cell_type'].tolist()):
        sele_markers = [m for m in df_ref_panel['Markers'].tolist() if marker_cell.get(m) == cell_type]
        if sele_markers:
            marker_genes_dict[cell_type] = sele_markers

    # Plot marker expression dot plot
    fig, ax = plt.subplots(figsize=(50, 8))
    sc.pl.dotplot(adata, marker_genes_dict, "clusters", dendrogram=True, show=False, ax=ax)
    plt.savefig(OUTPUT_PATH_PLOT / f'marker_expression_{output_prefix}.png', dpi=300, bbox_inches='tight')
    plt.show()

def perform_decoupler_annotation(adata):
    """
    Performs machine learning-based cell type annotation using decoupler.
    """
    print("\nStep 5: Performing decoupler-based cell type annotation.")
    
    # Load gene panel
    df_ref_panel = pd.read_csv(GENE_PANEL_PATH).iloc[1:, :2]
    df_ref_panel.columns = ['Markers', 'cell_type']

    # Run ORA (Over Representation Analysis) with decoupler
    dc.run_ora(
        mat=adata,
        net=df_ref_panel,
        source='cell_type',
        target='Markers',
        min_n=3,
        verbose=True,
        use_raw=False
    )
    
    acts = dc.get_acts(adata, obsm_key='ora_estimate')
    
    # Handle infinite values
    acts_v = acts.X.ravel()
    max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
    acts.X[~np.isfinite(acts.X)] = max_e
    
    # Rank sources (cell types) by cluster
    df = dc.rank_sources_groups(acts, groupby='clusters', reference='rest', method='t-test_overestim_var')
    n_ctypes = 1
    ctypes_dict = df.groupby('group').head(n_ctypes).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
    ctypes_dict_new = {k: v[0] for k, v in ctypes_dict.items()}

    # Map clusters to annotated cell types
    acts.obs["cell_type"] = acts.obs["clusters"].map(ctypes_dict_new).astype("category")

    # Get colors from the original adata object
    cluster_colors = {cluster: color for cluster, color in zip(adata.obs['clusters'].unique(), adata.uns['clusters_colors'])}
    
    # Map colors to the new cell types
    dict_cell_color = {
        ctypes_dict_new[cluster]: cluster_colors[cluster]
        for cluster in ctypes_dict_new
    }
    
    acts.uns['cell_type_colors'] = [dict_cell_color[ct] for ct in acts.obs['cell_type'].cat.categories]
    
    # Plot UMAP with annotated cell types
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        acts,
        color="cell_type",
        legend_fontsize=10,
        title="Annotation of cells",
        show=False,
        ax=ax,
    )
    plt.savefig(OUTPUT_PATH_PLOT / 'Annotation_umap_decoupler.png', dpi=600)
    acts.obs.to_csv(OUTPUT_PATH_PLOT / 'culster_obs_allgenes.csv')
    plt.show()

def main():
    """
    Main function to run the entire analysis workflow.
    """
    # ------------------------------------------------
    # Main workflow - Initial analysis
    # ------------------------------------------------
    print("Starting main analysis workflow...")
    # Filter and prepare data by removing cluster 5
    cellmeta_combined, cell_by_gene_combined = filter_and_prepare_data(cluster_to_remove=5)
    
    # Create and preprocess the AnnData object
    adata = create_and_process_adata(counts_file=NEW_CELL_BY_GENE, meta_file=NEW_CELL_META)
    
    # Perform initial clustering and visualization (resolution=0.5 as in your original script)
    adata = perform_clustering_and_visualization(adata, resolution=0.5, output_prefix='main')
    
    # Save the adata.obs DataFrame which contains the clusters
    adata.obs.to_csv(NEW_CELL_META_WITH_CLUSTER)

    # Find and plot marker genes using the main gene panel
    find_and_plot_marker_genes(adata, GENE_PANEL_PATH, 'all_clusters')
    
    # Perform decoupler annotation
    perform_decoupler_annotation(adata)

    # ------------------------------------------------
    # Re-clustering workflow for a specific cluster (LPM-1b, which was cluster 2)
    # ------------------------------------------------
    print("\nStarting re-clustering workflow for cluster 2...")
    
    # Filter for cells in cluster 2 (LPM-1b)
    adata_cluster = adata.obs
    adata_cluster_2 = adata_cluster[adata_cluster['clusters'] == '2']
    
    selected_entity_ids_cluster2 = adata_cluster_2.index.astype(int).tolist()
    
    # Filter the original dataframes to get only cluster 2 cells
    cellmeta_cluster2 = cellmeta_combined[cellmeta_combined['EntityID'].isin(selected_entity_ids_cluster2)]
    cell_by_gene_cluster2 = cell_by_gene_combined[cell_by_gene_combined['cell'].isin(selected_entity_ids_cluster2)]
    
    # Define new filenames for the reclustering data
    new_cell_meta_cluster = OUTPUT_PATH_DATA / 'cell_metadata_combine_iter3_cluster2.csv'
    new_cell_by_gene_cluster = OUTPUT_PATH_DATA / 'cell_by_gene_combine_iter3_cluster2.csv'
    
    # Save the reclustering data
    cellmeta_cluster2.to_csv(new_cell_meta_cluster, index=False)
    cell_by_gene_cluster2.to_csv(new_cell_by_gene_cluster, index=False)
    
    # Create and process new adata object for reclustering
    adata_recluster = create_and_process_adata(counts_file=new_cell_by_gene_cluster, meta_file=new_cell_meta_cluster)
    
    # Perform new clustering and visualization for LPM-1b
    adata_recluster = perform_clustering_and_visualization(adata_recluster, resolution=0.5, output_prefix='LPM-1b_recluster')

    # Find and plot marker genes for the re-clustered data
    find_and_plot_marker_genes(adata_recluster, GENE_PANEL_18_PATH, 'LPM-1b_recluster')

if __name__ == "__main__":
    main()
