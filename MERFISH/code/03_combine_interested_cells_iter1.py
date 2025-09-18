# %%
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq
import os
import subprocess
import decoupler as dc
import monkeybread as mb

# --- User Configuration ---
# Set the base directory for all data. This path should point to the root of
# your project's data, which contains folders like 'slide1region5' and
# '04_combined_cells_iter1'.
# Example: If your script is at '~/Projects/my_project/scripts/run.py' and your
# data is at '~/Projects/my_project/data/cell_seg_DAPI_mouse_rabbit/',
# you would set the path accordingly.
# Note: The `Path(__file__).resolve().parent` gets the directory of this script.
# The following assumes your data folder is a sibling of the script folder.
# Adjust the path as needed for your specific file structure.
base_data_dir = Path('../data')

# --- Functions ---
def deal_str(data):
    """Appends a tab character to a string representation of data."""
    return str(data) + '\t'

# --- Main Script ---
if __name__ == '__main__':
    filenames = ['slide1region5','slide1region7','slide1region9','slide2region2','slide2region5','slide1region4','slide1region6','slide1region8','slide2region0','slide2region4']
    output_path = base_data_dir / '04_combined_cells_iter1'
    
    dict_data = {
        "slide1region4": [11],
        "slide1region5": [11,5],
        "slide1region6": [11],
        "slide1region7": [11],
        "slide1region8": [11,5],
        "slide1region9": [11],
        "slide2region0": [11],
        "slide2region2": [11],
        "slide2region4": [11,5],
        "slide2region5": [11,5]  
    }

    dict_cells = {}
    datameta_combined = None

    # %% Combine selected cell data from different regions
    for index, filename in enumerate(filenames):
        print('===========================================')
        print(filename)
        
        vizgen_dir = base_data_dir / filename
        datapath_meta_pooled = vizgen_dir / 'cell_metadata_pooled.csv'
        
        if not datapath_meta_pooled.exists():
            print(f"Error: File not found at {datapath_meta_pooled}")
            continue

        datameta_pooled = pd.read_csv(datapath_meta_pooled, index_col=0)
        datameta_selected = datameta_pooled[datameta_pooled['cluster'].isin(dict_data[filename])]
        dict_cells[filename] = datameta_selected.shape[0]
        
        if datameta_combined is None:
            datameta_combined = datameta_selected
        else:
            datameta_combined = pd.concat([datameta_combined, datameta_selected], axis=0)

    if datameta_combined is not None:
        datameta_combined = datameta_combined.reset_index(drop=True)
        cell_meta = datameta_combined.iloc[:, 2:-2]     
        os.makedirs(output_path, exist_ok=True)
        cell_meta.to_csv(output_path / 'cell_metadata_combine.csv', index=False)  
    else:
        print("No combined data to process. Exiting.")
        exit()

    # %% Display combined cell counts and metadata
    print(datameta_combined)
    dict_cells_index = {'cell_number': dict_cells}
    print(pd.DataFrame(dict_cells_index))

    # %% Filter combined cell-by-gene data
    cell_by_gene_combine_path = base_data_dir / '01_combined_cells' / 'cell_by_gene_combine.csv'
    if not cell_by_gene_combine_path.exists():
        print(f"Error: File not found at {cell_by_gene_combine_path}")
        exit()

    cell_by_gene_combine_data = pd.read_csv(cell_by_gene_combine_path)
    selected_cell_lists = datameta_combined['EntityID'].tolist()

    cell_by_gene_com_new = cell_by_gene_combine_data[cell_by_gene_combine_data['cell'].isin(selected_cell_lists)]
    print(cell_by_gene_com_new)

    # %% Cell type annotation
    gene_panel = base_data_dir / 'marker_data' / 'makers_annotated_merfish.csv'
    if not gene_panel.exists():
        print(f"Error: File not found at {gene_panel}")
        exit()
    
    df_ref_panel_ini = pd.read_csv(gene_panel)
    df_ref_panel_all = df_ref_panel_ini.iloc[1:, :2]
    df_ref_panel_all.columns = ['Markers', 'cell_type']
    print(df_ref_panel_all)

    marker_cell = dict(zip(df_ref_panel_all['Markers'], df_ref_panel_all['cell_type']))

    # %% Filter gene panel to 127 markers
    list_127 = ['cell'] + df_ref_panel_all['Markers'].tolist()
    cell_by_gene_com_new_127 = cell_by_gene_com_new[list_127]
    print(cell_by_gene_com_new_127)
    
    # Save the filtered cell-by-gene data
    cell_by_gene_com_new_127.to_csv(output_path / 'cell_by_gene_combine.csv', index=False)

    # %% Use subprocess to create directory and copy file
    # Note: Using `subprocess` is not the most Pythonic way, but since it was in your
    # original code, we are simply modifying the paths. A better approach would be
    # to use `os.makedirs` and `shutil.copyfile`.
    images_dir = output_path / 'images'
    os.makedirs(images_dir, exist_ok=True)
    
    source_transform_file = base_data_dir / '..' / 'cell_seg_all' / 'images' / 'micron_to_mosaic_pixel_transform.csv'
    dest_transform_file = images_dir / 'micron_to_mosaic_pixel_transform.csv'
    
    if source_transform_file.exists():
        try:
            subprocess.run(f'cp "{source_transform_file}" "{dest_transform_file}"', shell=True, check=True)
            print(f"Successfully copied {source_transform_file} to {dest_transform_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error copying file: {e}")
            print(f"Command output: {e.stdout}")
            print(f"Command error: {e.stderr}")
    else:
        print(f"Source file not found: {source_transform_file}")

    # %% Scanpy and Squidpy analysis
    adata = sq.read.vizgen(
        path=str(output_path),
        counts_file="cell_by_gene_combine.csv",
        meta_file="cell_metadata_combine.csv",
        transformation_file=str(Path('images') / 'micron_to_mosaic_pixel_transform.csv'),
    )

    sc.pp.calculate_qc_metrics(adata, percent_top=(10, 50, 100, 120), inplace=True)
    print("before cell filter: {}".format(adata.X.shape))

    sc.pp.filter_cells(adata, min_counts=3, inplace=True)
    sc.pp.filter_genes(adata, min_cells=1, inplace=True)
    print("after cell filter: {}".format(adata.X.shape))
    
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
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

    len_cluster = len(list(set(adata.obs['clusters'].tolist())))

    def generate_distinct_colors(num_colors):
        colors = plt.cm.get_cmap('tab20', num_colors)
        return [colors(i) for i in range(num_colors)]

    def rgba_to_hex(rgba):
        return '#{:02x}{:02x}{:02x}'.format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))

    num_colors = len_cluster
    distinct_colors = generate_distinct_colors(num_colors)
    hex_colors = [rgba_to_hex(color) for color in distinct_colors]
    print(hex_colors)
    
    plt.figure(figsize=(10, 2))
    for i in range(num_colors):
        plt.plot(i, 0, marker='o', markersize=20, color=hex_colors[i])
    plt.xlim(-1, num_colors)
    plt.yticks([])
    plt.title(f'{len_cluster} Distinct Colors from tab20 Colormap')
    plt.show()

    custom_colors = dict(zip(list(range(len_cluster)), hex_colors))
    adata.uns['clusters_colors'] = list(custom_colors.values())

    output_path_plot = base_data_dir / '05_combined_cells_cluster_plot_iter1_umap2'
    os.makedirs(output_path_plot, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata,
        color=['clusters'],            
        legend_fontsize=10,
        legend_fontoutline=2,
        title="clustering of pooled cells",
        show=False,
        ax=ax,
    )
    plt.savefig(output_path_plot / 'cluster_umap_pooled_cells.png', dpi=600)

    sq.pl.spatial_scatter(
        adata,
        shape=None,
        color=["clusters"],
        wspace=0.4,
        size=8,
        return_ax=True,
    )
    plt.savefig(output_path_plot / 'Visualize annotation spatial coordinates.png', dpi=600)

    # %% Cell metadata and marker gene analysis
    cell_meta_with_cluster = pd.DataFrame(adata.obs)
    cell_meta_with_cluster['EntityID'] = cell_meta_with_cluster.index.map(deal_str)
    cell_meta_with_cluster.to_csv(output_path / 'cell_metadata_combine_with_cluster_iter1.csv')

    print(f"Shape of cell_meta_with_cluster: {cell_meta_with_cluster.shape}")

    marker_genes_dict = {}
    for cell_type in set(df_ref_panel_all['cell_type'].tolist()):
        sele_markers = [marker for marker in set(df_ref_panel_all['Markers'].tolist()) if marker_cell[marker] == cell_type]
        marker_genes_dict[cell_type] = sele_markers

    fig, ax = plt.subplots(figsize=(50, 8))
    sc.pl.dotplot(adata, marker_genes_dict, "clusters", dendrogram=True, show=False, ax=ax)
    plt.savefig(output_path_plot / 'marker_expression_top_gene.png', dpi=300)

    # %% Machine learning annotation
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
    acts_v = acts.X.ravel()
    max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
    acts.X[~np.isfinite(acts.X)] = max_e

    df = dc.rank_sources_groups(acts, groupby='clusters', reference='rest', method='t-test_overestim_var')
    n_ctypes = 3
    ctypes_dict = df.groupby('group').head(n_ctypes).groupby('group')['names'].apply(list).to_dict()

    ctypes_dict_new = {x: ctypes_dict[x][0] for x in ctypes_dict.keys()}
    print(ctypes_dict_new)

    n_top_genes = 1187
    dict_cell_color = {
        value: custom_colors[int(key)]
        for key, value in ctypes_dict_new.items()
    }
    print(dict_cell_color)

    acts.obs["cell type"] = acts.obs["clusters"].map(ctypes_dict_new).astype("category")
    acts.uns['cell type_colors'] = list(dict_cell_color.values())

    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        acts,
        color="cell type", 
        legend_fontsize=10,
        legend_fontoutline=2,
        title=f"Annotation of cells_{n_top_genes}",
        show=False,
        ax=ax,
        palette=dict_cell_color,   
    )
    plt.savefig(output_path_plot / f'Annotation_umap_top_gene_{n_top_genes}.png', dpi=600)
    
    acts.obs.to_csv(base_data_dir / 'culster_obs_allgenes.csv')
    plt.show()
