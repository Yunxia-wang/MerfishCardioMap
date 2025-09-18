import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
import scanpy as sc
import squidpy as sq
import seaborn as sns
import decoupler as dc
import monkeybread as mb
# ---
# NOTE: This script assumes a project structure where the script is
# located in a subdirectory (e.g., 'scripts/') and the data is
# located in a sibling directory at the project root (e.g., 'data/').
# Adjust the `base_data_dir` path below if your structure is different.
# The path is constructed to point to the 'cell_seg_DAPI_mouse_rabbit' folder.
# ---
base_data_dir = Path(__file__).resolve().parent.parent / 'data' / 'cell_seg_DAPI_mouse_rabbit'


# Define the paths for the cluster data and plotting results
cluster_datapath = base_data_dir / '01_combined_cells'
plot_result = base_data_dir / '03_plot_spatial_for_cluster'
os.makedirs(plot_result, exist_ok=True)

# Get the cell meta with cluster
annotated_cell = pd.read_csv(cluster_datapath / 'cell_metadata_combine_with_cluster.csv')
dict_cell_clu = annotated_cell.set_index('EntityID')['clusters'].to_dict()

len_cluster = len(list(set(annotated_cell['clusters'].tolist())))
cell_list = [int(i) for i in annotated_cell['EntityID']]

def generate_distinct_colors(num_colors):
    colors = plt.cm.get_cmap('tab20', num_colors)
    return [colors(i) for i in range(num_colors)]

def rgba_to_hex(rgba):
    return '#{:02x}{:02x}{:02x}'.format(int(rgba[0]*255), int(rgba[1]*255), int(rgba[2]*255))

num_colors = len_cluster
distinct_colors = generate_distinct_colors(num_colors)
hex_colors = [rgba_to_hex(color) for color in distinct_colors]
print(hex_colors)

custom_colors = dict(zip(list(range(len_cluster)), hex_colors))
custom_colors[num_colors] = '#f1f2f1'
print(custom_colors)

# Plot the colors to visualize them
plt.figure(figsize=(12, 2))
for i in range(len(custom_colors)):
    plt.plot(i, 0, marker='o', markersize=20, color=custom_colors[i])
plt.xlim(-1, len(custom_colors))
plt.yticks([])
plt.title(str(len_cluster) + ' Distinct Colors from tab20 Colormap')
plt.show()

filenames = ['slide1region1','slide1region5','slide1region7','slide1region9','slide2region2','slide2region5','slide1region4','slide1region6','slide1region8','slide2region0','slide2region4']

for filename in filenames:
    print('===========================================')
    print(filename)

    vizgen_dir = base_data_dir / filename
    resultpath = plot_result / filename
    os.makedirs(resultpath, exist_ok=True)
    cellpath = "cell_by_gene.csv"
    cellmetapath = "cell_metadata_pooled.csv"

    # generate new cell_metadata file
    datapath_meta = vizgen_dir / 'cell_metadata.csv'
    datameta = pd.read_csv(datapath_meta)
    cell_clusters = []
    
    count_cell = 0
    for i in datameta['EntityID'].tolist():
        if i in cell_list:
            count_cell += 1
            cell_clusters.append(dict_cell_clu[i])
        else:
            cell_clusters.append(num_colors)
    
    datameta.loc[:, ('cluster')] = cell_clusters
    
    datameta.to_csv(vizgen_dir / 'cell_metadata_pooled.csv', index=False)
    print('count_cell:{}'.format(count_cell))
    print('datameta.shape:{}'.format(datameta.shape))
    print(datameta.head())

    adata = sq.read.vizgen(
        path=str(vizgen_dir),
        counts_file=cellpath,
        meta_file=cellmetapath,
        transformation_file="micron_to_mosaic_pixel_transform.csv",
    )

    # QC   
    sc.pp.calculate_qc_metrics(adata, percent_top=(10, 50, 100, 130), inplace=True)
    fdr = adata.obsm["blank_genes"].to_numpy().sum() / adata.var["total_counts"].sum() * 100
    print("fdr:{}".format(fdr))

    # filter   
    sc.pp.filter_cells(adata, min_counts=3, inplace=True)
    sc.pp.filter_genes(adata, min_cells=1, inplace=True)
    print("after cell filter: {}".format(adata.X.shape))
    
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    
    # make the cluter as cell type
    cell_types = [str(i) for i in adata.obs['cluster'].tolist()]
    # Create an AnnData object
    adata.obs['cell_type'] = pd.Categorical(cell_types)
    
    custom_colors_region = {}
    for cell_type in cell_types:
        custom_colors_region[cell_type] = custom_colors[int(cell_type)]
    
    # Add the custom color palette to the AnnData object
    adata.uns['cell_type_colors'] = list(custom_colors_region.values())    

    # plot cell spatial figure for all clusters of each region
    fig, axs = plt.subplots(1, 1, figsize=(16, 20))
    sc.pl.embedding(
        adata,
        "spatial",
        color='cell_type',
        groups=cell_types,
        wspace=0.4,
        size=60,
        ax=axs,
        palette=custom_colors_region,
        show=False
    )
    axs.yaxis.set_inverted(True)
    plt.savefig(resultpath / f'{filename}_all_cluster_kernel_density_spatial.png', dpi=600)
    
    # plot cell spatial figure for each clusters of each region     
    plot_cluster = set(cell_types)
    plot_cluster.remove(str(num_colors))
    for i in list(plot_cluster):
        fig, axs = plt.subplots(1, 1, figsize=(16, 20))
        sc.pl.embedding(
            adata,
            "spatial",
            color='cell_type',
            groups=str(i),
            wspace=0.4,
            ax=axs,
            size=60,
            palette=custom_colors_region,
            show=False,
            na_color='#f1f2f1'
        )
        axs.yaxis.set_inverted(True)
        plt.savefig(resultpath / f'{filename}_cluster{i}_kernel_density_spatial.png', dpi=600)
    
# ---
# Original `deal_str` function and processing loop
# ---
def deal_str(data):
    data = str(data) + '\t'
    return data

for filename in filenames:
    print('===========================================')
    print(filename)
    
    vizgen_dir = base_data_dir / filename / 'cell_metadata_pooled.csv'
    cell_meta = pd.read_csv(vizgen_dir)
    print(cell_meta)
    cell_meta['Cell_ID'] = cell_meta['EntityID'].map(deal_str)
    print(cell_meta)
    cell_meta.to_csv(vizgen_dir, index=False)
