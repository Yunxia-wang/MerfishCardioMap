# --- Imports and Setup ---
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import squidpy as sq
import pandas as pd
import numpy as np
from paste2 import PASTE2, projection
import matplotlib.colors as mcolors
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import plotly.express as px
import plotly.io as pio
import plotly.offline as py
import seaborn as sns
import pandas as pd
import geopandas as gpd
import plotly.graph_objects as go
from shapely.geometry import MultiPolygon, Polygon
from shapely import wkb
from shapely.affinity import translate
from skimage.draw import disk
import importlib.metadata
import os
from pyevtk.hl import pointsToVTK, imageToVTK
import scipy.ndimage as ndi
import vtk
from vtk.util import numpy_support

print(f"PASTE2 version: {importlib.metadata.version('paste2')}")
print(f"PyEVTK version: {importlib.metadata.version('pyevtk')}")

# --- Set up relative paths for data input and output ---
# Assuming your data is in a sub-folder named 'data' relative to this script.
data_base_path = '../data'
output_base_path = '../result'

# Make sure the output directory exists
output_dir = os.path.join(output_base_path, 'STIM_3D')
os.makedirs(output_dir, exist_ok=True)
savepath = output_dir

# --- Read in Spatial Transcriptomics slices as AnnData objects ---
layer_to_color_map = {'cluster_0': '#00aa00', 'cluster_1': '#55007f', 'cluster_2': '#ffaa00', 'cluster_3': '#00aaff', 'cluster_4': '#ff557f','cluster_5': '#cccccc'}

dict_region_names = {
    '17_slide1region4':'3D_embryo4',
    '16_slide2region5':'3D_embryo2',
    '15_slide2region4':'3D_embryo2',
    '14_slide2region5':'3D_embryo4',
    '13_slide1region5':'3D_embryo4',
    '12_slide2region4':'3D_embryo4',
    '11_slide1region8':'3D_embryo4',
    '10_slide1region5':'3D_embryo2',
    '9_slide1region8':'3D_embryo2',
    '8_slide1region7':'3D_embryo2',
    '7_slide1region9':'3D_embryo2',
    '6_slide2region0':'3D_embryo1_right',
    '5_slide2region0':'3D_embryo1_left',
}

# Create a dictionary to store the AnnData objects
anndata_dict = {}

for index, region_name0 in enumerate(dict_region_names.keys()):
    print('----------------')
    print(region_name0)
    
    region_name = list(dict_region_names.keys())[index].split('_')[1]
    # Use os.path.join to build the relative path
    output_path_name = os.path.join(data_base_path, region_name)
    counts_file_name = 'cell_by_gene_'+ dict_region_names[region_name0] +'.csv'
    meta_file_name = 'cell_metadata_' + dict_region_names[region_name0] +'.csv'
    
    sliceX = sq.read.vizgen(
        path=output_path_name,
        counts_file=counts_file_name,
        meta_file=meta_file_name,
        transformation_file="micron_to_mosaic_pixel_transform.csv",
        )
    
    cell_types = [str(i) for i in sliceX.obs['cluster'].tolist()]
    
    custom_colors_region = {}
    for cell_type in cell_types:
        custom_colors_region[cell_type] = layer_to_color_map[cell_type]
    
    sliceX.uns['cluster_colors'] = list(custom_colors_region.values())
    sliceX.obs['cluster'] = sliceX.obs['cluster'].astype('category')
    
    anndata_dict[region_name0] = sliceX
    print(sliceX.obs.shape)
slice_lists = list(anndata_dict.values())    

# --- Pairwise alignment ---
pi_AB_dict = {}

for index, i in enumerate(range(len(slice_lists))):
    sliceA = slice_lists[i]
    sliceB = slice_lists[i+1]
    pi_AB = PASTE2.partial_pairwise_align(sliceA, sliceB, s=0.3)
    pi_AB_dict[i] = pi_AB
    
    if i == len(slice_lists) - 2:
        break

# --- Plotting the results of pairwise alignment ---
def largest_indices(ary, n):
    """Returns the n largest indices from a numpy array."""
    flat = ary.flatten()
    indices = np.argpartition(flat, -n)[-n:]
    indices = indices[np.argsort(-flat[indices])]
    return np.unravel_index(indices, ary.shape)


def plot2D_samples_mat(xs, xt, G, thr=1e-8, alpha=0.2, top=1000, weight_alpha=False, **kwargs):
    if ('color' not in kwargs) and ('c' not in kwargs):
        kwargs['color'] = '#cccccc'
    mx = G.max()
    idx = largest_indices(G, top)
    for l in range(len(idx[0])):
        plt.plot([xs[idx[0][l], 0], xt[idx[1][l], 0]], [xs[idx[0][l], 1], xt[idx[1][l], 1]],
                 alpha=alpha * (1 - weight_alpha) + (weight_alpha * G[idx[0][l], idx[1][l]] / mx), c='k')


def plot_slice_pairwise_alignment(slice1, slice2, pi, thr=1 - 1e-8, alpha=0.05, top=1000, name='',
                                  weight_alpha=False):
    coordinates1, coordinates2 = slice1.obsm['spatial'], slice2.obsm['spatial']
    offset = (coordinates1[:, 0].max() - coordinates2[:, 0].min()) * 1.1
    temp = np.zeros(coordinates2.shape)
    temp[:, 0] = offset
    plt.figure(figsize=(20, 10))
    plot2D_samples_mat(coordinates1, coordinates2 + temp, pi, thr=thr, c='#cccccc', alpha=alpha, top=top,
                       weight_alpha=weight_alpha)
    plt.scatter(coordinates1[:, 0], coordinates1[:, 1], linewidth=0, s=100, marker=".", color=list(
        slice1.obs['cluster'].map(
            dict(zip(slice1.obs['cluster'].cat.categories, slice1.uns['cluster_colors'])))))
    plt.scatter(coordinates2[:, 0] + offset, coordinates2[:, 1], linewidth=0, s=100, marker=".", color=list(
        slice2.obs['cluster'].map(
            dict(zip(slice2.obs['cluster'].cat.categories, slice2.uns['cluster_colors'])))))
    plt.gca().invert_yaxis()
    plt.axis('off')
    plt.show()

for i in range(len(pi_AB_dict)):
    print(f"Plotting alignment for slice {i} and {i+1}")
    plot_slice_pairwise_alignment(slice_lists[i], slice_lists[i+1], pi_AB_dict[i])


# --- Project all slices onto the same coordinate system ---
pis = list(pi_AB_dict.values())
slices = slice_lists

new_slices = projection.partial_stack_slices_pairwise(slices, pis)

def plot_slices_overlap(slices, layer_to_color_map=layer_to_color_map):
    plt.figure(figsize=(15,15))
    for i in range(len(slices)):
        adata = slices[i]
        colors = list(adata.obs['cluster'].astype('str').map(layer_to_color_map))
        plt.scatter(adata.obsm['spatial'][:,0],adata.obsm['spatial'][:,1],linewidth=0,s=100, marker=".",color=colors)
    plt.legend(handles=[mpatches.Patch(color=layer_to_color_map[adata.obs['cluster'].cat.categories[i]], label=adata.obs['cluster'].cat.categories[i]) for i in range(len(adata.obs['cluster'].cat.categories))],fontsize=10,title='Cortex layer',title_fontsize=15,bbox_to_anchor=(1, 1))
    plt.gca().invert_yaxis()
    plt.axis('off')
    plt.show()
    
plot_slices_overlap(new_slices[0:len(new_slices)])

# --- Export to VTK for 3D visualization ---
# Export as a point cloud
all_coords = []
all_clusters = []

for z, adata in enumerate(new_slices):
    coords = adata.obsm["spatial"]
    n_cells = coords.shape[0]
    z_coords = np.full((n_cells,1), z*20.0)
    coords_3d = np.hstack([coords, z_coords])
    all_coords.append(coords_3d)
    
    clusters = np.array([int(c.split("_")[1]) for c in adata.obs["cluster"]])
    all_clusters.append(clusters)

all_coords = np.vstack(all_coords)
all_clusters = np.concatenate(all_clusters)

x = np.ascontiguousarray(all_coords[:,0])
y = np.ascontiguousarray(all_coords[:,1])
z = np.ascontiguousarray(all_coords[:,2])
clusters = np.ascontiguousarray(all_clusters.astype(np.int32))

filename = os.path.join(savepath, "cells_points")
pointsToVTK(filename, x, y, z, data={"cluster": clusters})

print(f"Saved {filename}.vtk (point cloud)")

# Export as a volume image
# This assumes you have the 'cluster_volume' and 'img_shape' variables defined from prior steps.
# The user's original notebook snippet for this part was not included, so I'll create a dummy
# volume for demonstration purposes. In your full script, these would be populated correctly.
print("Generating dummy volume data for VTK export. Please replace this with your actual data.")
img_shape = (500, 500)
cluster_volume = np.zeros((len(new_slices), img_shape[0], img_shape[1]), dtype=np.uint8)
# Example dummy data
for z_idx, adata in enumerate(new_slices):
    coords = adata.obsm["spatial"]
    clusters = [int(c.split("_")[1]) + 1 for c in adata.obs["cluster"]] # +1 to avoid 0
    for coord, cluster in zip(coords, clusters):
        # A very basic rasterization
        rr, cc = disk((coord[1], coord[0]), 5, shape=img_shape)
        cluster_volume[z_idx, rr, cc] = cluster

# Create RGB volume
cluster_colors = np.array([
    [0, 170, 0],     # cluster_0
    [85, 0, 127],    # cluster_1
    [255, 170, 0],   # cluster_2
    [0, 170, 255],   # cluster_3
    [255, 85, 127],  # cluster_4
    [204, 204, 204]  # cluster_5
], dtype=np.uint8)

# The following two blocks were from the user's notebook and were not fully connected to the rest of the code.
# I've left them here for context, assuming you would integrate them into your workflow.
# They are correct for saving a volume to a VTI file.

def save_volume_to_vti(filename, volume, spacing=(1.0, 1.0, 1.0)):
    """
    Save a 3D numpy array as a .vti file using VTK Python bindings.
    """
    if volume.dtype != np.uint8:
        volume = (255 * (volume / volume.max())).astype(np.uint8)

    vtk_data = numpy_support.numpy_to_vtk(num_array=volume.ravel(order='F'),
                                          deep=True,
                                          array_type=vtk.VTK_UNSIGNED_CHAR)

    img = vtk.vtkImageData()
    nz, ny, nx = volume.shape
    img.SetDimensions(nx, ny, nz)
    img.SetSpacing(spacing)
    img.GetPointData().SetScalars(vtk_data)

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(img)
    writer.Write()
    print(f"Saved {filename}")

# Example usage for the dummy volume
save_volume_to_vti(os.path.join(savepath, "aligned_volume.vti"), cluster_volume, spacing=(10.0, 1.0, 1.0))

filename = os.path.join(savepath, "aligned_clusters_30")
spacing = (10.0, 1.0, 1.0)  # adjust ZXY voxel spacing
volume = np.ascontiguousarray(cluster_volume)
imageToVTK(filename, origin=(0,0,0), spacing=spacing, pointData={"cluster": volume})
print(f"Saved {filename}.vti")
