import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def process_embryo_region(data_path_relative, section_name, rotation_params=None, filter_coords=None):
    """
    Processes a single embryo region by filtering cells and saving the results.
    This consolidated function can handle both coordinate filtering and rotation.

    Args:
        data_path_relative (str): The specific relative path to the region's data.
        section_name (str): A descriptive name for the section (e.g., 'slide2region4 Embryo 1').
        rotation_params (tuple, optional): A tuple of (theta_deg, invert_xaxis, invert_yaxis)
                                          for applying rotation and axis inversion.
        filter_coords (tuple, optional): A tuple of (x_min, x_max, y_min, y_max)
                                         to filter cell coordinates.
    """
    print(f"\n--- Processing {section_name} ---")

    # Define the base directory relative to where this script is executed.
    base_path = '../data'
    datapath = os.path.join(base_path, 'cell_seg_DAPI_mouse_rabbit', data_path_relative)
    
    # Load shared metadata once
    cellmetapath_combined_iter3 = os.path.join(base_path, '04_combined_cells_iter1', 'cell_metadata_combine_with_cluster_iter3.csv')
    cellmetapath_combined = os.path.join(base_path, '01_combined_cells', 'cell_metadata_combine_with_cluster.csv')

    try:
        annotated_cell_1032 = pd.read_csv(cellmetapath_combined_iter3)
        cells_list_1032 = annotated_cell_1032['EntityID'].tolist()
        dict_cell_clu_1032 = annotated_cell_1032.set_index('EntityID')['clusters'].to_dict()

        annotated_cell = pd.read_csv(cellmetapath_combined)
        annotated_cell = annotated_cell[(annotated_cell['transcript_count'] > 30) | (annotated_cell['EntityID'].isin(cells_list_1032))]
        dict_cell_clu = annotated_cell.set_index('EntityID')['clusters'].to_dict()
    except FileNotFoundError as e:
        print(f"Error loading common metadata: {e}")
        return

    # Load data for the specific region
    try:
        cell_meta_path = os.path.join(datapath, 'cell_metadata_pooled_iter3.csv')
        df1 = pd.read_csv(cell_meta_path)
    except FileNotFoundError as e:
        print(f"Error loading region data: {e}")
        return

    # Filter by common cells
    df1 = df1[df1['EntityID'].isin(dict_cell_clu.keys())].copy()

    # Apply coordinate filtering if provided
    if filter_coords:
        x_min, x_max, y_min, y_max = filter_coords
        df1 = df1[(df1['center_x'] > x_min) & (df1['center_x'] < x_max) &
                  (df1['center_y'] > y_min) & (df1['center_y'] < y_max)]
        
    # Get cell centroid coordinates
    xI = np.array(df1['center_x'])
    yI = np.array(df1['center_y'])

    # Plot raw data
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(xI, yI, s=1, alpha=1, label='source')
    ax.legend(markerscale=10)
    ax.set_title(f"Raw Cell Coordinates for {section_name}")
    if rotation_params:
        if rotation_params[1]:
            ax.invert_xaxis()
        if rotation_params[2]:
            ax.invert_yaxis()
    plt.show()

    # Apply rotation if specified
    if rotation_params and rotation_params[0] is not None:
        theta_deg = rotation_params[0]
        theta0 = (np.pi / 180) * -theta_deg
        L = np.array([[np.cos(theta0), -np.sin(theta0)],
                      [np.sin(theta0), np.cos(theta0)]])
        source_L = np.matmul(L, np.array([xI, yI]))
        df1['center_x'] = source_L[0]
        df1['center_y'] = source_L[1]
        print(f"Applied a rotation of {theta_deg} degrees.")

    # Get cell cluster of selected cell in specific regions
    cluster_iter3 = []
    for i in df1['EntityID'].tolist():
        if i in cells_list_1032:
            cluster_iter3.append('cluster_' + str(dict_cell_clu_1032[i]))
        else:
            cluster_iter3.append('cluster_5')

    df1['cluster'] = cluster_iter3
    df1['cluster'] = df1['cluster'].astype('category')
    print(f"Filtered cell metadata shape for saving: {df1.shape}")

    # Save results
    file_name_suffix = section_name.replace(' ', '_').replace('-', '').replace('Embryo', '3D_embryo')
    
    output_meta_path = os.path.join(datapath, f'cell_metadata_{file_name_suffix}.csv')
    df1.to_csv(output_meta_path, index=False)
    print(f"Saved cell metadata to: {output_meta_path}")

    cell_by_gene_path = os.path.join(datapath, 'cell_by_gene.csv')
    try:
        cell_by_gene = pd.read_csv(cell_by_gene_path)
    except FileNotFoundError as e:
        print(f"Error loading gene data: {e}")
        return

    cell_by_gene = cell_by_gene[cell_by_gene['cell'].isin(df1['EntityID'].tolist())]
    print(f"Filtered cell by gene shape for saving: {cell_by_gene.shape}")
    
    output_gene_path = os.path.join(datapath, f'cell_by_gene_{file_name_suffix}.csv')
    cell_by_gene.to_csv(output_gene_path, index=False)
    print(f"Saved cell by gene data to: {output_gene_path}")

    print("-" * 30)

if __name__ == '__main__':
    # Define a list of all processing tasks from your original scripts
    processing_tasks = [
        # Original Script 1 tasks
        {'path': 'slide2region4', 'name': 'slide2region4 Embryo 1', 'rot_params': (None, True, False)},
        {'path': 'slide1region8', 'name': 'slide1region8 Embryo 1', 'rot_params': (None, False, True)},
        {'path': 'slide1region7', 'name': 'slide1region7 Embryo 1', 'rot_params': (None, False, False)},
        {'path': 'slide1region9', 'name': 'slide1region9 Embryo 1', 'rot_params': (None, False, True)},
        {'path': 'slide2region0', 'name': 'slide2region0 Embryo 1 left', 'rot_params': (-160, False, True)},
        
        # Original Script 2 tasks
        {'path': 'slide2region0', 'name': 'slide2region0 Embryo 1 right', 'filter_coords': (9500, 10180, 4500, 5300)},
        {'path': 'slide1region9', 'name': 'slide1region9 Embryo 2', 'filter_coords': (1400, 2300, 7800, 8550)},
        {'path': 'slide1region7', 'name': 'slide1region7 Embryo 2', 'filter_coords': (11000, 11700, 8160, 8900)},
        {'path': 'slide1region8', 'name': 'slide1region8 Embryo 2', 'filter_coords': (4300, 5000, 9500, 10250)},
        {'path': 'slide1region5', 'name': 'slide1region5 Embryo 2', 'filter_coords': (5400, 6200, 5400, 6180)},
        {'path': 'slide1region8', 'name': 'slide1region8 Embryo 4', 'filter_coords': (4800, 5700, 11000, 12000)},
        {'path': 'slide2region4', 'name': 'slide2region4 Embryo 4', 'filter_coords': (8800, 9500, 4300, 5000)},
        
        # Original Script 3 tasks
        {'path': 'slide2region4', 'name': 'slide2region4 Embryo 4', 'rot_params': (60, False, False), 'filter_coords': (2450, 3300, 5200, 5670)},
        {'path': 'slide1region5', 'name': 'slide1region5 Embryo 4', 'rot_params': (165, False, False), 'filter_coords': (5900, 6600, 6900, 7600)},
        {'path': 'slide2region5', 'name': 'slide2region5 Embryo 4', 'rot_params': (60, False, False), 'filter_coords': (4350, 5000, 7750, 8300)},
        {'path': 'slide2region4', 'name': 'slide2region4 Embryo 2', 'rot_params': (100, False, False), 'filter_coords': (1080, 1800, 5600, 6400)},
        {'path': 'slide2region5', 'name': 'slide2region5 Embryo 2', 'rot_params': (None, False, False), 'filter_coords': (3100, 3850, 8700, 9280)},
        {'path': 'slide1region4', 'name': 'slide1region4 Embryo 4', 'rot_params': (None, False, False), 'filter_coords': (500, 1200, 3000, 3800)},
    ]

    for task in processing_tasks:
        process_embryo_region(
            data_path_relative=task['path'],
            section_name=task['name'],
            rotation_params=task.get('rot_params'),
            filter_coords=task.get('filter_coords')
        )
