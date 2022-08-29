# %%
import xml.etree.ElementTree as TE
import numpy as np
import pandas as pd
import scipy.io as sio
from pathlib import Path

def read_MultiCellDS_xml(file_name, output_dir = '.'):
    
    output_path = Path(output_dir)
    xml_file = output_path / 'output00000001.xml'
    try:
        tree = TE.parse(xml_file)
    except:
        print('Data File:', xml_file, 'not found!')
        exit(0)

    root = tree.getroot()
    MCDS = {}

    # Get current simulated time
    metadata_node = root.find('metadata')
    time_node = metadata_node.find('current_time')
    MCDS['metadata'] = {}
    MCDS['metadata']['current_time'] = float(time_node.text)
    MCDS['metadata']['time_units'] = time_node.get('units')

    # Get current runtime
    time_node = metadata_node.find('current_runtime')
    MCDS['metadata']['current_runtime'] = float(time_node.text)
    MCDS['metadata']['runtime_units'] = time_node.get('units')

    # find the microenvironment node
    me_node = root.find('microenvironment')
    me_node = me_node.find('domain')

    # find the mesh node
    mesh_node = me_node.find('mesh')
    MCDS['metadata']['spatial_units'] = mesh_node.get('units')
    MCDS['mesh'] = {}

    # check for cartesian mesh
    cartesian = False
    mesh_type = mesh_node.get('type')
    if mesh_type == 'Cartesian':
        cartesian = True

    # while we're at it, find the mesh
    coord_str = mesh_node.find('x_coordinates').text
    delimiter = mesh_node.find('x_coordinates').get('delimiter')
    x_coords = np.array(coord_str.split(delimiter), dtype=np.float)

    coord_str = mesh_node.find('y_coordinates').text
    delimiter = mesh_node.find('y_coordinates').get('delimiter')
    y_coords = np.array(coord_str.split(delimiter), dtype=np.float)

    coord_str = mesh_node.find('z_coordinates').text
    delimiter = mesh_node.find('z_coordinates').get('delimiter')
    z_coords = np.array(coord_str.split(delimiter), dtype=np.float)

    # reshape into a mesh grid
    xx, yy, zz = np.meshgrid(x_coords, y_coords, z_coords)

    MCDS['mesh']['x_coordinates'] = xx
    MCDS['mesh']['y_coordinates'] = yy
    MCDS['mesh']['z_coordinates'] = zz

    # Voxel data must be loaded from .mat file
    voxel_file = mesh_node.find('voxels').find('filename').text
    voxel_path = output_path / voxel_file
    try:
        initial_mesh = sio.loadmat(voxel_path)['mesh']
    except:
        print('Data file', voxel_path, 'missing!')
        print('Referenced in', xml_file)
        exit(0)

    # center of voxel specified by first three rows [ x, y, z ]
    # volume specified by fourth row
    MCDS['mesh']['voxels'] = {}
    MCDS['mesh']['voxels']['centers'] = initial_mesh[:3, :]
    MCDS['mesh']['voxels']['volumes'] = initial_mesh[3, :]

    # Continuum_variables, unlike in the matlab version the individual chemical
    # species will be primarily accessed through their names e.g.
    # MCDS['continuum_variables']['oxygen']['units']
    # MCDS['continuum_variables']['glucose']['data']
    MCDS['continuum_variables'] = {}
    variables_node = me_node.find('variables')
    file_node = me_node.find('data').find('filename')

    # micro environment data is shape [4+n, len(voxels)] where n is the number
    # of species being tracked. the first 3 rows represent (x, y, z) of voxel centers. 
    # The fourth row contains the voxel volume. The 5th row and up will contain values
    # for that species in that voxel.
    me_file = file_node.text
    me_path = output_path / me_file
    try:
        me_data = sio.loadmat(me_path)['multiscale_microenvironment']
    except:
        print('Data file', me_path, 'missing!')
        print('Referenced in', xml_file)
        exit(0)

    var_children = variables_node.findall('variable')

    for i, species in enumerate(var_children):
        species_name = species.get('name')
        MCDS['continuum_variables'][species_name] = {}
        MCDS['continuum_variables'][species_name]['units'] = species.get('units')
        
        # travel down one level on tree
        species = species.find('physical_parameter_set')

        # diffusion data for each species
        MCDS['continuum_variables'][species_name]['diffusion_coefficient'] = {}
        MCDS['continuum_variables'][species_name]['diffusion_coefficient']['value'] \
            = float(species.find('diffusion_coefficient').text)
        MCDS['continuum_variables'][species_name]['diffusion_coefficient']['units'] \
            = species.find('diffusion_coefficient').get('units')

        # decay data for each species
        MCDS['continuum_variables'][species_name]['decay_rate'] = {}
        MCDS['continuum_variables'][species_name]['decay_rate']['value'] \
            = float(species.find('decay_rate').text)
        MCDS['continuum_variables'][species_name]['decay_rate']['units'] \
            = species.find('decay_rate').get('units')
        
        # store data from microenvironment file as numpy array
        MCDS['continuum_variables'][species_name]['data'] \
            = me_data[4+i, :].reshape(xx.shape)

    # in order to get to the good stuff we have to pass through a few different
    # hierarchal levels
    cell_node = root.find('cellular_information')
    cell_node = cell_node.find('cell_populations')
    cell_node = cell_node.find('cell_population')
    cell_node = cell_node.find('custom')
    # we want the PhysiCell data, there is more of it
    for child in cell_node.findall('simplified_data'):
        if child.get('source') == 'PhysiCell':
            cell_node = child
            break

    MCDS['discrete_cells'] = {}
    df_labels = []
    # iterate over 'label's which are children of 'labels' these will be used to
    # label dataframe columns
    for label in cell_node.find('labels').findall('label'):
        if int(label.get('size')) > 1:
            dir_label = ['_x', '_y', '_z']
            for i in range(int(label.get('size'))):
                df_labels.append(label.text + dir_label[i])
        else:
            df_labels.append(label.text)
            
    # load the file, we want it transposed from the matlab format
    cell_file = cell_node.find('filename').text
    cell_path = output_path / cell_file
    try:
        cell_data = sio.loadmat(cell_path)['cells'].T
    except:
        print('Data file', cell_path, 'missing!')
        print('Referenced in', xml_file)
        exit(0)
        
    # we will save the cell information as a dataframe in the MCDS object
    cell_df = pd.DataFrame(cell_data, columns=df_labels)
    MCDS['discrete_cells']['data_df'] = cell_df

    return MCDS
