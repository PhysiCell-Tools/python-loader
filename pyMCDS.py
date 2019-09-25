import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import scipy.io as sio
import sys
from pathlib import Path

class pyMCDS:
    '''
    This class contains a dictionary of dictionaries that contains all of the 
    output from a single time step of a PhysiCell Model. This class assumes that
    all output files are stored in the same directory. Data is loaded by reading
    the .xml file for a particular timestep.

    Parameters
    ----------
    xml_name : string
        String containing the name of the xml file without the path
    output_path : string (default: '.')
        String containing the path (relative or absolute) to the directory
        where PhysiCell output files are stored
    
    Attributes
    ----------
    data : dictionary
        Hierarchical container for all of the data retrieved by parsing the xml
        file and the files referenced therein.
    '''
    def __init__(self, xml_file, output_path='.'):
        self.data = self._read_xml(xml_file, output_path)

    ## METADATA RELATED FUNCTIONS

    def get_time(self):
        return self.data['metadata']['current_time']

    ## MESH RELATED FUNCTIONS

    def get_mesh(self, flat=False):
        '''
        Return a meshgrid of the computational domain. Can return either full
        3D or a 2D plane for contour plots.

        Parameters
        ----------
        flat : bool
            If flat is set to true, we return only the x and y meshgrid.
            Otherwise we return x, y, and z

        Returns
        -------
        xx, yy, zz : three arrays, each [nx_voxel, ny_voxel, nz_voxel]
            Each array contains coordinates to describe the x, y, and z values 
            of every voxel center in the computational domain.
        
        if flat == True
        xx, yy : two arrays, each [nx_voxel, ny_voxel]
            Each array contains coordinates to describe the x and y values of
            flattened "voxels" on a plane
        '''
        if flat == True:
            xx = self.data['mesh']['x_coordinates'][:, :, 0]
            yy = self.data['mesh']['y_coordinates'][:, :, 0]

            return xx, yy

        # if we dont want a plane just return appropriate values
        else:
            xx = self.data['mesh']['x_coordinates']
            yy = self.data['mesh']['y_coordinates']
            zz = self.data['mesh']['z_coordinates']

            return xx, yy, zz

    def get_2D_mesh(self):
        '''
        This functional returns the x, y meshgrid as two numpy arrays. It is 
        identical to get_mesh with the option flat=True

        Returns
        -------
        xx : ndarray [nx_voxels, ny_voxels]
            x-coordinates for every voxel center in the computational domain
        
        yy : ndarray [nx_voxels, ny_voxels]
            y-coordinates for every voxel center in the computational domain
        '''
        xx = self.data['mesh']['x_coordinates'][:, :, 0]
        yy = self.data['mesh']['y_coordinates'][:, :, 0]

        return xx, yy

    def get_linear_coordinates(self):
        '''
        Return all of the unique coordinates of voxel centers for x-, y-, and,
        z-directions in the mesh.

        Returns
        -------
        X : ndarray [nx_voxel]
            ordered x values of voxel centers
        Y : ndarray [ny_voxel]
            ordered y values of voxel centers
        Z : ndarray [nz_voxels]
            ordered z values of voxel centers
        '''
        xx, yy, zz = self.get_mesh()

        X = np.unique(xx)
        Y = np.unique(yy)
        Z = np.unique(zz)

        return X, Y, Z

    def get_mesh_spacing(self):
        '''
        Returns the space in between voxel centers for the mesh in terms of the
        mesh's spatial units.

        Returns
        -------
        dx : float
            Distance between voxel centers in the same units as the other 
            spatial measurements
        '''
        X, Y, Z = self.get_linear_coordinates()
        dx = (X.max() - X.min()) / X.shape[0]
        dy = (Y.max() - Y.min()) / Y.shape[0]
        dz = (Z.max() - Z.min()) / Z.shape[0]
        
        if np.abs(dx - dy) > 1e-10 or np.abs(dy - dz) > 1e-10 \
            or np.abs(dx - dz) > 1e-10:
            print('Warning: grid spacing may be axis dependent.')
        
        return round(dx)

    ## MICROENVIRONMENT RELATED FUNCTIONS

    def get_substrate_names(self):
        '''
        Returns list of chemical species in microenvironment

        Returns
        -------
        species_list : array-like (str) [n_species,]
            Contains names of chemical species in microenvironment
        '''
        species_list = []
        for name in self.data['continuum_variables']:
            species_list.append(name)

        return species_list
    
    def get_concentrations(self, species_name, z_slice=None):
        '''
        Returns the concentration array for the specified chemical species
        in the microenvironment. Can return either the whole 3D picture, or
        a 2D plane of concentrations

        Parameters
        ----------
        species_name : string
            Name of the chemical species to get concentrations for
        
        z_slice : float
            z-axis position to use as plane for 2D output
        Returns
        -------
        conc_arr : ndarray (np.float) [x_voxels, y_voxels, z_voxels]
            Contains the concentration of the specified chemical in each voxel.
            The array spatially maps to a meshgrid of the voxel centers.
        '''
        if z_slice is not None:
            # check to see that z_slice is a valid plane
            zz = self.data['mesh']['z_coordinates']
            assert z_slice in zz, 'Specified z_slice {} not in z_coordinates'.format(z_slice)

            # do the processing if its ok
            mask = zz == z_slice
            full_conc = self.data['continuum_variables'][species_name]['data']
            conc_arr = full_conc[mask].reshape((zz.shape[0], zz.shape[1]))
        # its much simpler if we just want the whole thing
        else:
            conc_arr = self.data['continuum_variables'][species_name]['data']

        return conc_arr

    ## CELL RELATED FUNCTIONS

    def get_cell_df(self):
        '''
        Builds DataFrame from data['discrete_cells']

        Returns
        -------
        cells_df : pd.Dataframe [n_cells, n_variables]
            Dataframe containing the cell data for all cells at this time step
        '''
        cells_df = pd.DataFrame(self.data['discrete_cells'])
        return cells_df
    
    def get_cell_variables(self):
        '''
        Returns the names of all of the cell variables tracked in 'diccrete cells'
        dictionary

        Returns
        -------
        var_list : list (string) [n_variables]
            Contains the names of the cell variables
        '''
        var_list = []
        for name in self.data['discrete_cells']:
            var_list.append(name)
        return var_list
        
    def _read_xml(self, xml_file, output_path='.'):
        '''
        Does the actual work of initializing MultiCellDS by parsing the xml
        '''

        output_path = Path(output_path)
        xml_file = output_path / xml_file
        tree = ET.parse(xml_file)

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
        # THIS DOESNT ACTUALLY DO ANYTHING YET
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
            raise FileNotFoundError(
                "No such file or directory:\n'{}' referenced in '{}'".format(voxel_path, xml_file))
            sys.exit(1)

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
        # of species being tracked. the first 3 rows represent (x, y, z) of voxel
        # centers. The fourth row contains the voxel volume. The 5th row and up will
        # contain values for that species in that voxel.
        me_file = file_node.text
        me_path = output_path / me_file
        # Changes here
        try:
            me_data = sio.loadmat(me_path)['multiscale_microenvironment']
        except:
            raise FileNotFoundError(
                "No such file or directory:\n'{}' referenced in '{}'".format(me_path, xml_file))
            sys.exit(1)

        var_children = variables_node.findall('variable')

        # we're going to need the linear x, y, and z coordinates later
        # but we dont need to get them in the loop
        X, Y, Z = np.unique(xx), np.unique(yy), np.unique(zz)

        for si, species in enumerate(var_children):
            species_name = species.get('name')
            MCDS['continuum_variables'][species_name] = {}
            MCDS['continuum_variables'][species_name]['units'] = species.get(
                'units')

            # initialize array for concentration data
            MCDS['continuum_variables'][species_name]['data'] = np.zeros(xx.shape)

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
            # iterate over each voxel
            for vox_idx in range(MCDS['mesh']['voxels']['centers'].shape[1]):
                # find the center
                center = MCDS['mesh']['voxels']['centers'][:, vox_idx]

                i = np.where(np.abs(center[0] - X) < 1e-10)[0][0]
                j = np.where(np.abs(center[1] - Y) < 1e-10)[0][0]
                k = np.where(np.abs(center[2] - Z) < 1e-10)[0][0]

                MCDS['continuum_variables'][species_name]['data'][j, i, k] \
                    = me_data[4+si, vox_idx]

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
        data_labels = []
        # iterate over 'label's which are children of 'labels' these will be used to
        # label data arrays
        for label in cell_node.find('labels').findall('label'):
            # I don't like spaces in my dictionary keys
            fixed_label = label.text.replace(' ', '_')
            if int(label.get('size')) > 1:
                # tags to differentiate repeated labels (usually space related)
                dir_label = ['_x', '_y', '_z']
                for i in range(int(label.get('size'))):
                    data_labels.append(fixed_label + dir_label[i])
            else:
                data_labels.append(fixed_label)

        # load the file
        cell_file = cell_node.find('filename').text
        cell_path = output_path / cell_file
        try:
            cell_data = sio.loadmat(cell_path)['cells']
        except:
            raise FileNotFoundError(
                "No such file or directory:\n'{}' referenced in '{}'".format(cell_path, xml_file))
            sys.exit(1)

        for col in range(len(data_labels)):
            MCDS['discrete_cells'][data_labels[col]] = cell_data[col, :]

        return MCDS
