import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy import io
import sys
import warnings
from pathlib import Path

class pyMCDS:
    """
    This class contains a dictionary of dictionaries that contains all of the
    output from a single time step of a PhysiCell Model. This class assumes that
    all output files are stored in the same directory. Data is loaded by reading
    the .xml file for a particular timestep.

    Parameters
    ----------
    xml_name: str
        String containing the name of the xml file without the path
        or whith path (in this case output_path has to be the default).

    output_path: str, optional
        String containing the path (relative or absolute) to the directory
        where PhysiCell output files are stored (default= ".")

    microenv: boole, optional
        Should the microenvironment be extracted too?
        Set microenv to False will speed up processing.
        Set microenv to False will behave similar to the original pyMCDS_cells.py script.
        Default setting is True.


    Attributes
    ----------
    data : dict
        Hierarchical container for all of the data retrieved by parsing the xml
        file and the files referenced therein.
    """
    def __init__(self, xml_file, output_path='.', microenv=True):
        self.data = self._read_xml(xml_file, output_path, microenv)

    ## METADATA RELATED FUNCTIONS

    def get_time(self):
        """
        Return BUE simulated time in secounds or minutes or hour?

        Returns
        -------
        BUE
        """
        return self.data['metadata']['current_time']

    ## MESH RELATED FUNCTIONS

    def get_mesh(self, flat=False):
        """
        Return a meshgrid of the computational domain. Can return either full
        3D or a 2D plane for contour plots.

        Parameters
        ----------
        flat : bool
            If flat is set to True, only the x and y meshgrid will be returned.
            Otherwise the x, y, and z meshgrid will be returned.
            Default setting is False.

        Returns
        -------
        splitting : list length=2 if flat=True, else length=3
            Contains arrays of voxel center coordinates as meshgrid with shape
            [nx_voxel, ny_voxel, nz_voxel] or [nx_voxel, ny_voxel] if flat=True.
        """
        if flat:
            ar_m = self.data['mesh']['x_coordinates'][:, :, 0]
            ar_n = self.data['mesh']['y_coordinates'][:, :, 0]

            return [ar_m, ar_n]

        # if we dont want a plane just return appropriate values
        else:
            ar_m = self.data['mesh']['x_coordinates']
            ar_n = self.data['mesh']['y_coordinates']
            ar_p = self.data['mesh']['z_coordinates']

            return [ar_m, ar_n, ar_p]

    def get_mesh_2D(self):
        """
        This function returns the x, y meshgrid as two numpy arrays. It is
        identical to get_mesh with the option flat=True

        Returns
        -------
        splitting : list length=2
            Contains arrays of voxel center coordinates in x and y dimensions
            as meshgrid with shape [nx_voxel, ny_voxel]
        """
        return self.get_mesh(flat=True)

    def get_linear_voxels(self):
        """
        Helper function to quickly grab voxel centers array stored linearly as
        opposed to meshgrid-style.

        Note: 
        Ther oder of the flattened coordinates is neither C (row major) nor Fortran (column major) style.
        However, the coordinate set is complete.

        Returns
        -------
        flattend array of voxel position
        """
        # nmp
        ar_n, ar_p, ar_m, = self.data['mesh']['voxels']['centers'] 
        return [ar_m, ar_n, ar_p]

    def get_mesh_spacing(self):
        """
        Returns the space in between voxel centers for the mesh in terms of the
        mesh's spatial units. Assumes that voxel centers fall on integer values.

        Returns
        -------
        dm : float
            Distance between voxel centers in the same units as the other
            spatial measurements
        """
        ar_m, ar_n, ar_p = self.get_linear_voxels()

        dm = np.round((ar_m.max() - ar_m.min()) / ar_m.shape[0])
        dn = np.round((ar_n.max() - ar_n.min()) / ar_n.shape[0])
        dp = np.round((ar_p.max() - ar_p.min()) / ar_p.shape[0])

        #if np.abs(dm - dn) > 1e-10 or np.abs(dn - dp) > 1e-10 or np.abs(dm - dp) > 1e-10:
        #    print('Warning: grid spacing may be axis dependent.')

        return [dm, dn, dp]

    def get_containing_voxel_ijk(self, x, y, z):
        """
        Internal function to get the meshgrid indices for the center of a voxel
        that contains the given position.

        Note that pyMCDS stores meshgrids as 'cartesian'
        (indexing='xy' in np.meshgrid) which means that we will have
        to use these indices as [j, i, k] on the actual meshgrid objects

        Parameters
        ----------
        x : float
            x-coordinate for the position
        y : float
            y-coordinate for the position
        z : float
            z-coordinate for the position

        Returns
        -------
        ijk : list length=3
            contains the i, j, and k indices for the containing voxel's center
        """
        ar_m, ar_n, ar_p = self.get_mesh()
        dm, dn, dp = self.get_mesh_spacing()

        if x > ar_m.max():
            warnings.warn(f'Position out of bounds: x out of bounds in pyMCDS._get_voxel_idx({x}, {y}, {z}). Setting x = x_max!')
            x = ar_m.max()
        elif x < ar_m.min():
            warnings.warn(f'Position out of bounds: x out of bounds in pyMCDS._get_voxel_idx({x}, {y}, {z}). Setting x = x_min!')
            x = ar_m.min()
        elif y > ar_n.max():
            warnings.warn(f'Position out of bounds: y out of bounds in pyMCDS._get_voxel_idx({x}, {y}, {z}). Setting y = y_max!')
            y = ar_n.max()
        elif y < ar_n.min():
            warnings.warn(f'Position out of bounds: y out of bounds in pyMCDS._get_voxel_idx({x}, {y}, {z}). Setting y = y_min!')
            y = ar_n.min()
        elif z > ar_p.max():
            warnings.warn(f'Position out of bounds: z out of bounds in pyMCDS._get_voxel_idx({x}, {y}, {z}). Setting z = z_max!')
            z = ar_p.max()
        elif z < ar_p.min():
            warnings.warn(f'Position out of bounds: z out of bounds in pyMCDS._get_voxel_idx({x}, {y}, {z}). Setting z = z_min!')
            z = ar_p.min()

        i = int(np.round((x - ar_m.min()) / dm))
        j = int(np.round((y - ar_n.min()) / dn))
        k = int(np.round((z - ar_p.min()) / dp))

        return [i, j, k]


    ## MICROENVIRONMENT RELATED FUNCTIONS

    def get_substrate_names(self):
        """
        Returns list of chemical species in microenvironment

        Returns
        -------
        ls_species : list (str)
            Contains names of chemical species in microenvironment
        """
        ls_species = sorted(self.data['continuum_variables'].keys())
        return ls_species

    def get_concentrations(self, species_name, z_slice=None):
        """
        Returns the concentration array for the specified chemical species
        in the microenvironment. Can return either the whole 3D picture, or
        a 2D plane of concentrations.

        Parameters
        ----------
        species_name : str
            Name of the chemical species for which to get concentrations

        z_slice : float
            z-axis position to use as plane for 2D output.
            This value must match a plane of voxel centers in the z-axis.
            BUE: what about the None case.

        Returns
        -------
        ar_conc : array (np.float) shape=[nx_voxels, ny_voxels, nz_voxels]
            Contains the concentration of the specified chemical in each voxel.
            The array spatially maps to a meshgrid of the voxel centers.
        """
        ar_conc = self.data['continuum_variables'][species_name]['data']

        if not (z_slice is None):
            # check to see that z_slice is a valid plane
            ar_p = self.data['mesh']['z_coordinates']
            assert z_slice in ar_p, 'Specified z_slice {} not in z_coordinates'.format(z_slice)

            # do the processing if its ok
            mask = ar_p == z_slice
            ar_conc = ar_conc[mask].reshape((ar_p.shape[0], ar_p.shape[1]))

        return ar_conc

    def get_concentrations_at(self, x, y, z=0):
        """
        Return concentrations of each chemical species inside a particular voxel
        that contains the point described in the arguments.

        Parameters
        ----------
        x : float
            x-position for the point of interest
        y : float
            y_position for the point of interest
        z : float
            z_position for the point of interest

        Returns
        -------
        ar_concs : array, shape=[n_substrates,]
            array of concentrations in the order given by get_substrate_names()
        """
        i, j, k = self.get_containing_voxel_ijk(x, y, z)
        sub_name_list = self.get_substrate_names()
        ar_concs = np.zeros(len(sub_name_list))

        for ix, sub_name in enumerate(sub_name_list):
            ar_concs[ix] = self.get_concentrations(sub_name)[j, i, k]

        return ar_concs

    def get_concentrations_df(self, z_slice=None):
        """
        BUE

        Parameters
        ----------
        z_slice : float
            BUE check out get_concentrations

        Returns
        -------
        df_conc : DataFrame shape=[]
            BUE
        """
        # flatten mesh coordnates
        ar_m, ar_n, ar_p = self.get_mesh()
        ar_m = ar_m.flatten(order='C')
        ar_n = ar_n.flatten(order='C')
        ar_p = ar_p.flatten(order='C')
        # get min
        m_min = ar_m.min()
        n_min = ar_n.min()
        p_min = ar_p.min()
        # get voxel spacing
        dm, dn, dp = self.get_mesh_spacing()
        # get voxel coordinates
        ai_i = ((ar_m - m_min) / ds).astype(int)
        ai_j = ((ar_n - n_min) / ds).astype(int)
        ai_k = ((ar_p - p_min) / ds).astype(int)

        # handle coordinates
        ls_column = [
            'voxel_i','voxel_j','voxel_k',
            'voxel_position_m','voxel_position_n','voxel_position_p'
        ]
        la_data = [ai_i, ai_j, ai_k, ar_m, ar_n, ar_p]
        # handle concentraions
        for substrat_name in self.get_substrate_names():
            ls_column.append(substrat_name)
            la_data.append(ar_conc.flatten(order='C'))

        # generate data frame 
        aa_data  = np.array(la_data)
        df_conc = pd.DataFrame(aa_data.T, columns=ls_column)
        #d_dtype = {'voxel_i': int, 'voxel_j': int, 'voxel_k': float}  # bue: voxel_positions are all real.
        #df_conc = conc_df.astype(d_dtype)

        # filter
        if not (z_slice is None):
            assert z_slice in ar_p, 'Error: specifier z_slice {} not in z_coordinates'.format(z_slice)
            df_conc = df_conc.loc[df_conc.voxel_position_p == z_slice, :]

        # output
        return df_conc



    ## CELL RELATED FUNCTIONS

    def get_cell_variables(self):
        """
        Returns the names of all of the cell variables tracked in ['discrete cells']
        dictionary

        Returns
        -------
        ls_variables : list, shape=[n_variables]
            Contains the names of the cell variables
        """
        ls_variables = sorted(self.data['discrete_cells'].keys())
        return ls_variables

    def get_cell_df(self):
        """
        Builds DataFrame from data['discrete_cells']

        Returns
        -------
        df_cell : DataFrame, shape=[n_cells, n_variables]
            Dataframe containing the cell data for all cells at this time step
        """
        # get cell position and more
        df_cell = pd.DataFrame(self.data['discrete_cells'])
        df_voxel = df_cell.loc[:,['position_x','position_y','position_z']].copy()

        # get mesh spacing
        dm, dn, dp = self.get_mesh_spacing()

        # get mesh and voxel min max values
        ar_m, ar_n, ar_p = self.get_mesh()
        m_min = ar_m.min()
        n_min = ar_n.min()
        p_min = ar_p.min()
        i_min = 0
        j_min = 0
        k_min = 0
        i_max = ar_m.flatten().shape[0]
        j_max = ar_n.flatten().shape[0]
        k_max = ar_p.flatten().shape[0]

        # get voxel for each cell
        df_voxel.loc[:,'voxel_i'] = np.round((df_voxel.loc[:,'position_x'] - m_min) / dm)
        df_voxel.loc[:,'voxel_j'] = np.round((df_voxel.loc[:,'position_y'] - n_min) / dn)
        df_voxel.loc[:,'voxel_k'] = np.round((df_voxel.loc[:,'position_z'] - p_min) / dp)
        df_voxel.loc[(df_voxel.voxel_i > i_max), 'voxel_i'] =  i_max
        df_voxel.loc[(df_voxel.voxel_i < i_min), 'voxel_i'] =  i_min
        df_voxel.loc[(df_voxel.voxel_j > j_max), 'voxel_j'] =  j_max
        df_voxel.loc[(df_voxel.voxel_j < j_min), 'voxel_j'] =  j_min
        df_voxel.loc[(df_voxel.voxel_k > k_max), 'voxel_k'] =  k_max
        df_voxel.loc[(df_voxel.voxel_k < k_min), 'voxel_k'] =  k_min

        # get voxel position for each cell
        df_voxel['voxel_position_m'] = df_voxel.loc[:,'voxel_i'] * ds
        df_voxel['voxel_position_n'] = df_voxel.loc[:,'voxel_j'] * ds
        df_voxel['voxel_position_p'] = df_voxel.loc[:,'voxel_k'] * ds

        # output
        df_cell = pd.merge(cells_df, df_voxel, left_index=True, right_index=True)
        df_cell = df_cell.loc[:, sorted(df_cell.columns)]
        df_cell = df_cell.astype({
            'ID': int,
            'voxel_i': int, 'voxel_j': int, 'voxel_k': int
            # 'voxel_positions are all real.
        })
        df_cell.set_index('ID', inplace=True)
        return df_cell

    def get_cell_df_at(self, x, y, z=0):
        """
        Returns a dataframe for cells in the same voxel as the position given by
        x, y, and z.

        Parameters
        ----------
        x : float
            x-position for the point of interest
        y : float
            y_position for the point of interest
        z : float
            z_position for the point of interest

        Returns
        -------
        vox_df : DataFrame, shape=[n_cell_in_voxel, n_variables]
            cell dataframe containing only cells in the same voxel as the point
            specified by x, y, and z.
        """
        # get mesh and mesh spacing
        dm, dn, dp = self.get_mesh_spacing()
        ar_m, ar_n, ar_p = self.get_mesh()

        # get voxel coordinate
        i, j, k = self.get_containing_voxel_ijk(x, y, z)
        m = ar_m[j, i, k]
        n = ar_n[j, i, k]
        p = ar_p[j, i, k]

        # get voxel
        df_cell = self.get_cell_df()
        inside_voxel = (
            (df_cell['position_x'] < m + dm / 2) &
            (df_cell['position_x'] > m - dm / 2) &
            (df_cell['position_y'] < n + dn / 2) &
            (df_cell['position_y'] > n - dn / 2) &
            (df_cell['position_z'] < p + dp / 2) &
            (df_cell['position_z'] > p - dp / 2)
        )
        df_voxel = df_cell[inside_voxel]
        return df_voxel

    def _read_xml(self, xml_file, output_path='.', microenv=True):
        """
        Internal function. Does the actual work of initializing MultiCellDS by parsing the xml

        Parameters
        ----------
            BUE

        Returns
        -------
        MCDS: obj
            BUE
        """
        # file and path manipulation
        xml_file = xml_file.replace('\\','/')
        if (xml_file.find('/') > -1) and (output_path == '.'):
            ls_xml_file = xml_file.split('/')
            xml_file = ls_xml_file.pop(-1)
            output_path = '/'.join(ls_xml_file)
        output_path = Path(output_path)
        xml_pathfile = output_path / xml_file

        # read xml path/file
        tree = ET.parse(xml_pathfile)
        print('Reading {}'.format(xml_pathfile))

        root = tree.getroot()
        MCDS = {}

        # get current simulated time
        metadata_node = root.find('metadata')
        time_node = metadata_node.find('current_time')
        MCDS['metadata'] = {}
        MCDS['metadata']['current_time'] = float(time_node.text)
        MCDS['metadata']['time_units'] = time_node.get('units')

        # get current runtime
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

        # while we're at it, find the mesh
        coord_str = mesh_node.find('x_coordinates').text
        delimiter = mesh_node.find('x_coordinates').get('delimiter')
        x_coords = np.array(coord_str.split(delimiter), dtype=np.float64)

        coord_str = mesh_node.find('y_coordinates').text
        delimiter = mesh_node.find('y_coordinates').get('delimiter')
        y_coords = np.array(coord_str.split(delimiter), dtype=np.float64)

        coord_str = mesh_node.find('z_coordinates').text
        delimiter = mesh_node.find('z_coordinates').get('delimiter')
        z_coords = np.array(coord_str.split(delimiter), dtype=np.float64)

        # reshape into a mesh grid
        xx, yy, zz = np.meshgrid(x_coords, y_coords, z_coords)

        MCDS['mesh']['x_coordinates'] = xx
        MCDS['mesh']['y_coordinates'] = yy
        MCDS['mesh']['z_coordinates'] = zz

        # get microenvironment data
        if microenv:

            # Voxel data must be loaded from .mat file
            voxel_file = mesh_node.find('voxels').find('filename').text
            voxel_pathfile = output_path / voxel_file
            try:
                initial_mesh = io.loadmat(voxel_pathfile)['mesh']
            except:
                raise FileNotFoundError(
                    "No such file or directory:\n'{}' referenced in '{}'".format(voxel_pathfile, xml_pathfile))
                sys.exit(1)

            print('Reading {}'.format(voxel_pathfile))

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
            me_pathfile = output_path / me_file
            # Changes here
            try:
                me_data = io.loadmat(me_pathfile)['multiscale_microenvironment']
            except:
                raise FileNotFoundError(
                    "No such file or directory:\n'{}' referenced in '{}'".format(me_pathfile, xml_pathfile))
                sys.exit(1)

            print('Reading {}'.format(me_pathfile))

            var_children = variables_node.findall('variable')

            # we're going to need the linear x, y, and z coordinates later
            # but we dont need to get them in the loop
            X, Y, Z = np.unique(xx), np.unique(yy), np.unique(zz)

            for si, species in enumerate(var_children):
                species_name = species.get('name')
                MCDS['continuum_variables'][species_name] = {}
                MCDS['continuum_variables'][species_name]['units'] = species.get('units')

                print('Parsing {:s} data'.format(species_name))

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

        # get cell data
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

        print( 'working on discrete cell data...\n')

        MCDS['discrete_cells'] = {}
        data_labels = []
        # iterate over 'label's which are children of 'labels' these will be used to
        # label data arrays
        n = 0;
        for label in cell_node.find('labels').findall('label'):
            # I don't like spaces in my dictionary keys
            fixed_label = label.text.replace(' ', '_')
            nlabels = int(label.get('size'))
            if nlabels > 1:
                # tags to differentiate repeated labels (usually space related)
                # print("n=",n)
                spatial_type = False
                if( fixed_label == 'position' ):
                    spatial_type = True
                elif( fixed_label == 'orientation' ):
                    spatial_type = True
                elif( fixed_label == 'velocity' ):
                    spatial_type = True
                elif( fixed_label == 'migration_bias_direction' ):
                    spatial_type = True
                elif( fixed_label == 'motility_vector' ):
                    spatial_type = True

                if( nlabels == 3 and spatial_type == True ):
                    dir_label = ['_x', '_y', '_z']
                else:
                    dir_label = [];
                    for nn in range(100):
                        dir_label.append( '_%u' % nn )
                # print( dir_label )
                for i in range(int(label.get('size'))):
                    # print( fixed_label + dir_label[i] )
                    data_labels.append(fixed_label + dir_label[i])
            else:
                data_labels.append(fixed_label)
            # print(fixed_label)
            n += 1
        # load the file
        cell_file = cell_node.find('filename').text
        cell_pathfile = output_path / cell_file
        try:
            cell_data = io.loadmat(cell_pathfile)['cells']
        except:
            raise FileNotFoundError(
                "No such file or directory:\n'{}' referenced in '{}'".format(cell_pathfile, xml_pathfile))
            sys.exit(1)

        print('Reading {}'.format(cell_pathfile))

        for col in range(len(data_labels)):
            MCDS['discrete_cells'][data_labels[col]] = cell_data[col, :]

        return MCDS
