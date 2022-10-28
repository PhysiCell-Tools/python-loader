###########
#
#
#
###########


# load library
import numpy as np
import pandas as pd
import pathlib
from scipy import io
import sys
import xml.etree.ElementTree as ET

# functions
def graphfile_parser(s_pathfile):
    """
    input:
        s_pathfile

    output:
        dei_graph

    descriptio:

    """
    # processing
    dei_graph = {}
    f = open(s_pathfile)
    for i, s_line in enumerate(f):
        #print('processing line:', s_line.strip())
        s_key, s_value = s_line.strip().split(':')
        ei_value = set()
        if len(s_value.strip()) :
            ei_value = set([int(s_id) for s_id in s_value.split(',')])
        dei_graph.update({int(s_key): ei_value})
    f.close()

    # output
    return(dei_graph)


# object classes
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
    def __init__(self, xmlfile, output_path='.', graph=True, microenv=True, verbose=True):
        self.microenv = microenv
        self.graph = graph
        self.verbose = verbose
        self.data = self._read_xml(xmlfile, output_path)


    ## METADATA RELATED FUNCTIONS

    def get_physicell_version(self):
        """
        """
        return self.data['metadata']['physicell_version']


    def get_timestamp(self):
        """
        Return BUE simulated time in secounds or minutes or hour?

        Returns
        -------
        BUE
        """
        return self.data['metadata']['created']


    def get_time(self):
        """
        Return BUE simulated time in secounds or minutes or hour?

        Returns
        -------
        BUE
        """
        return self.data['metadata']['current_time']


    def get_runtime(self):
        """
        Return BUE simulated time in secounds or minutes or hour?

        Returns
        -------
        BUE
        """
        return self.data['metadata']['current_runtime']


    ## GRAPH RELATED FUNCTIONS

    def get_attached_graph_dict(self):
        """
        BUE
        """
        return self.data['discrete_cells']['graph']['attached_cells']


    def get_neighbor_graph_dict(self):
        """
        BUE
        """
        return self.data['discrete_cells']['graph']['neighbor_cells']


    ## MESH RELATED FUNCTIONS

    def get_x_range(self):
        '''
        '''
        return self.data['mesh']['x_range']


    def get_y_range(self):
        '''
        '''
        return self.data['mesh']['y_range']


    def get_z_range(self):
        '''
        '''
        return self.data['mesh']['z_range']


    def get_mesh_m_range(self):
        '''
        '''
        return self.data['mesh']['m_range']


    def get_mesh_n_range(self):
        '''
        '''
        return self.data['mesh']['n_range']


    def get_mesh_p_range(self):
        '''
        '''
        return self.data['mesh']['p_range']


    def get_voxel_i_range(self):
        '''
        '''
        return self.data['mesh']['i_range']


    def get_voxel_j_range(self):
        '''
        '''
        return self.data['mesh']['j_range']


    def get_voxel_k_range(self):
        '''
        '''
        return self.data['mesh']['k_range']


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
            # if we dont want a plane just return appropriate values

        Returns
        -------
        splitting : list length=2 if flat=True, else length=3
            Contains arrays of voxel center coordinates as meshgrid with shape
            [nx_voxel, ny_voxel, nz_voxel] or [nx_voxel, ny_voxel] if flat=True.
        """
        if flat:
            ar_m = self.data['mesh']['m_coordinates'][:, :, 0]
            ar_n = self.data['mesh']['n_coordinates'][:, :, 0]

            return [ar_m, ar_n]

        else:
            ar_m = self.data['mesh']['m_coordinates']
            ar_n = self.data['mesh']['n_coordinates']
            ar_p = self.data['mesh']['p_coordinates']

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


    def get_mesh_linear(self):
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
        ar_m, ar_n, ar_p, = self.data['mesh']['mnp_coordinates']  # bue ['voxels']['centers']
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
        ar_m, ar_n, ar_p = self.get_mesh_linear()

        dm = np.round((ar_m.max() - ar_m.min()) / (len(set(ar_m)) - 1))
        dn = np.round((ar_n.max() - ar_n.min()) / (len(set(ar_n)) - 1))
        if (len(set(ar_p)) == 1):
            dp = 1
        else:
            dp = np.round((ar_p.max() - ar_p.min()) / (len(set(ar_p)) - 1))
        return [dm, dn, dp]


    def get_voxel_volume(self):
        """
        """
        ar_volume = np.unique(self.data['mesh']['volumes'])  # bue ['voxels']['volumes']
        if ar_volume.shape != (1,):
            sys.exit(f'Error @ pyMCDS.get_voxel_volume : mesh is not built out of a unique voxel volume {ar_volume}.')
        r_volume = ar_volume[0]
        return(r_volume)


    # bue def get_containing_voxel_ijk(self, x, y, z):
    def get_voxel_ijk(self, x, y, z):
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
        tr_m = self.get_mesh_m_range()
        tr_n = self.get_mesh_n_range()
        tr_p = self.get_mesh_p_range()
        dm, dn, dp = self.get_mesh_spacing()

        i = int(np.round((x - tr_m[0]) / dm))
        j = int(np.round((y - tr_n[0]) / dn))
        k = int(np.round((z - tr_p[0]) / dp))

        return [i, j, k]

    def is_in_mesh(self, x, y, z, halt=True):
        '''
        #print("_read_xml: bbox_coords= ",bbox_coords)  # rwh: [xmin ymin zmin xmax ymax zmax]
        '''
        b_isinmesh = True

        # check againt boundary box
        tr_xrange = self.get_x_range()
        tr_yrange = self.get_y_range()
        tr_zrange = self.get_z_range()

        if (x < tr_xrange[0]) or (x > tr_xrange[1]):
            print(f'Warning @ pyMCDS.is_in_mesh : x = {x} out of bounds: x-range is {tr_xrange}.')
            b_isinmesh = False
        elif (y < tr_yrange[0]) or (y > tr_yrange[1]):
            print(f'Warning @ pyMCDS.is_in_mesh : y = {y} out of bounds: y-range is {tr_yrange}.')
            b_isinmesh = False
        elif (z < tr_zrange[0]) or (z > tr_zrange[1]):
            print(f'Warning @ pyMCDS.is_in_mesh : z = {z} out of bounds: z-range is {tr_zrange}.')
            b_isinmesh = False

        # output
        if halt and not b_isinmesh:
            sys.exit('Processing stopped!')
        return(b_isinmesh)


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


    def get_substrate_df(self):
        """
        """
        # extract data
        dr_sub = {}
        for s_substrate in self.get_substrate_names():

            s_diffusion_key = f'{s_substrate}_diffusion_coefficient'
            s_diffusion_value = self.data['continuum_variables'][s_substrate]['diffusion_coefficient']['value']
            dr_sub.update({s_diffusion_key: [s_diffusion_value]})

            s_decay_key = f'{s_substrate}_decay_rate'
            s_decay_value = self.data['continuum_variables'][s_substrate]['decay_rate']['value']
            dr_sub.update({s_decay_key: [s_decay_value]})

        # output
        df_substrate = pd.DataFrame(dr_sub, index=['value']).T
        df_substrate.index.name = 'parameter'
        return(df_substrate)


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
        ar_p = self.data['mesh']['p_coordinates']

        if not (z_slice is None): 
            if (z_slice in ar_p):
                # filter by z_slice
                mask = ar_p == z_slice
                ar_conc = ar_conc[mask].reshape((ar_p.shape[0], ar_p.shape[1]))
            else:
                sys.exit(f'Error @ pyMCDS.get_concentrations : specified z_slice {z_slice} not in z_coordinates {np.unique(ar_p)}.')

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
        i, j, k = self.get_voxel_ijk(x, y, z)
        ls_substrate = self.get_substrate_names()
        ar_concs = np.zeros(len(ls_substrate))

        for n, s_substrate in enumerate(ls_substrate):
            ar_concs[n] = self.get_concentrations(s_substrate)[j, i, k]
            if self.verbose:
                print(f'pyMCD.get_concentrations_at(x={x},y={y},z={z}) > jkl: [{i},{j},{k}] > substrate: {s_substrate} {ar_concs[n]}')

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

        # check to see that z_slice is a valid plane
        if not (z_slice is None) and not (z_slice in ar_p):
            sys.exit(f'Error @ pyMCDS.get_concentrations_df : specified z_slice {z_slice} not in z_coordinates {np.unique(ar_p)}.')

        # get voxel spacing
        dm, dn, dp = self.get_mesh_spacing()

        # get voxel coordinates
        ai_i = ((ar_m - ar_m.min()) / dm)
        ai_j = ((ar_n - ar_n.min()) / dn)
        ai_k = ((ar_p - ar_p.min()) / dp)

        # handle coordinates
        ls_column = [
            'voxel_i','voxel_j','voxel_k',
            'mesh_center_m','mesh_center_n','mesh_center_p'
        ]
        la_data = [ai_i, ai_j, ai_k, ar_m, ar_n, ar_p]
        # handle concentraions
        for s_substrate in self.get_substrate_names():
            ls_column.append(s_substrate)
            ar_conc = self.get_concentrations(species_name=s_substrate, z_slice=None)
            la_data.append(ar_conc.flatten(order='C'))

        # generate data frame
        aa_data  = np.array(la_data)
        df_conc = pd.DataFrame(aa_data.T, columns=ls_column)
        d_dtype = {'voxel_i': int, 'voxel_j': int, 'voxel_k': int}  # bue: mesh_center are all real.
        df_conc = df_conc.astype(d_dtype)

        # filter
        if not (z_slice is None):
           df_conc = df_conc.loc[df_conc.mesh_center_p == z_slice, :]

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
        ls_variables = sorted(self.data['discrete_cells']['data'].keys())
        return ls_variables


    def get_cell_df(self):
        """
        Builds DataFrame from data['discrete_cells']['data']

        Returns
        -------
        df_cell : DataFrame, shape=[n_cells, n_variables]
            Dataframe containing the cell data for all cells at this time step
        """

        # get cell position and more
        df_cell = pd.DataFrame(self.data['discrete_cells']['data'])
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
        i_max = np.unique(ar_m).shape[0] - 1
        j_max = np.unique(ar_n).shape[0] - 1
        k_max = np.unique(ar_p).shape[0] - 1

        # get voxel for each cell
        df_voxel.loc[:,'voxel_i'] = np.round((df_voxel.loc[:,'position_x'] - m_min) / dm).astype(int)
        df_voxel.loc[:,'voxel_j'] = np.round((df_voxel.loc[:,'position_y'] - n_min) / dn).astype(int)
        df_voxel.loc[:,'voxel_k'] = np.round((df_voxel.loc[:,'position_z'] - p_min) / dp).astype(int)
        df_voxel.loc[(df_voxel.voxel_i > i_max), 'voxel_i'] =  i_max
        df_voxel.loc[(df_voxel.voxel_i < i_min), 'voxel_i'] =  i_min
        df_voxel.loc[(df_voxel.voxel_j > j_max), 'voxel_j'] =  j_max
        df_voxel.loc[(df_voxel.voxel_j < j_min), 'voxel_j'] =  j_min
        df_voxel.loc[(df_voxel.voxel_k > k_max), 'voxel_k'] =  k_max
        df_voxel.loc[(df_voxel.voxel_k < k_min), 'voxel_k'] =  k_min

        # merge voxel (inner join)
        df_cell = pd.merge(df_cell, df_voxel, on=['position_x', 'position_y', 'position_z'])

        # merge cell_density (left join)
        df_cellcount = df_cell.loc[:,['voxel_i','voxel_j','voxel_k','ID']].groupby(['voxel_i','voxel_j','voxel_k']).count().reset_index()
        ls_column = list(df_cellcount.columns)
        ls_column[-1] = 'cell_count_voxel'
        df_cellcount.columns = ls_column
        s_density = f"cell_density_{self.data['metadata']['spatial_units']}3"
        df_cellcount[s_density] = df_cellcount.loc[:,'cell_count_voxel'] / self.get_voxel_volume()
        df_cell = pd.merge(
            df_cell,
            df_cellcount,
            on = ['voxel_i', 'voxel_j', 'voxel_k'],
            how = 'left',
        )

        # microenviroment
        if self.microenv:
            # merge substrate (left join)
            df_sub = self.get_substrate_df()
            for s_index in df_sub.index:
                df_cell[s_index] = df_sub.loc[s_index,'value']

            # merge concentration (left join)
            df_conc = self.get_concentrations_df(z_slice=None)
            df_cell = pd.merge(
                df_cell,
                df_conc,
                on = ['voxel_i', 'voxel_j', 'voxel_k'],
                how = 'left',
            )

        # output
        df_cell = df_cell.loc[:, sorted(df_cell.columns)]
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
        i, j, k = self.get_voxel_ijk(x, y, z)
        m = ar_m[j, i, k]
        n = ar_n[j, i, k]
        p = ar_p[j, i, k]

        # get voxel
        df_cell = self.get_cell_df()
        inside_voxel = (
            (df_cell['position_x'] <= m + dm / 2) &
            (df_cell['position_x'] >= m - dm / 2) &
            (df_cell['position_y'] <= n + dn / 2) &
            (df_cell['position_y'] >= n - dn / 2) &
            (df_cell['position_z'] <= p + dp / 2) &
            (df_cell['position_z'] >= p - dp / 2)
        )
        df_voxel = df_cell[inside_voxel]
        return df_voxel


    ## UNIT OVERVIEW RELATED FUNCTION

    def get_unit_df(self):
        """
        """
        # extract data
        ds_unit = {}
        # units for metadata parameters
        ds_unit.update({'time': [self.data['metadata']['time_units']]})
        ds_unit.update({'runtime': [self.data['metadata']['runtime_units']]})
        ds_unit.update({'spatial_unit': [self.data['metadata']['spatial_units']]})

        # microenvironment
        if self.microenv:
            for s_substrate in self.get_substrate_names():
                # unit from substrate paramaters
                s_unit = self.data['continuum_variables'][s_substrate]['units']
                ds_unit.update({s_substrate: [s_unit]})

                # units from microenvironment paramaters
                s_diffusion_key = f'{s_substrate}_diffusion_coefficient'
                s_diffusion_unit = self.data['continuum_variables'][s_substrate]['diffusion_coefficient']['units']
                ds_unit.update({s_diffusion_key: [s_diffusion_unit]})

                s_decay_key = f'{s_substrate}_decay_rate'
                s_decay_unit = self.data['continuum_variables'][s_substrate]['decay_rate']['units']
                ds_unit.update({s_decay_key: [s_decay_unit]})

        # units from cell parameters
        ds_unit.update(self.data['discrete_cells']['units'])

        # output
        del ds_unit['ID']
        df_unit= pd.DataFrame(ds_unit, index=['unit']).T
        df_unit.index.name = 'parameter'
        df_unit.sort_index(inplace=True)
        return(df_unit)


    ## LOAD DATA

    def _read_xml(self, xmlfile, output_path='.'):
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
        xmlfile = xmlfile.replace('\\','/')
        if (xmlfile.find('/') > -1) and (output_path == '.'):
            ls_xmlfile = xmlfile.split('/')
            xmlfile = ls_xmlfile.pop(-1)
            output_path = '/'.join(ls_xmlfile)
        output_path = pathlib.Path(output_path)
        xmlpathfile = output_path / xmlfile

        # read xml path/file
        # 20221027 juliano: d = xmltodict.parse(open('PhysiCell_settings.xml').read(), process_namespaces=True)
        tree = ET.parse(xmlpathfile)
        if self.verbose:
            print(f'Reading: {xmlpathfile}')

        root = tree.getroot()
        MCDS = {}


        ####################
        # handle meta data #
        ####################

        if self.verbose:
            print('working on meta data ...')

        ### find the metadata node ###
        metadata_node = root.find('metadata')
        MCDS['metadata'] = {}

        # get physicell software version
        software_node = metadata_node.find('software')
        physicellv_node = software_node.find('version')
        MCDS['metadata']['physicell_version'] = physicellv_node.text

        # get timestamp
        time_node = metadata_node.find('created')
        MCDS['metadata']['created'] = time_node.text

        # get current simulated time
        time_node = metadata_node.find('current_time')
        MCDS['metadata']['current_time'] = float(time_node.text)
        MCDS['metadata']['time_units'] = time_node.get('units')

        # get current runtime
        time_node = metadata_node.find('current_runtime')
        MCDS['metadata']['current_runtime'] = float(time_node.text)
        MCDS['metadata']['runtime_units'] = time_node.get('units')

        # find the microenvironment node
        me_node = root.find('microenvironment')
        me_node = me_node.find('domain')


        ####################
        # handle mesh data #
        ####################

        if self.verbose:
            print('working on mesh data ...')

        ### find the mesh node ###
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
        ar_mmm, ar_nnn, ar_ppp = np.meshgrid(x_coords, y_coords, z_coords)
        
        MCDS['mesh']['m_coordinates'] = ar_mmm
        MCDS['mesh']['n_coordinates'] = ar_nnn
        MCDS['mesh']['p_coordinates'] = ar_ppp
        
        # get voxel range
        MCDS['mesh']['i_range'] = (0, len(set(ar_mmm.flatten())))
        MCDS['mesh']['j_range'] = (0, len(set(ar_nnn.flatten())))
        MCDS['mesh']['k_range'] = (0, len(set(ar_ppp.flatten())))

        # get mesh center range
        MCDS['mesh']['m_range'] = (ar_mmm.min(), ar_mmm.max())
        MCDS['mesh']['n_range'] = (ar_nnn.min(), ar_nnn.max())
        MCDS['mesh']['p_range'] = (ar_ppp.min(), ar_ppp.max())

        # get mesh bounding box range [xmin ymin zmin xmax ymax zmax]
        bboxcoor_str = mesh_node.find('bounding_box').text
        delimiter = mesh_node.find('bounding_box').get('delimiter')
        ar_bboxcoor = np.array(bboxcoor_str.split(delimiter), dtype=np.float64) 
        
        MCDS['mesh']['x_range'] = (ar_bboxcoor[0], ar_bboxcoor[3])
        MCDS['mesh']['y_range'] = (ar_bboxcoor[1], ar_bboxcoor[4])
        MCDS['mesh']['z_range'] = (ar_bboxcoor[2], ar_bboxcoor[5])

        # voxel data must be loaded from .mat file
        voxelfile = mesh_node.find('voxels').find('filename').text
        voxelpathfile = output_path / voxelfile
        try:
            initial_mesh = io.loadmat(voxelpathfile)['mesh']
            if self.verbose:
                print(f'Reading: {voxelpathfile}')
        except:
            raise FileNotFoundError(f'Error @ pyMCDS._read_xml : no such file or directory: {voxelpathfile}\nreferenced in: {xmlpathfile}.')
            sys.exit(1)


        # center of voxel specified by first three rows [ x, y, z ]
        # volume specified by fourth row
        #MCDS['mesh']['voxels'] = {}
        MCDS['mesh']['mnp_coordinates'] = initial_mesh[:3, :]   # bue ['voxels']['centers']
        MCDS['mesh']['volumes'] = initial_mesh[3, :]  # bue ['voxels']['volumes']


        ################################
        # handle microenvironment data #
        ################################

        if self.microenv:
            if self.verbose:
                print('working on microenvironment data ...')

            MCDS['continuum_variables'] = {}

            # Continuum_variables, unlike in the matlab version the individual chemical
            # species will be primarily accessed through their names e.g.
            # MCDS['continuum_variables']['oxygen']['units']
            # MCDS['continuum_variables']['glucose']['data']
            variables_node = me_node.find('variables')
            file_node = me_node.find('data').find('filename')

            # micro environment data is shape [4+n, len(voxels)] where n is the number
            # of species being tracked. the first 3 rows represent (x, y, z) of voxel
            # centers. The fourth row contains the voxel volume. The 5th row and up will
            # contain values for that species in that voxel.
            mefile = file_node.text
            mepathfile = output_path / mefile
            # Changes here
            try:
                me_data = io.loadmat(mepathfile)['multiscale_microenvironment']
                if self.verbose:
                    print(f'Reading: {mepathfile}')
            except:
                raise FileNotFoundError(f'Error @ pyMCDS._read_xml : no such file or directory: {mepathfile,}\nreferenced in: {xmlpathfile}.')
                sys.exit(1)


            var_children = variables_node.findall('variable')

            # we're going to need the linear x, y, and z coordinates later
            # but we dont need to get them in the loop
            ar_m = np.unique(ar_mmm)
            ar_n = np.unique(ar_nnn)
            ar_p = np.unique(ar_ppp)

            for i_s, species in enumerate(var_children):
                # i don't like spaces in species names!
                species_name = species.get('name').replace(' ', '_')

                MCDS['continuum_variables'][species_name] = {}
                MCDS['continuum_variables'][species_name]['units'] = species.get('units')

                if self.verbose:
                    print(f'Parsing: {species_name} data')

                # initialize array for concentration data
                MCDS['continuum_variables'][species_name]['data'] = np.zeros(ar_mmm.shape)

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
                for vox_idx in range(MCDS['mesh']['mnp_coordinates'].shape[1]):

                    # find the voxel coordinate
                    ar_center = MCDS['mesh']['mnp_coordinates'][:, vox_idx]
                    i = np.where(np.abs(ar_center[0] - ar_m) < 1e-10)[0][0]
                    j = np.where(np.abs(ar_center[1] - ar_n) < 1e-10)[0][0]
                    k = np.where(np.abs(ar_center[2] - ar_p) < 1e-10)[0][0]

                    # store value
                    MCDS['continuum_variables'][species_name]['data'][j, i, k] = me_data[4+i_s, vox_idx]


        ####################
        # handel cell data #
        ####################

        if self.verbose:
            print('working on discrete cell data ...')

        # in order to get to the good stuff we have to pass through a few different hierarchal levels
        cell_node = root.find('cellular_information')
        cell_node = cell_node.find('cell_populations')
        cell_node = cell_node.find('cell_population')
        cell_node = cell_node.find('custom')
        # we want the PhysiCell data, there is more of it
        for child in cell_node.findall('simplified_data'):
            if child.get('source') == 'PhysiCell':
                cellchild_node = child
                break

        MCDS['discrete_cells'] = {}

        # iterate over 'label's which are children of 'labels' these will be used to
        # label data arrays
        data_labels = []
        ds_unit = {}
        for label in cellchild_node.find('labels').findall('label'):
            # I don't like spaces in my dictionary keys!
            fixed_label = label.text.replace(' ', '_')
            nlabels = int(label.get('size'))
            s_unit = label.get('units')
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
                    s_label = fixed_label + dir_label[i]
                    # print(s_label)
                    data_labels.append(s_label)
                    ds_unit.update({s_label : s_unit})
            else:
                # print(fixed_label)
                data_labels.append(fixed_label)
                ds_unit.update({fixed_label : s_unit})

        # store unit
        MCDS['discrete_cells']['units'] = ds_unit

        # load the file
        cellfile = cellchild_node.find('filename').text
        cellpathfile = output_path / cellfile
        try:
            cell_data = io.loadmat(cellpathfile)['cells']
            if self.verbose:
                print(f'Reading: {cellpathfile}')
        except:
            raise FileNotFoundError(f'Error @ pyMCDS._read_xml : no such file or directory: {cellpathfile}\nreferenced in: {xmlpathfile}.')
            sys.exit(1)

        # store data
        MCDS['discrete_cells']['data'] = {}
        for col in range(len(data_labels)):
            MCDS['discrete_cells']['data'][data_labels[col]] = cell_data[col, :]


        ##################
        # get graph data #
        ##################

        if self.graph:

            if self.verbose:
                print('working on graph data ...')

            MCDS['discrete_cells']['graph'] = {}

            # neighborhod cell graph
            cellgraph_node = cell_node.find('neighbor_graph')
            cellfile = cellgraph_node.find('filename').text
            cellpathfile = output_path / cellfile
            try:
                dei_graph = graphfile_parser(s_pathfile=cellpathfile)
                if self.verbose:
                    print(f'Reading: {cellpathfile}')
            except:
                raise FileNotFoundError(f'Error @ pyMCDS._read_xml : no such file or directory: {cellpathfile}\nreferenced in: {xmlpathfile}.')
                sys.exit(1)

            # store data
            MCDS['discrete_cells']['graph'].update({'neighbor_cells': dei_graph})

            # attached cell graph
            cellgraph_node = cell_node.find('attached_cells_graph')
            cellfile = cellgraph_node.find('filename').text
            cellpathfile = output_path / cellfile
            try:
                dei_graph = graphfile_parser(s_pathfile=cellpathfile)
                if self.verbose:
                    print(f'Reading: {cellpathfile}')
            except:
                raise FileNotFoundError(f'Error @ pyMCDS._read_xml : no such file or directory: {cellpathfile}\nreferenced in: {xmlpathfile}.')
                sys.exit(1)

            # store data
            MCDS['discrete_cells']['graph'].update({'attached_cells': dei_graph})

        # output
        if self.verbose:
            print('done!')
        return MCDS
