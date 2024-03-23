#########
# title: pyMCDS.py
#
# language: python3
# date: 2022-08-22
# license: BSD-3-Clause
# authors: Patrick Wall, Randy Heiland, Furkan Kurtoglu, Paul Macklin, Elmar Bucher
#
# description:
#     pyMCDS.py definds an object class, able to load and access
#     within python a single time step from the PhysiCell model output folder.
#     pyMCDS.py was forked from the original PhysiCell-Tools python-loader
#     implementation and further developed.
#########


# load library
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import numpy as np
import os
import pandas as pd
from pcdl import pdplt
from scipy import io
import sys
import vtk
from vtkmodules.vtkCommonCore import vtkPoints
import xml.etree.ElementTree as ET
from pcdl.VERSION import __version__


# const physicell codec
# implemation based on PhysiCell/core/PhysiCell_constants.h.
ds_cycle_model = {
    '0' : 'advanced_Ki67_cycle_model',
    '1' : 'basic_Ki67_cycle_model',
    '2' : 'flow_cytometry_cycle_model',
    '3' : 'live_apoptotic_cycle_model',
    '4' : 'total_cells_cycle_model',
    '5' : 'live_cells_cycle_model',
    '6' : 'flow_cytometry_separated_cycle_model',
    '7' : 'cycling_quiescent_model',
}
ds_death_model = {
    '100' : 'apoptosis_death_model',
    '101' : 'necrosis_death_model',
    '102' : 'autophagy_death_model',
    '9999' : 'custom_cycle_model',
}

ds_cycle_phase = {
    '0' : 'Ki67_positive_premitotic',
    '1' : 'Ki67_positive_postmitotic',
    '2' : 'Ki67_positive',
    '3' : 'Ki67_negative',
    '4' : 'G0G1_phase',
    '5' : 'G0_phase',
    '6' : 'G1_phase',
    '7' : 'G1a_phase',
    '8' : 'G1b_phase',
    '9' : 'G1c_phase',
    '10' : 'S_phase',
    '11' : 'G2M_phase',
    '12' : 'G2_phase',
    '13' : 'M_phase',
    '14' : 'live',
    '15' : 'G1pm_phase',
    '16' : 'G1ps_phase',
    '17' : 'cycling',
    '18' : 'quiescent',
    '9999' : 'custom_phase',
}
ds_death_phase = {
    '100' : 'apoptotic',
    '101' : 'necrotic_swelling',
    '102' : 'necrotic_lysed',
    '103' : 'necrotic',
    '104' : 'debris',
}

# const physicell variable names
es_var_subs = {
    'chemotactic_sensitivities',
    'fraction_released_at_death',
    'fraction_transferred_when_ingested',
    'internalized_total_substrates',
    'net_export_rates',
    'saturation_densities',
    'secretion_rates',
    'uptake_rates',
}
es_var_death = {
    'death_rates',
}
es_var_cell = {
    'attack_rates',
    'cell_adhesion_affinities',
    'fusion_rates',
    'live_phagocytosis_rates',
    'transformation_rates',
}
es_var_spatial = {
    'migration_bias_direction',
    'motility_vector',
    'orientation',
    'position',
    'velocity',
}

# const physicell variable types
do_var_type = {
    # integer
    'ID': int,
    'cell_count_voxel': int,
    'chemotaxis_index': int,
    'maximum_number_of_attachments': int,
    'number_of_nuclei': int,
    # boolean
    'contact_with_basement_membrane': bool,
    'dead': bool,
    'is_motile': bool,
    # categorical
    'cell_type': str,  # id mapping
    'current_death_model': str,  # codec mapping
    'current_phase': str,  # codec mapping
    'cycle_model': str,  # codec mapping
}

# const coordinate variable names
es_coor_conc = {
    'ID',
    'voxel_i','voxel_j','voxel_k',
    'mesh_center_m','mesh_center_n','mesh_center_p',
    'time', 'runtime',
    'xmlfile',
}
es_coor_cell = {
    'ID',
    'voxel_i', 'voxel_j', 'voxel_k',
    'mesh_center_m', 'mesh_center_n', 'mesh_center_p',
    'position_x', 'position_y', 'position_z',
    'time', 'runtime',
    'xmlfile',
}


# functions
def graphfile_parser(s_pathfile):
    """
    input:
        s_pathfile: string
            path to and file name from graph.txt file.

    output:
        dei_graph: dictionary of sets of integers.
            object maps each cell ID to connected cell IDs.

    description:
        code parses PhysiCell's own graphs format and
        returns the content in a dictionary object.
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
    return dei_graph


# object classes
class pyMCDS:
    def __init__(self, xmlfile, output_path='.', custom_type={}, microenv=True, graph=True, physiboss=True, settingxml='PhysiCell_settings.xml', verbose=True):
        """
        input:
            xmlfile: string
                name of the xml file with or without path.
                in the with path case, output_path has to be set to the default!

            output_path: string; default '.'
                relative or absolute path to the directory where
                the PhysiCell output files are stored.

            custom_type: dictionary; default is {}
                variable to specify custom_data variable types
                other than float (int, bool, str) like this: {var: dtype, ...}.
                downstream float and int will be handled as numeric,
                bool as Boolean, and str as categorical data.

            microenv: boole; default True
                should the microenvironment be extracted?
                setting microenv to False will use less memory and speed up
                processing, similar to the original pyMCDS_cells.py script.

            graph: boole; default True
                should the graphs be extracted?
                setting graph to False will use less memory and speed up processing.

            physiboss: boole; default True
                should physiboss state data extracted, if found?
                setting physiboss to False will use less memory and speed up processing.

            settingxml: string; default PhysiCell_settings.xml
                from which settings.xml should the cell type ID label mapping,
                parameters, and units be extracted?
                set to None or False if the xml file is missing!

            verbose: boole; default True
                setting verbose to False for less text output, while processing.

        output:
            mcds: pyMCDS class instance
                all fetched content is stored at mcds.data.

        description:
            pyMCDS.__init__ will generate a class instance with a
            dictionary of dictionaries data structure that contains all
            output from a single PhysiCell model time step. furthermore,
            this class, and as such it's instances, offers functions
            to access the stored data.
            the code assumes that all related output files are stored in
            the same directory. data is loaded by reading the xml file
            for a particular time step and the therein referenced files.
        """
        self.path = None
        self.xmlfile = None
        self.custom_type = custom_type
        self.microenv = microenv
        self.graph = graph
        self.physiboss = physiboss
        if type(settingxml) is str:
            settingxml = settingxml.replace('\\','/').split('/')[-1]
        self.settingxml = settingxml
        self.verbose = verbose
        self.data = self._read_xml(xmlfile, output_path)
        self.get_concentration_df = self.get_conc_df

    def set_verbose_false(self):
        """
        input:

        output:
            set verbose false.

        description:
            function to set verbosity.
        """
        self.verbose = False

    def set_verbose_true(self):
        """
        input:

        output:
            set verbose true.

        description:
            function to set verbosity.
        """
        self.verbose = True


    ## METADATA RELATED FUNCTIONS ##

    def get_multicellds_version(self):
        """
        input:

        output:
            s_version : sting
                MultiCellDS xml version which stored the data.

        description:
            function returns as a string the MultiCellDS xml version
            that was used to store this data.
        """
        return self.data['metadata']['multicellds_version']


    def get_physicell_version(self):
        """
        input:

        output:
            s_version : sting
                PhysiCell version which generated the data.

        description:
            function returns as a string the PhysiCell version
            that was used to generate this data.
        """
        return self.data['metadata']['physicell_version']


    def get_timestamp(self):
        """
        input:

        output:
            s_timestap : sting
                timestamp from when this data was generated.

        description:
            function returns as a string the timestamp from when
            this data was generated.
        """
        return self.data['metadata']['created']


    def get_time(self):
        """
        input:

        output:
            r_time : floating point number
                simulation time in [min].

        description:
            function returns as a real number
            the simulation time in minutes.
        """
        return self.data['metadata']['current_time']


    def get_runtime(self):
        """
        input:

        output:
            r_time : floating point number
                wall time in [sec].

        description:
            function returns as a real number, the wall time in seconds
            the simulation took to run up to this time step.
        """
        return self.data['metadata']['current_runtime']


    ## MESH RELATED FUNCTIONS  ##

    def get_voxel_ijk_range(self):
        """
        input:

        output:
            lti_i : list of tuple of 2 integer numbers
                i-axis, j-aixs, and k-axis voxel range.

        decritpion:
            function returns in a list of tuples the lowest and highest
            i-axis, j-axis, and k-axis voxel value.
        """
        return self.data['mesh']['ijk_range'].copy()


    def get_mesh_mnp_range(self):
        """
        input:

        output:
            ltr_mnp : list of tuple of 2 floating point numbers
                m-axis, n-axis, and p-axis  mesh center range.

        decritpion:
            function returns in a list of tuples the lowest and highest
            m-axis, n-axis, and p-axis mesh center value.
        """
        return self.data['mesh']['mnp_range'].copy()


    def get_xyz_range(self):
        """
        input:

        output:
            ltr_xyz : list of tuple of 2 floating point numbers
                x-axis, y-axis, and z-axis position range.

        decritpion:
            function returns in a list of tuples the lowest and highest
            x-axis, y-axis, and z-axis position value.
        """
        return self.data['mesh']['xyz_range'].copy()


    def get_voxel_ijk_axis(self):
        """
        input:

        output:
            lai_ijk : list of 3 numpy arrays of integer numbers
                i-axis, j-axis, and k-axis voxel coordinates axis.

        description:
            function returns a list of voxel coordinate vectors,
            one for the i-axis, j-axis, and k-axis.
        """
        return self.data['mesh']['ijk_axis'].copy()


    def get_mesh_mnp_axis(self):
        """
        input:

        output:
            lar_mnp : list of 3 numpy arrays of floating point numbers
                m-axis, n-axis, and p-axis mesh center axis coordinates.

        description:
            function returns a list of mesh center vectors,
            one for the m-axis, n-axis, and p-axis.
        """
        return self.data['mesh']['mnp_axis'].copy()


    def get_mesh(self, flat=False):
        """
        input:
            flat : bool; default False
                if flat is True, only the m-axis mesh center
                and n-axis mesh center meshgrids will be returned.
                else the m, n, and p mesh center meshgrids will be returned.

        output:
            aar_meshgrid : 4-way (3D) or 3-way (2D) numpy arrays tensor of floating point numbers
                meshgrid shaped object, with the mesh center
                coordinate values from the m, n, p-axis or m, n-axis.

        description:
            function returns a numpy array of meshgrids each stores the
            mesh center coordinate values from one particular axis.
            the function can either return meshgrids for the full
            m, n, p 3D cube, or only the 2D planes along the p-axis.
        """
        if flat:
            ar_m = self.data['mesh']['mnp_grid'][0][:, :, 0]
            ar_n = self.data['mesh']['mnp_grid'][1][:, :, 0]
            return np.array([ar_m, ar_n]).copy()

        else:
            return self.data['mesh']['mnp_grid'].copy()


    def get_mesh_2D(self):
        """
        input:

        output:
            aar_meshgrid : 3-way numpy arrays tensor of floating point numbers
                meshgrid shaped objects, with the mesh center
                coordinate values from the m and n-axis.

        description:
            function is identical to the self.get_mesh(self, flat=True)
            function call.
        """
        return self.get_mesh(flat=True)


    def get_mesh_coordinate(self):
        """
        input:

        output:
            aar_meshaxis : numpy array of 3 one dimensional numpy floating point number arrays
                n, m, and p-axis mesh center coordinate vectors.

        description:
            function returns three vectors with mesh center coordinate values,
            one for each axis.
        """
        return self.data['mesh']['mnp_coordinate'].copy()


    def get_mesh_spacing(self):
        """
        input:

        output:
            lr_mnp_spacing: list of 3 floating point numbers
                mesh spacing in m, n, and p direction.

        description:
            function returns the distance in between mesh centers,
            in the spacial unit defined in the PhysiCell_settings.xml file.
        """
        tr_m_range, tr_n_range, tr_p_range = self.get_mesh_mnp_range()
        ar_m_axis, ar_n_axis, ar_p_axis = self.get_mesh_mnp_axis()

        dm = (tr_m_range[1] - tr_m_range[0]) / (ar_m_axis.shape[0] - 1)
        dn = (tr_n_range[1] - tr_n_range[0]) / (ar_n_axis.shape[0] - 1)
        if (len(set(tr_p_range)) == 1):
            dp = np.float64(1.0)
        else:
            dp = (tr_p_range[1] - tr_p_range[0]) / (ar_p_axis.shape[0] - 1)
        return [dm, dn, dp]


    def is_in_mesh(self, x, y, z, halt=False):
        """
        input:
            x: floating point number
                position x-coordinate.

            y: floating point number
                position y-coordinate.

            z: floating point number
                position z-coordinate.

            halt: boolean; default is False
                should program execution break or just spit out a warning,
                if position is not in mesh?

        output:
            b_isinmesh: boolean
                declares if the given coordinate is inside the mesh.

        description:
            function evaluates, if the given position coordinate
            is inside the boundaries. if the coordinate is outside the
            mesh, a warning will be printed. if additionally
            halt is set to True, program execution will break.
        """
        b_isinmesh = True

        # check against boundary box
        tr_x, tr_y, tr_z = self.get_xyz_range()

        if (x < tr_x[0]) or (x > tr_x[1]):
            print(f'Warning @ pyMCDS.is_in_mesh : x = {x} out of bounds: x-range is {tr_x}.')
            b_isinmesh = False
        elif (y < tr_y[0]) or (y > tr_y[1]):
            print(f'Warning @ pyMCDS.is_in_mesh : y = {y} out of bounds: y-range is {tr_y}.')
            b_isinmesh = False
        elif (z < tr_z[0]) or (z > tr_z[1]):
            print(f'Warning @ pyMCDS.is_in_mesh : z = {z} out of bounds: z-range is {tr_z}.')
            b_isinmesh = False

        # output
        if halt and not b_isinmesh:
            sys.exit('Processing stopped!')
        return b_isinmesh


    def get_voxel_spacing(self):
        """
        input:

        output:
            lr_ijk_spacing: list of 3 floating point numbers
                voxel spacing in i, j, and k direction.

        description:
            function returns the voxel width, height, depth measurement,
            in the spacial unit defined in the PhysiCell_settings.xml file.
        """
        r_volume = self.get_voxel_volume()
        dm, dn, _ = self.get_mesh_spacing()
        dp  = r_volume / (dm * dn)
        return [dm, dn, dp]


    def get_voxel_volume(self):
        """
        input:

        output:
            r_volume: floating point number
                voxel volume value related to the spacial unit
                defined in the PhysiCell_settings.xml file.

        description:
            function returns the volume value for a single voxel, related
            to the spacial unit defined in the PhysiCell_settings.xml file.
        """
        ar_volume = np.unique(self.data['mesh']['volumes'])
        if ar_volume.shape != (1,):
            sys.exit(f'Error @ pyMCDS.get_voxel_volume : mesh is not built out of a unique voxel volume {ar_volume}.')
        r_volume = ar_volume[0]
        return r_volume


    def get_voxel_ijk(self, x, y, z, is_in_mesh=True):
        """
        input:
            x: floating point number
                position x-coordinate.

            y: floating point number
                position y-coordinate.

            z: floating point number
                position z-coordinate.

            is_in_mesh: boolean; default is True
                should function check, if the given coordinate is in the mesh,
                and only calculate ijk values if is so?

        output:
            li_ijk : list of 3 integers
                i, j, k indices for the voxel
                containing the x, y, z position.

        description:
            function returns the meshgrid indices i, j, k
            for the given position x, y, z.
        """
        lr_ijk = None
        b_calc = True

        if is_in_mesh:
            b_calc = self.is_in_mesh(x=x, y=y, z=z, halt=False)

        if b_calc:
            tr_m, tr_n, tr_p = self.get_mesh_mnp_range()
            dm, dn, dp = self.get_voxel_spacing()

            i = int(np.round((x - tr_m[0]) / dm))
            j = int(np.round((y - tr_n[0]) / dn))
            k = int(np.round((z - tr_p[0]) / dp))

            lr_ijk = [i, j, k]

        return lr_ijk


    ## MICROENVIRONMENT RELATED FUNCTIONS ##

    def get_substrate_names(self):
        """
        input:

        output:
            ls_substrate: list of stings
                alphabetically ordered list of all tracked substrates.

        description:
            function returns all chemical species names,
            modeled in the microenvironment.
        """
        ls_substrate = sorted(self.data['continuum_variables'].keys())
        return ls_substrate


    def get_substrate_dict(self):
        """
        input:

        output:
            ds_substrate: dictionary of stings
                dictionary that maps substrate IDs to labels.

        description:
            function returns a dictionary that maps ID and name from all
            microenvironment_setup variables,
            specified in the PhysiCell_settings.xml file.
        """
        return self.data['metadata']['substrate']


    def get_substrate_df(self):
        """
        input:

        output:
            df_substrae: pandas dataframe
                one substrate per row and decay_rate and difusion_coefficient
                factors as columns.

        description:
            function returns a dataframe with each substrate's
            decay_rate and difusion_coefficient.
        """
        # extract data
        ls_column = ['substrate','decay_rate','diffusion_coefficient']
        ll_sub = []
        for s_substrate in self.get_substrate_names():
            s_decay_value = self.data['continuum_variables'][s_substrate]['decay_rate']['value']
            s_diffusion_value = self.data['continuum_variables'][s_substrate]['diffusion_coefficient']['value']
            ll_sub.append([s_substrate, s_decay_value, s_diffusion_value])

        # generate dataframe
        df_substrate = pd.DataFrame(ll_sub, columns=ls_column)
        df_substrate.set_index('substrate', inplace=True)
        df_substrate.columns.name = 'feature'

        # output
        return df_substrate


    def get_concentration(self, substrate, z_slice=None, halt=False):
        """
        input:
            substrate: string
                substrate name.

            z_slice: floating point number; default is None
                z-axis position to slice a 2D xy-plain out of the
                3D substrate concentration mesh. if None the
                whole 3D mesh will be returned.

            halt: boolean; default is False
                should program execution break or just spit out a warning,
                if z_slice position is not an exact mesh center coordinate?
                if False, z_slice will be adjusted to the nearest
                mesh center value, the smaller one, if the coordinate
                lies on a saddle point.

        output:
            ar_conc: numpy array of floating point numbers
                substrate concentration meshgrid or xy-plain slice
                through the meshgrid.

        description:
            function returns the concentration meshgrid, or a xy-plain slice
            out of the whole meshgrid, for the specified chemical species.
        """
        ar_conc = self.data['continuum_variables'][substrate]['data'].copy()

        # check if z_slice is a mesh center or None
        if not (z_slice is None):
            _, _, ar_p_axis = self.get_mesh_mnp_axis()
            if not (z_slice in ar_p_axis):
                print(f'Warning @ pyMCDS.get_concentration : specified z_slice {z_slice} is not an element of the z-axis mesh centers set {ar_p_axis}.')
                if halt:
                    sys.exit('Processing stopped!')
                else:
                    z_slice = ar_p_axis[abs(ar_p_axis - z_slice).argmin()]
                    print(f'z_slice set to {z_slice}.')

            # filter by z_slice
            _, _, ar_p_grid = self.get_mesh()
            mask = ar_p_grid == z_slice
            ar_conc = ar_conc[mask].reshape((ar_p_grid.shape[0], ar_p_grid.shape[1]))

        # output
        return ar_conc


    def get_concentration_at(self, x, y, z=0):
        """
        input:
            x: floating point number
                position x-coordinate.

            y: floating point number
                position y-coordinate.

            z: floating point number; default is 0
                position z-coordinate.

        output:
            ar_concs: numpy array of floating point numbers
                array of substrate concentrations in the order
                given by get_substrate_names().

        description:
            function return concentrations of each chemical species
            inside a particular voxel that contains the point specified
            in the arguments.
        """
        ar_concs = None

        # is coordinate inside the domain?
        b_calc = self.is_in_mesh(x=x, y=y, z=z, halt=False)
        if b_calc:

            # get voxel coordinate and substrate names
            i, j, k = self.get_voxel_ijk(x, y, z, is_in_mesh=False)
            ls_substrate = self.get_substrate_names()
            ar_concs = np.zeros(len(ls_substrate))

            # get substrate concentrations
            for n, s_substrate in enumerate(ls_substrate):
                ar_concs[n] = self.get_concentration(s_substrate)[j, i, k]
                if self.verbose:
                    print(f'pyMCD.get_concentration_at(x={x},y={y},z={z}) | jkl: [{i},{j},{k}] | substrate: {s_substrate} {ar_concs[n]}')

        # output
        return ar_concs


    def get_conc_df(self, z_slice=None, halt=False, values=1, drop=set(), keep=set()):
        """
        input:
            z_slice: floating point number; default is None
                z-axis position to slice a 2D xy-plain out of the
                3D substrate concentration mesh. if None the
                whole 3D mesh will be returned.

            halt: boolean; default is False
                should program execution break or just spit out a warning,
                if z_slice position is not an exact mesh center coordinate?
                if False, z_slice will be adjusted to the nearest
                mesh center value, the smaller one, if the coordinate
                lies on a saddle point.

            values: integer; default is 1
                minimal number of values a variable has to have to be outputted.
                variables that have only 1 state carry no information.
                None is a state too.

            drop: set of strings; default is an empty set
                set of column labels to be dropped for the dataframe.
                don't worry: essential columns like ID, coordinates
                and time will never be dropped.
                Attention: when the keep parameter is given, then
                the drop parameter has to be an empty set!

            keep: set of strings; default is an empty set
                set of column labels to be kept in the dataframe.
                set values=1 to be sure that all variables are kept.
                don't worry: essential columns like ID, coordinates
                and time will always be kept.

        output:
            df_conc : pandas dataframe
                dataframe stores all substrate concentrations in each voxel.

        description:
            function returns a dataframe with concentration values
            for all chemical species in all voxels. additionally, this
            dataframe lists voxel and mesh center coordinates.
        """
        # check keep and drop
        if (len(keep) > 0) and (len(drop) > 0):
            sys.exit(f"Error @ pyMCDS.get_conc_df : when keep is given {keep}, then drop has to be an empty set {drop}!")

        # check if z_slice is a mesh center or None
        if not (z_slice is None):
            _, _, ar_p_axis = self.get_mesh_mnp_axis()
            if not (z_slice in ar_p_axis):
                print(f'Warning @ pyMCDS.get_conc_df : specified z_slice {z_slice} is not an element of the z-axis mesh centers set {ar_p_axis}.')
                if halt:
                    sys.exit('Processing stopped!')
                else:
                    z_slice = ar_p_axis[abs(ar_p_axis - z_slice).argmin()]
                    print(f'z_slice set to {z_slice}.')

        # flatten mesh coordnates
        ar_m, ar_n, ar_p = self.get_mesh()
        ar_m = ar_m.flatten(order='C')
        ar_n = ar_n.flatten(order='C')
        ar_p = ar_p.flatten(order='C')

        # get mesh spacing
        dm, dn, dp = self.get_voxel_spacing()

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

        # handle concentrations
        for s_substrate in self.get_substrate_names():
            ls_column.append(s_substrate)
            ar_conc = self.get_concentration(substrate=s_substrate, z_slice=None)
            la_data.append(ar_conc.flatten(order='C'))

        # generate dataframe
        aa_data  = np.array(la_data)
        df_conc = pd.DataFrame(aa_data.T, columns=ls_column)
        df_conc['time'] = self.get_time()
        df_conc['runtime'] = self.get_runtime() / 60  # in min
        df_conc['xmlfile'] = self.xmlfile
        d_dtype = {'voxel_i': int, 'voxel_j': int, 'voxel_k': int}
        df_conc = df_conc.astype(d_dtype)

        # filter z_slice
        if not (z_slice is None):
           df_conc = df_conc.loc[df_conc.mesh_center_p == z_slice, :]

        # filter
        es_feature = set(df_conc.columns).difference(es_coor_conc)
        if (len(keep) > 0):
            es_delete = es_feature.difference(keep)
        else:
            es_delete = es_feature.intersection(drop)

        if (values > 1):
            for s_column in set(df_conc.columns).difference(es_coor_conc):
                if len(set(df_conc.loc[:,s_column])) < values:
                    es_delete.add(s_column)
        if self.verbose and (len(es_delete) > 0):
            print('es_delete:', es_delete)
        df_conc.drop(es_delete, axis=1, inplace=True)

        # output
        df_conc.sort_values(['voxel_i', 'voxel_j', 'voxel_k', 'time'], inplace=True)
        return df_conc


    def plot_contour(self, substrate, z_slice=0, vmin=None, vmax=None, alpha=1, fill=True, cmap='viridis', title=None, grid=True, xlim=None, ylim=None, xyequal=True, figsize=None, ax=None):
        """
        input:
            substrate: string
                substrate name.

            z_slice: floating point number; default is 0
                z-axis position to slice a 2D xy-plain out of the
                3D substrate concentration mesh. if z_slice position
                is not an exact mesh center coordinate, then z_slice
                will be adjusted to the nearest mesh center value,
                the smaller one, if the coordinate lies on a saddle point.

            vmin: floating point number; default is None
                color scale min value.
                None will take the min value found in the data.

            vmax: floating point number; default is None
                color scale max value.
                None will take the max value found in the data.

            alpha: floating point number; default is 1
                alpha channel transparency value
                between 1 (not transparent at all) and 0 (totally transparent).

            fill: boolean; default is True
                True generates a matplotlib contourf plot.
                False generates a matplotlib contour plot.

            cmap: string; default is viridis
                matplotlib color map color label.
                https://matplotlib.org/stable/tutorials/colors/colormaps.html

            title: string; default None
                possible plot title string.

            grid: boolean; default True
                should be plotted on  a grid or on a blank page?
                True will plot on a grid.

            xlim: tuple of two floating point numbers; default is None
                to specify min and max x axis value.
                None will extract agreeable values from the data.

            ylim: tuple of two floating point numbers; default is None
                to specify min and max y axis value.
                None will extract agreeable values from the data.

            xyequal: boolean; default True
                to specify equal axis spacing for x and y axis.

            figsize: tuple of floating point numbers; default is None
                the specif the figure x and y measurement in inch.
                None result in the default matplotlib setting, which is [6.4, 4.8].

            ax: matplotlib axis object; default setting is None
                the ax object, which will be used as a canvas for plotting.
                None will generate a figure and ax object from scratch.

        output:
            fig: matplotlib figure, containing the ax axis object,
                with contour plot and color bar.

        description:
            function returns a matplotlib contour (or contourf) plot,
            inclusive color bar, for the substrate specified.
        """
        # handle z_slice input
        _, _, ar_p_axis = self.get_mesh_mnp_axis()
        if not (z_slice in ar_p_axis):
            z_slice = ar_p_axis[abs(ar_p_axis - z_slice).argmin()]
            if self.verbose:
                print(f'z_slice set to {z_slice}.')

        # get data z slice
        df_conc = self.get_conc_df(values=1, drop=set(), keep=set())
        df_conc = df_conc.loc[(df_conc.mesh_center_p == z_slice),:]
        # extend to x y domain border
        df_mmin = df_conc.loc[(df_conc.mesh_center_m == df_conc.mesh_center_m.min()), :].copy()
        df_mmin.mesh_center_m = self.get_xyz_range()[0][0]
        df_mmax = df_conc.loc[(df_conc.mesh_center_m == df_conc.mesh_center_m.max()), :].copy()
        df_mmax.mesh_center_m = self.get_xyz_range()[0][1]
        df_conc = pd.concat([df_conc, df_mmin, df_mmax], axis=0)
        df_nmin = df_conc.loc[(df_conc.mesh_center_n == df_conc.mesh_center_n.min()), :].copy()
        df_nmin.mesh_center_n =self.get_xyz_range()[1][0]
        df_nmax = df_conc.loc[(df_conc.mesh_center_n == df_conc.mesh_center_n.max()), :].copy()
        df_nmax.mesh_center_n = self.get_xyz_range()[1][1]
        df_conc = pd.concat([df_conc, df_nmin, df_nmax], axis=0)
        # sort dataframe
        df_conc.sort_values(['mesh_center_m', 'mesh_center_n', 'mesh_center_p'], inplace=True)

        # meshgrid shape
        ti_shape = (self.get_voxel_ijk_axis()[0].shape[0]+2, self.get_voxel_ijk_axis()[1].shape[0]+2)
        x = (df_conc.loc[:,'mesh_center_m'].values).reshape(ti_shape)
        y = (df_conc.loc[:,'mesh_center_n'].values).reshape(ti_shape)
        z = (df_conc.loc[:,substrate].values).reshape(ti_shape)

        # handle vmin and vmax input
        if (vmin is None):
            vmin = np.floor(df_conc.loc[:,substrate].min())
        if (vmax is None):
            vmax = np.ceil(df_conc.loc[:,substrate].max())

        # get figure and axis orbject
        if (ax is None):
            # handle figsize
            if (figsize is None):
                figsize = (6.4, 4.8)
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = plt.gcf()

        # set equal axis spacing
        if xyequal:
            ax.axis('equal')

        # get contour plot
        if fill:
            ax.contourf(x,y,z, vmin=vmin, vmax=vmax, alpha=alpha, cmap=cmap)
        else:
            ax.contour(x,y,z, vmin=vmin, vmax=vmax, alpha=alpha, cmap=cmap)

        # set title
        if not (title is None):
            ax.set_title(title)

        # set grid
        ax.grid(visible=grid)

        # set axis lim
        if not (xlim is None):
            ax.set_xlim(xlim[0], xlim[1])
        if not (ylim is None):
            ax.set_ylim(ylim[0], ylim[1])

        # get colorbar
        fig.colorbar(
            mappable=cm.ScalarMappable(norm=colors.Normalize(vmin=vmin, vmax=vmax), cmap=cmap),
            label=substrate,
            ax=ax
        )

        # output
        return fig


    def make_conc_vtk(self):
        """
        input:

        output:
            vtk: vtk file that contains 3D distributions of all substrates
                over microenvironment with corresponding time stamp.

        description:
            function creates rectilinear grid vtk file contains distribution of
            substrates over microenvironment. You can post-process this file
            in other software like Paraview.
        """
        # Get microenviornment data frame
        df_micenv = self.get_conc_df()

        # Create a rectilinear grid
        vr_grid = vtk.vtkRectilinearGrid()

        # Define dimensions of the grid
        t_dims = (len(pd.unique(df_micenv['voxel_i'])), len(pd.unique(df_micenv['voxel_j'])), len(pd.unique(df_micenv['voxel_k'])))

        # Define coordinates for the grid
        vf_x = vtk.vtkFloatArray()
        vf_y = vtk.vtkFloatArray()
        vf_z = vtk.vtkFloatArray()

        # Assign coordinates for the grid
        vf_x.SetNumberOfTuples(t_dims[0])
        vf_y.SetNumberOfTuples(t_dims[1])
        vf_z.SetNumberOfTuples(t_dims[2])

        # Populate the coordinates
        for i in range(t_dims[0]):
            vf_x.SetValue(i, pd.unique(df_micenv['mesh_center_m'])[i])

        for j in range(t_dims[1]):
            vf_y.SetValue(j, pd.unique(df_micenv['mesh_center_n'])[j])

        for k in range(t_dims[2]):
            vf_z.SetValue(k, pd.unique(df_micenv['mesh_center_p'])[k])

        # Grid dimensions
        vr_grid.SetDimensions(t_dims)
        vr_grid.SetXCoordinates(vf_x)
        vr_grid.SetYCoordinates(vf_y)
        vr_grid.SetZCoordinates(vf_z)

        # Learn substrate names
        l_substrate_names = self.get_substrate_names()

        # For loop to fill rectilinear grid
        for name_index, name in enumerate(l_substrate_names):
            vf_values = vtk.vtkFloatArray()
            vf_values.SetNumberOfComponents(1)
            vf_values.SetNumberOfTuples(t_dims[0] * t_dims[1] * t_dims[2])
            vf_values.SetName(name)  # Set the name of the array
            # Populate the substrate values
            for k in range(t_dims[2]):
                for j in range(t_dims[1]):
                    for i in range(t_dims[0]):
                        idx = i + t_dims[0] * (j + t_dims[1] * k)
                        conc_at_specific_position = df_micenv[name][idx]
                        vf_values.SetValue(idx, conc_at_specific_position)
            if name_index == 0:
                vr_grid.GetPointData().SetScalars(vf_values)
            else:
                vr_grid.GetPointData().AddArray(vf_values)
            del vf_values


        # Save vtk file
        s_vtkpathfile = self.path + '/' + self.xmlfile.replace('.xml','_conc.vtk')
        vw_writer = vtk.vtkXMLRectilinearGridWriter()
        vw_writer.SetFileName(s_vtkpathfile)
        vw_writer.SetInputData(vr_grid)
        vw_writer.Write()
        return s_vtkpathfile


    ## CELL RELATED FUNCTIONS ##

    def get_cell_variables(self):
        """
        input:

        output:
            ls_variables: list of strings
                alphabetically ordered list of all tracked cell variable names.

        description:
            function returns all modeled cell variable names.
        """
        ls_variables = sorted(self.data['discrete_cells']['data'].keys())
        return ls_variables


    def get_celltype_dict(self):
        """
        input:

        output:
            ds_celltype: dictionary of stings
                dictionary that maps cell_type IDs to labels.

        description:
            function returns a dictionary that maps ID and name from all
            cell_definitions, specified in the PhysiCell_settings.xml file.
        """
        return self.data['metadata']['cell_type']


    def get_cell_df(self, values=1, drop=set(), keep=set()):
        """
        input:
            values: integer; default is 1
                minimal number of values a variable has to have to be outputted.
                variables that have only 1 state carry no information.
                None is a state too.

            drop: set of strings; default is an empty set
                set of column labels to be dropped for the dataframe.
                don't worry: essential columns like ID, coordinates
                and time will never be dropped.
                Attention: when the keep parameter is given, then
                the drop parameter has to be an empty set!

            keep: set of strings; default is an empty set
                set of column labels to be kept in the dataframe.
                set values=1 to be sure that all variables are kept.
                don't worry: essential columns like ID, coordinates,
                time and runtime (wall time) will always be kept.

        output:
            df_cell: pandas dataframe
                dataframe lists, one cell per row, all tracked variables
                values related to this cell. the variables are cell_position,
                mesh_center, and voxel coordinates, all cell_variables,
                all substrate rates and concentrations, and additional
                the surrounding cell density.

        description:
            function returns a dataframe with a cell centric view
            of the simulation.
        """
        # check keep and drop
        if (len(keep) > 0) and (len(drop) > 0):
            sys.exit(f"Error @ pyMCDS.get_cell_df : when keep is given {keep}, then drop has to be an empty set {drop}!")

        # get cell position and more
        df_cell = pd.DataFrame(self.data['discrete_cells']['data'])
        df_cell['time'] = self.get_time()
        df_cell['runtime'] = self.get_runtime() / 60  # in min
        df_cell['xmlfile'] = self.xmlfile
        df_voxel = df_cell.loc[:,['position_x','position_y','position_z']].copy()

        # get mesh spacing
        dm, dn, dp = self.get_voxel_spacing()

        # get mesh and voxel min max values
        tr_m_range, tr_n_range, tr_p_range = self.get_mesh_mnp_range()
        tr_i_range, tr_j_range, tr_k_range = self.get_voxel_ijk_range()

        # get voxel for each cell
        df_voxel.loc[:,'voxel_i'] = np.round((df_voxel.loc[:,'position_x'] - tr_m_range[0]) / dm).astype(int)
        df_voxel.loc[:,'voxel_j'] = np.round((df_voxel.loc[:,'position_y'] - tr_n_range[0]) / dn).astype(int)
        df_voxel.loc[:,'voxel_k'] = np.round((df_voxel.loc[:,'position_z'] - tr_p_range[0]) / dp).astype(int)
        df_voxel.loc[(df_voxel.voxel_i > tr_i_range[1]), 'voxel_i'] = tr_i_range[1]  # i_max
        df_voxel.loc[(df_voxel.voxel_i < tr_i_range[0]), 'voxel_i'] = tr_i_range[0]  # i_min
        df_voxel.loc[(df_voxel.voxel_j > tr_j_range[1]), 'voxel_j'] = tr_j_range[1]  # j_max
        df_voxel.loc[(df_voxel.voxel_j < tr_j_range[0]), 'voxel_j'] = tr_j_range[0]  # j_min
        df_voxel.loc[(df_voxel.voxel_k > tr_k_range[1]), 'voxel_k'] = tr_k_range[1]  # k_max
        df_voxel.loc[(df_voxel.voxel_k < tr_k_range[0]), 'voxel_k'] = tr_k_range[0]  # k_min

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

        # get column label set
        es_column = set(df_cell.columns)

        # get vector length
        for s_var_spatial in es_var_spatial:
            es_vector = es_column.intersection({f'{s_var_spatial}_x',f'{s_var_spatial}_y',f'{s_var_spatial}_z'})
            if len(es_vector) > 0:
                # linear algebra
                #a_vector = df_cell.loc[:,ls_vector].values
                #a_length = np.sqrt(np.diag(np.dot(a_vector, a_vector.T)))
                # pythoagoras
                a_length = None
                for s_vector in es_vector:
                    a_vectorsq = df_cell.loc[:,s_vector].values**2
                    if (a_length is None):
                        a_length = a_vectorsq
                    else:
                        a_length += a_vectorsq
                a_length = a_length**(1/2)
                # result
                df_cell[f'{s_var_spatial}_vectorlength'] = a_length

        # physicell
        if not (self.data['discrete_cells']['physiboss'] is None):
            df_cell = pd.merge(
                df_cell,
                self.data['discrete_cells']['physiboss'],
                left_index = True,
                right_index = True,
                how = 'left',
            )


        # microenvironment
        if self.microenv:
            # merge substrate (left join)
            df_sub = self.get_substrate_df()
            for s_sub in df_sub.index:
                 for s_rate in df_sub.columns:
                     s_var = f'{s_sub}_{s_rate}'
                     df_cell[s_var] = df_sub.loc[s_sub,s_rate]

        # merge concentration (left join)
        df_conc = self.get_conc_df(z_slice=None, values=1, drop=set(), keep=set())
        df_conc.drop({'time', 'runtime','xmlfile'}, axis=1, inplace=True)
        df_cell = pd.merge(
            df_cell,
            df_conc,
            on = ['voxel_i', 'voxel_j', 'voxel_k'],
            how = 'left',
        )

        # variable typing
        do_type = {}
        [do_type.update({k:v}) for k,v in do_var_type.items() if k in es_column]
        do_type.update(self.custom_type)
        do_int = do_type.copy()
        [do_int.update({k:int}) for k in do_int.keys()]
        ls_int = sorted(do_int.keys())
        df_cell.loc[:,ls_int] = df_cell.loc[:,ls_int].round()
        df_cell = df_cell.astype(do_int)
        df_cell = df_cell.astype(do_type)

        # categorical translation
        df_cell.loc[:,'current_death_model'] = df_cell.loc[:,'current_death_model'].replace(ds_death_model)  # bue 20230614: this column looks like an artefact to me
        df_cell.loc[:,'cycle_model'] = df_cell.loc[:,'cycle_model'].replace(ds_cycle_model)
        df_cell.loc[:,'cycle_model'] = df_cell.loc[:,'cycle_model'].replace(ds_death_model)
        df_cell.loc[:,'current_phase'] = df_cell.loc[:,'current_phase'].replace(ds_cycle_phase)
        df_cell.loc[:,'current_phase'] = df_cell.loc[:,'current_phase'].replace(ds_death_phase)
        df_cell.loc[:,'cell_type'] = df_cell.loc[:,'cell_type'].replace(self.data['metadata']['cell_type'])

        # filter
        es_feature = set(df_cell.columns).difference(es_coor_cell)
        if (len(keep) > 0):
            es_delete = es_feature.difference(keep)
        else:
            es_delete = es_feature.intersection(drop)

        if (values > 1):  # by minimal number of states
            for s_column in set(df_cell.columns).difference(es_coor_cell):
                if len(set(df_cell.loc[:,s_column])) < values:
                    es_delete.add(s_column)
        if self.verbose and (len(es_delete) > 0):
            print('es_delete:', es_delete)
        df_cell.drop(es_delete, axis=1, inplace=True)

        # output
        df_cell = df_cell.loc[:,sorted(df_cell.columns)]
        df_cell.set_index('ID', inplace=True)
        df_cell = df_cell.copy()
        return df_cell


    def get_cell_df_at(self, x, y, z=0, values=1, drop=set(), keep=set()):
        """
        input:
            x: floating point number
                position x-coordinate.

            y: floating point number
                position y-coordinate.

            z: floating point number; default is 0
                position z-coordinate.

            values: integer; default is 1
                minimal number of values a variable has to have to be outputted.
                variables that have only 1 state carry no information.
                None is a state too.

            drop: set of strings; default is an empty set
                set of column labels to be dropped for the dataframe.
                don't worry: essential columns like ID, coordinates
                and time will never be dropped.
                Attention: when the keep parameter is given, then
                the drop parameter has to be an empty set!

            keep: set of strings; default is an empty set
                set of column labels to be kept in the dataframe.
                set values=1 to be sure that all variables are kept.
                don't worry: essential columns like ID, coordinates
                and time will always be kept.

        output:
            df_voxel: pandas dataframe
                x, y, z voxel filtered cell dataframe.

        description:
            function returns the cell dataframe for the voxel
            specified with the x, y, z position coordinate.
        """
        df_voxel = None

        # is coordinate inside the domain?
        b_calc = self.is_in_mesh(x=x, y=y, z=z, halt=False)
        if b_calc:

            # get mesh and mesh spacing
            dm, dn, dp = self.get_voxel_spacing()
            ar_m, ar_n, ar_p = self.get_mesh()

            # get voxel coordinate
            i, j, k = self.get_voxel_ijk(x, y, z, is_in_mesh=False)
            m = ar_m[j, i, k]
            n = ar_n[j, i, k]
            p = ar_p[j, i, k]

            # get voxel
            df_cell = self.get_cell_df(values=values, drop=drop, keep=keep)
            inside_voxel = (
                (df_cell['position_x'] <= m + dm / 2) &
                (df_cell['position_x'] >= m - dm / 2) &
                (df_cell['position_y'] <= n + dn / 2) &
                (df_cell['position_y'] >= n - dn / 2) &
                (df_cell['position_z'] <= p + dp / 2) &
                (df_cell['position_z'] >= p - dp / 2)
            )
            df_voxel = df_cell[inside_voxel]

        # output
        return df_voxel


    def plot_scatter(self, focus='cell_type', z_slice=0, z_axis=None, alpha=1, cmap='viridis', title=None, grid=True, legend_loc='lower left', xlim=None, ylim=None, xyequal=True, s=None, figsize=None, ax=None):
        """
        input:
            focus: string; default is 'cell_type'
                column name within cell dataframe.

            z_slice: floating point number; default is 0
                z-axis position to slice a 2D xy-plain out of the
                3D substrate concentration mesh. if z_slice position
                is not an exact mesh center coordinate, then z_slice
                will be adjusted to the nearest mesh center value,
                the smaller one, if the coordinate lies on a saddle point.

            z_axis: for a categorical focus: set of labels;
               for a numeric focus: tuple of two floats; default is None
               depending on the focus column variable dtype, default extracts
               labels or min and max values from data.

            alpha: floating point number; default is 1
                alpha channel transparency value
                between 1 (not transparent at all) and 0 (totally transparent).

            cmap: dictionary of strings or string; default viridis.
                dictionary that maps labels to colors strings.
                matplotlib colormap string.
                https://matplotlib.org/stable/tutorials/colors/colormaps.html

            title: string; default None
                possible plot title string.

            grid: boolean default True.
                plot axis grid lines.

            legend_loc: string; default is 'lower left'.
                the location of the categorical legend, if applicable.
                possible strings are: best,
                upper right, upper center, upper left, center left,
                lower left, lower center, lower right, center right,
                center.

            xlim: tuple of two floats; default is None
                x axis min and max value.
                default takes min and max from mesh x axis range.

            ylim: tuple of two floats; default is None
                y axis min and max value.
                default takes min and max from mesh y axis range.

           xyequal: boolean; default True
                to specify equal axis spacing for x and y axis.

            s: integer; default is None
                scatter plot dot size in pixel.
                typographic points are 1/72 inch.
                the marker size s is specified in points**2.
                plt.rcParams['lines.markersize']**2 is in my case 36.
                None tries to take the value from the initial.svg file.
                fall back setting is 36.

            figsize: tuple of floating point numbers; default is None
                the specif the figure x and y measurement in inch.
                None result in the default matplotlib setting, which is [6.4, 4.8].

        output:
            fig: matplotlib figure, containing the ax axis object,
                with scatter plot and color bar (numerical data)
                or color legend (categorical data).

        description:
            function returns a (pandas) matplotlib scatter plot,
            inclusive color bar, for the substrate specified.
        """
        # handle z_slice
        _, _, ar_p_axis = self.get_mesh_mnp_axis()
        if not (z_slice in ar_p_axis):
            z_slice = ar_p_axis[abs(ar_p_axis - z_slice).argmin()]
            if self.verbose:
                print(f'z_slice set to {z_slice}.')

        # get data z slice
        df_cell = self.get_cell_df(values=1, drop=set(), keep=set())
        df_cell = df_cell.loc[(df_cell.mesh_center_p == z_slice),:]

        # handle z_axis categorical cases
        if (str(df_cell.loc[:,focus].dtype) in {'bool', 'object'}):
            lr_extrema = [None, None]
            if (z_axis is None):
                # extract set of labels from data
                es_category = set(df_cell.loc[:,focus])
                if (str(df_cell.loc[:,focus].dtype) in {'bool'}):
                    es_category = es_category.union({True, False})
            else:
                es_category = z_axis

        # handle z_axis numerical cases
        else:  # df_cell.loc[:,focus].dtype is numeric
            es_category = None
            if (z_axis is None):
                # extract min and max values from data
                r_min = df_cell.loc[:,focus].min()
                r_max = df_cell.loc[:,focus].max()
                lr_extrema = [r_min, r_max]
            else:
                lr_extrema = z_axis

        # handle z_axis summary
        if self.verbose:
            print(f'categories found: {es_category}.')
            print(f'min max extrema set to: {lr_extrema}.')

        # handle xlim and ylim
        if (xlim is None):
            xlim = self.get_xyz_range()[0]
            if self.verbose:
                print(f'xlim set to: {xlim}.')
        if (ylim is None):
            ylim = self.get_xyz_range()[1]
            if self.verbose:
                print(f'ylim set to: {ylim}.')

        # get figure and axis orbject
        if (ax is None):
            # handle figsize
            if (figsize is None):
                figsize = (6.4, 4.8)
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = plt.gcf()

        # layout the canavas
        if xyequal:
            ax.axis('equal')

        # handle categorical variable
        if not (es_category is None):
            s_focus_color = focus + '_color'
            # use specified category color dictionary
            if type(cmap) is dict:
                ds_color = cmap
                df_cell[s_focus_color] = 'gray'
                for s_category, s_color in ds_color.items():
                    df_cell.loc[(df_cell.loc[:,focus] == s_category), s_focus_color] = s_color
            # generate category color dictionary
            else:
                ds_color = pdplt.df_label_to_color(
                    df_abc = df_cell,
                    s_focus = focus,
                    es_label = es_category,
                    s_nolabel = 'gray',
                    s_cmap = cmap,
                    b_shuffle = False,
                )
            # generate color list
            c = list(df_cell.loc[:, s_focus_color].values)
            s_cmap = None

        # handle numeric variable
        else:
            c = focus
            s_cmap = cmap

        # plot scatter
        df_cell.plot(
            kind = 'scatter',
            x = 'position_x',
            y = 'position_y',
            c = c,
            vmin = lr_extrema[0],
            vmax = lr_extrema[1],
            alpha = alpha,
            cmap = s_cmap,
            title = title,
            xlim = xlim,
            ylim = ylim,
            s = s,
            grid = grid,
            ax = ax,
        )

        # plot categorical data legen
        if not (es_category is None):
            pdplt.ax_colorlegend(
                ax = ax,
                ds_color = ds_color,
                s_loc = legend_loc,
                s_fontsize = 'small',
            )

        # output
        return fig


    def make_cell_vtk(self, l_attributes=['cell_type'], visualize='True'):
        """
        input:

            attributes: list of strings; default is 'cell_type'
                column name within cell dataframe.

            visualize: boolean; default is True
                visualize cells using vtk Renderer.

        output:
            vtk: 3D Glyph VTK that contains cells.


        description:
            function that 3D Glyph VTK file for cells. Cells can have
            specificed attributes like 'cell_type', 'pressure', 'dead', etc.
            You can post-process this file in other software like Paraview.
        """
        # Get cell data frame
        df_cell = self.get_cell_df(values=1, drop=set(), keep=set())
        df_cell = df_cell.reset_index()

        # Get positions and radii
        se_x_pos = df_cell['position_x']
        se_y_pos = df_cell['position_y']
        se_z_pos = df_cell['position_z']
        se_radius =  df_cell['radius']

        # Create VTK instances to fill for positions and radii
        vp_points = vtkPoints()
        vf_radii = vtk.vtkFloatArray()
        vf_radii.SetName("radius")

        # Fill VTK instance with positions and radii
        for i in range(len(se_x_pos)):
            vp_points.InsertNextPoint(se_x_pos[i], se_y_pos[i], se_z_pos[i])
            vf_radii.InsertNextValue(se_radius[i])


        # Create Data Instances
        vf_data = vtk.vtkFloatArray()
        vf_data.SetNumberOfComponents(2)
        vf_data.SetNumberOfTuples(len(se_x_pos))
        vf_data.CopyComponent(0, vf_radii, 0)
        vf_data.SetName("positions_and_radii")

        # Create Unstructred Grid for Data
        vu_grid = vtk.vtkUnstructuredGrid()
        vu_grid.SetPoints(vp_points)
        vu_grid.GetPointData().AddArray(vf_data)
        vu_grid.GetPointData().SetActiveScalars("positions_and_radii")

        # Fill This grid with given attributes
        for name_index, name in enumerate(l_attributes):
            custom_data_i = df_cell[name]
            if (type(custom_data_i[0]) in {str, np.str_}):
                custom_data_vtk = vtk.vtkStringArray()
                custom_data_vtk.SetName(name)
            elif(type(custom_data_i[0]) in {float, np.float_, np.float16, np.float32, np.float64, int}):
                custom_data_vtk = vtk.vtkFloatArray()
                custom_data_vtk.SetName(name)
            elif (type(custom_data_i[0]) in {bool, np.bool_}):
                custom_data_vtk =vtk.vtkStringArray()
                custom_data_vtk.SetName(name)

            for i in range(len(custom_data_i)):
                if (type(custom_data_i[0]) in {bool, np.bool_}):
                    if (custom_data_i[i] == True):
                        custom_data_vtk.InsertNextValue('True')
                    else:
                        custom_data_vtk.InsertNextValue('False')
                else:
                    custom_data_vtk.InsertNextValue(custom_data_i[i])

            vu_grid.GetPointData().AddArray(custom_data_vtk)
            del custom_data_vtk


        # Create sphere source
        vsp_sphere = vtk.vtkSphereSource()
        vsp_sphere.SetRadius(1.0)
        vsp_sphere.SetPhiResolution(16)
        vsp_sphere.SetThetaResolution(32)

        # Create Glyph to save
        vg_glyph = vtk.vtkGlyph3D()
        vg_glyph.SetInputData(vu_grid)
        vg_glyph.SetSourceConnection(vsp_sphere.GetOutputPort())

        # Define important preferences for VTK
        vg_glyph.ClampingOff()
        vg_glyph.SetScaleModeToScaleByScalar()
        vg_glyph.SetScaleFactor(1.0)
        vg_glyph.SetColorModeToColorByScalar()
        vg_glyph.Update()

        # Write VTK
        s_vtkpathfile = self.path + '/' + self.xmlfile.replace('.xml','_cells.vtk')
        vw_writer = vtk.vtkXMLPolyDataWriter()
        vw_writer.SetFileName(s_vtkpathfile)
        vw_writer.SetInputData(vg_glyph.GetOutput())
        vw_writer.Write()


        # Visualize if needed
        if (visualize):
            # visualization
            # set up the mapper
            mapper = vtk.vtkPolyDataMapper()
            # mapper.SetInput(glyph.GetOutput())
            mapper.SetInputConnection(vg_glyph.GetOutputPort())

            mapper.ScalarVisibilityOn()
            mapper.ColorByArrayComponent("data", 1)

            # set up the actor
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)

            # do renderer setup stuff
            ren = vtk.vtkRenderer()
            renWin = vtk.vtkRenderWindow()
            renWin.AddRenderer(ren)
            renWin.SetSize(640, 480)
            iren = vtk.vtkRenderWindowInteractor()
            iren.SetRenderWindow(renWin)

            # add the actor to the renderer
            ren.AddActor(actor)

            # render
            iren.Initialize()
            renWin.Render()
            iren.Start()

        return s_vtkpathfile


    ## GRAPH RELATED FUNCTIONS ##

    def get_attached_graph_dict(self):
        """
        input:

        output:
            dei_graph: dictionary of sets of integers
                maps each cell ID to the attached connected cell IDs.

        description:
            function returns the attached cell graph as a dictionary object.
        """
        return self.data['discrete_cells']['graph']['attached_cells']


    def get_neighbor_graph_dict(self):
        """
        input:

        output:
            dei_graph: dictionary of sets of integers
                maps each cell ID to the connected neighbor cell IDs.

        description:
            function returns the cell neighbor graph as a dictionary object.
        """
        return self.data['discrete_cells']['graph']['neighbor_cells']


    def make_graph_gml(self, graph_type='neighbor', edge_attr=True, node_attr=[]):
        """
        input:
            graph_type: string; default is neighbor
                to specify which physicell output data should be processed.
                attached: processes mcds.get_attached_graph_dict dictionary.
                neighbor: processes mcds.get_neighbor_graph_dict dictionary.

            edge_attr: boolean; default True
                specifies if the spatial Euclidean distance is used for
                edge attribute, to generate a weighted graph.

            node_attr: list of strings; default is empty list
                list of mcds.get_cell_df dataframe columns, used for
                node attributes.

        output:
            gml file, generated under the returned path.

        description:
            function to generate graph files in the gml graph modelling language
            standard format.

            gml was the outcome of an initiative that started at
            the international symposium on graph drawing 1995 in Passau
            and ended at Graph Drawing 1996 in Berkeley. the networkx python
            and igraph C and python libraries for graph analysis are
            gml compatible and can as such read and write this file format.

            https://en.wikipedia.org/wiki/Graph_Modelling_Language
            https://networkx.org/
            https://igraph.org/
        """
        # load dataframe for celltype information
        df_cell = self.get_cell_df()
        ds_unit = self.get_unit_dict()
        s_unit_simtime = ds_unit["time"]
        r_simtime = self.get_time()
        if (graph_type in {'attached'}):
            dei_graph = self.get_attached_graph_dict()
        elif (graph_type in {'neighbor'}):
            dei_graph = self.get_neighbor_graph_dict()
        else:
            sys.exit(f'Erro @ make_graph_gml : unknowen graph_type {graph_type}. knowen are attached and neighbor.')

        # generate filename
        s_gmlpathfile = self.path + '/' + self.xmlfile.replace('.xml',f'_{graph_type}.gml')

        # open result gml file
        f = open(s_gmlpathfile, 'w')
        f.write(f'Creator "pcdl_v{__version__}"\ngraph [\n')
        f.write(f'  id {int(r_simtime)}\n  comment "time_{s_unit_simtime}"\n  label "{graph_type}_graph"\n  directed 0\n')
        for i_src, ei_dst in dei_graph.items():
            #print(f'{i_src} {sorted(ei_dst)}')
            # node
            f.write(f'  node [\n    id {i_src}\n    label "node_{i_src}"\n')
            # node attributes
            for s_attr in node_attr:
                o_attr = df_cell.loc[i_src, s_attr]
                if (type(o_attr) in {bool, np.bool_, int, np.int_, np.int8, np.int16, np.int32, np.int64}):
                    f.write(f'    {s_attr} {int(o_attr)}\n')
                elif (type(o_attr) in {float, np.float_, np.float16, np.float32, np.float64}):  # np.float128
                    f.write(f'    {s_attr} {o_attr}\n')
                elif (type(o_attr) in {str, np.str_}):
                    f.write(f'    {s_attr} "{o_attr}"\n')
                else:
                    sys.exit(f'Error @ make_graph_gml : attr {o_attr}; type {type (o_attr)}; type seems not to be bool, int, float, or string.')
            f.write(f'  ]\n')
            # edge
            for i_dst in ei_dst:
                if (i_src < i_dst):
                    f.write(f'  edge [\n    source {i_src}\n    target {i_dst}\n    label "edge_{i_src}_{i_dst}"\n')
                    if (edge_attr):
                        # edge distance attribute
                        x = df_cell.loc[i_src, 'position_x'] - df_cell.loc[i_dst, 'position_x']
                        y = df_cell.loc[i_src, 'position_y'] - df_cell.loc[i_dst, 'position_y']
                        z = df_cell.loc[i_src, 'position_z'] - df_cell.loc[i_dst, 'position_z']
                        r_distance = (x**2 + y**2 + z**2)**(1/2)
                        f.write(f'    distance_{ds_unit["position_y"]} {round(r_distance)}\n')
                    f.write(f'  ]\n')
            # development
            #if (i_src > 16):
            #    break
        # close result gml file
        f.write(']\n')
        f.close()

        # output
        return s_gmlpathfile


    ## MODEL PARAMETER SETTING RELATED FUNCTIONS ##

    def get_parameter_dict(self):
        """
        input:

        output:
            d_parameter: dictionary
            dictionary, mapping parameter and input value.

        description:
            function retunes a dictionary that maps the models
            input parameters and values fund in the settings.xml file.
            only parameters compatible with the PhysiCell Studio settings.xml
            version are listed. other parameters might be missing!
        """
        return self.data['setting']['parameters']


    def get_rule_df(self):
        """
        input:

        output:
            df_rules: pandas dataframe
            rules.csv loaded as datafarme.

        description:
            function returns the rule csv, linked in the settings.xml file,
            as a datafram.
        """
        return self.data['setting']['rules']


    def get_unit_dict(self):
        """
        input:

        output:
            ds_unit: dictionary
                dictionary, which tracks units from cell and microenvironment
                variables and settings parameters.

        description:
            function returns a dictionary that stores all tracked variables
            and their units and all parameters and their units found in the
            settings.xml file.
        """
        return self.data['setting']['units']


    ## LOAD DATA  ##
    def _read_xml(self, xmlfile, output_path='.'):
        """
        input:
            self: pyMCDS class instance.

            xmlfile: string
                name of the xml file with or without path
                in the with path case, output_path has to be set to the default!

            output_path: string; default '.'
                relative or absolute path to the directory where
                the PhysiCell output files are stored.

        output:
            self: pyMCDS class instance with loaded data.

        description:
            internal function to load the data from the PhysiCell output files
            into the pyMCDS instance.
        """
        #####################
        # path and filename #
        #####################

        # file and path manipulation
        s_xmlfile = xmlfile.replace('\\','/')
        if (xmlfile.find('/') > -1) and (output_path == '.'):
            ls_xmlfile = xmlfile.split('/')
            s_xmlfile = ls_xmlfile.pop(-1)
            output_path = '/'.join(ls_xmlfile)
        self.path = output_path
        self.xmlfile = s_xmlfile

        # generate output dictionary
        d_mcds = {}
        d_mcds['metadata'] = {}
        d_mcds['metadata']['substrate'] = {}
        d_mcds['metadata']['cell_type'] = {}
        d_mcds['setting'] = {}
        d_mcds['setting']['parameters'] = {}
        d_mcds['setting']['units'] = {}
        d_mcds['setting']['rules'] = None  # df_ruleset
        d_mcds['mesh'] = {}
        d_mcds['continuum_variables'] = {}
        d_mcds['discrete_cells'] = {}
        d_mcds['discrete_cells']['units'] = {}


        #####################################
        # PhysiCell_settings.xml extraction #
        #####################################
        # bue 2024-03-11: this part tries to be compatible with the PhysiCell Studio settings.xml file only,
        # older and other settings.xml versions are disregarded.

        if not ((self.settingxml is None) or (self.settingxml is False)):
            # load Physicell_settings xml file
            s_xmlpathfile_setting = self.path + '/' + self.settingxml
            x_tree = ET.parse(s_xmlpathfile_setting)
            if self.verbose:
                print(f'reading: {s_xmlpathfile_setting}')
            x_root = x_tree.getroot()

            # skip <domain> node for mesh

            # find <overall> node for dt
            x_overall = x_root.find('overall')
            if not (x_overall is None):
                # dt_diffusion
                try:
                    d_mcds['setting']['parameters'].update({'dt_diffusion': float(x_overall.find('dt_diffusion').text)})
                    d_mcds['setting']['units'].update({'dt_diffusion': x_overall.find('dt_diffusion').get('units')})
                except AttributeError:
                    print(f'Warning @ pyMCDS._read_setting_xml : <overall><dt_diffusion> node missing.')
                # dt_mechanics
                try:
                    d_mcds['setting']['parameters'].update({'dt_mechanics': float(x_overall.find('dt_mechanics').text)})
                    d_mcds['setting']['units'].update({'dt_mechanics': x_overall.find('dt_mechanics').get('units')})
                except AttributeError:
                    print(f'Warning @ pyMCDS._read_setting_xml : <overall><dt_mechanics> node missing.')

                # dt_phenotype
                try:
                    d_mcds['setting']['parameters'].update({'dt_phenotype': float(x_overall.find('dt_phenotype').text)})
                    d_mcds['setting']['units'].update({'dt_phenotype': x_overall.find('dt_phenotype').get('units')})
                except AttributeError:
                    print(f'Warning @ pyMCDS._read_setting_xml : <overall><dt_phenotype> node missing.')
            else:
                print(f'Warning @ pyMCDS._read_setting_xml : <overall> node missing.')

            # skip <parallel> node
            # skip <save> node

            # find <options> node
            x_options = x_root.find('options')
            if not (x_options is None):
                try:
                    d_mcds['setting']['parameters'].update({'legacy_random_points_on_sphere_in_divide': str(x_options.find('legacy_random_points_on_sphere_in_divide').text).lower == 'true'})
                except AttributeError:
                    print(f'Warning @ pyMCDS._read_setting_xml : <options><legacy_random_points_on_sphere_in_divide> node missing.')
                try:
                    d_mcds['setting']['parameters'].update({'virtual_wall_at_domain_edge': str(x_options.find('virtual_wall_at_domain_edge').text).lower == 'true'})
                except AttributeError:
                    print(f'Warning @ pyMCDS._read_setting_xml : <options><virtual_wall_at_domain_edge> node missing.')
                try:
                    d_mcds['setting']['parameters'].update({'disable_automated_spring_adhesions': str(x_options.find('disable_automated_spring_adhesions').text).lower == 'true'})
                except AttributeError:
                    print(f'Warning @ pyMCDS._read_setting_xml : <options><disable_automated_spring_adhesions> node missing.')
            else:
                print(f'Warning @ pyMCDS._read_setting_xml : <options> node missing.')

            # find <microenvironment_setup> node
            x_microenvironment = x_root.find('microenvironment_setup')
            # substrate loop
            for x_variable in x_microenvironment.findall('variable'):

                # <variable>
                s_id = str(x_variable.get('ID'))
                s_substrate = x_variable.get('name').replace(' ', '_')
                d_mcds['metadata']['substrate'].update({s_id : s_substrate})

                # <physical_parameter_set>
                x_microenvironment_physical = x_variable.find('physical_parameter_set')
                if not (x_microenvironment_physical is None):
                    # <diffusion_coefficient>
                    try:
                        d_mcds['setting']['parameters'].update({f'{s_substrate}_diffusion_coefficient': float(x_microenvironment_physical.find('diffusion_coefficient').text)})
                        d_mcds['setting']['units'].update({f'{s_substrate}_diffusion_coefficient': x_microenvironment_physical.find('diffusion_coefficient').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <variable name="{s_substrate}" ID="{s_id}"><physical_parameter_set><diffusion_coefficient> node missing.')
                    # <decay_rate>
                    try:
                        d_mcds['setting']['parameters'].update({f'{s_substrate}_decay_rate': float(x_microenvironment_physical.find('decay_rate').text)})
                        d_mcds['setting']['units'].update({f'{s_substrate}_decay_rate': x_microenvironment_physical.find('decay_rate').get('units')})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <variable name="{s_substrate}" ID="{s_id}"><physical_parameter_set><decay_rate> node missing.')
                else:
                    print(f'Warning @ pyMCDS._read_setting_xml : <variable name="{s_substrate}" ID="{s_id}"><physical_parameter_set> node missing.')
                # <initial_condition>
                try:
                    d_mcds['setting']['parameters'].update({f'{s_substrate}_initial_condition': float(x_variable.find('initial_condition').text)})
                    d_mcds['setting']['units'].update({f'{s_substrate}_initial_condition': x_variable.find('initial_condition').get('units')})
                except AttributeError:
                    print(f'Warning @ pyMCDS._read_setting_xml : <variable name="{s_substrate}" ID="{s_id}"><initial_condition> node missing.')
                # <Dirichlet_boundary_condition>
                try:
                    d_mcds['setting']['parameters'].update({f'{s_substrate}_dirichlet_boundary_condition_enabled': str(x_variable.find('Dirichlet_boundary_condition').get('enabled')).lower == 'true'})
                    d_mcds['setting']['parameters'].update({f'{s_substrate}_dirichlet_boundary_condition': float(x_variable.find('Dirichlet_boundary_condition').text)})
                    d_mcds['setting']['units'].update({f'{s_substrate}_dirichlet_boundary_condition': x_variable.find('Dirichlet_boundary_condition').get('units')})
                except AttributeError:
                    print(f'Warning @ pyMCDS._read_setting_xml : <variable name="{s_substrate}" ID="{s_id}"><Dirichlet_boundary_condition> node missing.')
                # <Dirichlet_options>
                x_microenvironment_dirichlet = x_variable.find('Dirichlet_options')
                if not (x_microenvironment_dirichlet is None):
                    for x_dirichlet in x_microenvironment_dirichlet:
                        try:
                            s_coor = x_dirichlet.get('ID')
                            d_mcds['setting']['parameters'].update({f'{s_substrate}_dirichlet_boundary_value_{s_coor}_enabled': str(x_dirichlet.text).lower == 'true'})
                            d_mcds['setting']['parameters'].update({f'{s_substrate}_dirichlet_boundary_value_{s_coor}': float(x_dirichlet.text)})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <variable name="{s_substrate}" ID="{s_id}"><Dirichlet_options></boundary_value> node missing.')
                        except TypeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <variable name="{s_substrate}" ID="{s_id}"><Dirichlet_options></boundary_value> node value missing.')
                else:
                    print(f'Warning @ pyMCDS._read_setting_xml : <variable name="{s_substrate}" ID="{s_id}"><Dirichlet_options> node missing.')
                # <options>
                x_microenvironment_options = x_microenvironment.find('options')
                if not (x_microenvironment_options is None):
                    # <options><calculate_gradients>
                    try:
                        d_mcds['setting']['parameters'].update({f'calculate_gradients': str(x_microenvironment_options.find('calculate_gradients').text).lower == 'true'})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <variable name="{s_substrate}" ID="{s_id}"><options><calculate_gradients> node missing.')
                    # <options><track_internalized_substrates_in_each_agent>
                    try:
                        d_mcds['setting']['parameters'].update({f'track_internalized_substrates_in_each_agent': str(x_microenvironment_options.find('track_internalized_substrates_in_each_agent').text).lower == 'true'})
                    except AttributeError:
                        print(f'Warning @ pyMCDS._read_setting_xml : <variable name="{s_substrate}" ID="{s_id}"><options><calculate_gradients> node missing.')
                    # skip <options><initial_condition>
                    # skip <options><dirichlet_nodes>
                else:
                    print(f'Warning @ pyMCDS._read_setting_xml : <variable name="{s_substrate}" ID="{s_id}"><microenvironment_setup><options> node missing.')

            # find <cell_definitions> node
            es_customdata = set()
            # cell loop
            for x_celltype in x_root.find('cell_definitions').findall('cell_definition'):

                # <cell_definition>
                s_id = str(x_celltype.get('ID'))
                s_celltype = x_celltype.get('name').replace(' ', '_')
                d_mcds['metadata']['cell_type'].update({s_id : s_celltype})

                # <cell_definition><phenotype>
                x_phenotype = x_celltype.find('phenotype')
                if not (x_phenotype is None):

                    # find <phenotype><cycle> node
                    x_cycle = x_phenotype.find('cycle')
                    if not (x_cycle is None):
                        d_mcds['setting']['parameters'].update({f'{s_celltype}_cycle_model': ds_cycle_model[str(x_cycle.get('code'))]})
                        # <phase_durations>
                        try:
                            for x_phase in x_cycle.find('phase_durations').findall('duration'):
                                s_index = ds_cycle_phase[str(x_phase.get('index'))]
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_cycle_phase_duration_{s_index}_fixed': str(x_phase.get('fixed_duration')).lower() == 'true'})
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_cycle_phase_duration_{s_index}': float(x_phase.text)})
                                d_mcds['setting']['units'].update({f'{s_celltype}_cycle_phase_duration_{s_index}': x_phase.get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><cycle><phase_durations> node missing.')
                        # <phase_transition_rates>
                        try:
                            for x_phase in x_cycle.find('phase_transition_rates').findall('rate'):
                                s_index = ds_cycle_phase[str(x_phase.get('start_index'))] + '_' + ds_cycle_phase[str(x_phase.get('end_index'))]
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_cycle_phase_transition_rate_{s_index}_fixed': str(x_phase.get('fixed_duration')).lower() == 'true'})
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_cycle_phase_transition_rate_{s_index}': float(x_phase.text)})
                                d_mcds['setting']['units'].update({f'{s_celltype}_cycle_phase_transition_rate_{s_index}': x_phase.get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><cycle><phase_transition_rates> node missing.')
                    else:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><cycle> node missing.')

                    # find <phenotype><death> node
                    x_death = x_phenotype.find('death')
                    if not (x_death is None):
                        # <model>
                        for x_model in x_death.findall('model'):
                            s_code = str(x_model.get('code'))
                            s_model = ds_death_model[s_code]  # apoptosis oder necrosis
                            # <death_rate>
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_death_{s_model}_rate': float(x_model.find('death_rate').text)})
                            # <phase_transition_rates>
                            try:
                                for x_phase in x_model.findall('duration'):
                                    s_index = str(x_phase.get('index'))
                                    d_mcds['setting']['parameters'].update({f'{s_celltype}_death_{s_model}_phase_{s_index}_fixed': str(x_phase.get('fixed_duration')).lower() == 'true'})
                                    d_mcds['setting']['parameters'].update({f'{s_celltype}_death_{s_model}_phase_{s_index}': float(x_phase.text)})
                            except AttributeError:
                                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><death><phase_durations> node missing.')
                            try:
                                for x_phase in x_model.findall('rate'):
                                    s_index = str(x_phase.get('start_index')) + '_' + str(x_phase.get('end_index'))
                                    d_mcds['setting']['parameters'].update({f'{s_celltype}_death_{s_model}_phase_{s_index}_fixed': str(x_phase.get('fixed_duration')).lower() == 'true'})
                                    d_mcds['setting']['parameters'].update({f'{s_celltype}_death_{s_model}_phase_{s_index}': float(x_phase.text)})
                            except AttributeError:
                                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><death><phase_transition_rates> node missing.')
                            # <parameters>
                            x_death_model_parameter = x_model.find('parameters')
                            if not (x_death_model_parameter is None):
                                # <unlysed_fluid_change_rate>
                                try:
                                    d_mcds['setting']['parameters'].update({f'{s_celltype}_death_{s_model}_unlysed_fluid_change_rate': float(x_death_model_parameter.find('unlysed_fluid_change_rate').text)})
                                    d_mcds['setting']['units'].update({f'{s_celltype}_death_{s_model}_unlysed_fluid_change_rate': x_death_model_parameter.find('unlysed_fluid_change_rate').get('units')})
                                except AttributeError:
                                    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><death><model code={s_code} name={s_model}><parameters><unlysed_fluid_change_rate> node missing.')
                                # <lysed_fluid_change_rate>
                                try:
                                    d_mcds['setting']['parameters'].update({f'{s_celltype}_death_{s_model}_lysed_fluid_change_rate': float(x_death_model_parameter.find('lysed_fluid_change_rate').text)})
                                    d_mcds['setting']['units'].update({f'{s_celltype}_death_{s_model}_lysed_fluid_change_rate': x_death_model_parameter.find('lysed_fluid_change_rate').get('units')})
                                except AttributeError:
                                    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><death><model code={s_code} name={s_model}><parameters><lysed_fluid_change_rate> node missing.')
                                # <cytoplasmic_biomass_change_rate>
                                try:
                                    d_mcds['setting']['parameters'].update({f'{s_celltype}_death_{s_model}_cytoplasmic_biomass_change_rate': float(x_death_model_parameter.find('cytoplasmic_biomass_change_rate').text)})
                                    d_mcds['setting']['units'].update({f'{s_celltype}_death_{s_model}_cytoplasmic_biomass_change_rate': x_death_model_parameter.find('cytoplasmic_biomass_change_rate').get('units')})
                                except AttributeError:
                                    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><death><model code={s_code} name={s_model}><parameters><cytoplasmic_biomass_change_rate> node missing.')
                                # <nuclear_biomass_change_rate>
                                try:
                                    d_mcds['setting']['parameters'].update({f'{s_celltype}_death_{s_model}_nuclear_biomass_change_rate': float(x_death_model_parameter.find('nuclear_biomass_change_rate').text)})
                                    d_mcds['setting']['units'].update({f'{s_celltype}_death_{s_model}_nuclear_biomass_change_rate': x_death_model_parameter.find('nuclear_biomass_change_rate').get('units')})
                                except AttributeError:
                                    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><death><model code={s_code} name={s_model}><parameters><nuclear_biomass_change_rate> node missing.')
                                # <calcification_rate>
                                try:
                                    d_mcds['setting']['parameters'].update({f'{s_celltype}_death_{s_model}_calcification_rate': float(x_death_model_parameter.find('calcification_rate').text)})
                                    d_mcds['setting']['units'].update({f'{s_celltype}_death_{s_model}_calcification_rate': x_death_model_parameter.find('calcification_rate').get('units')})
                                except AttributeError:
                                    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><death><model code={s_code} name={s_model}><parameters><calcification_rate> node missing.')
                                # <relative_rupture_volume>
                                try:
                                    d_mcds['setting']['parameters'].update({f'{s_celltype}_death_{s_model}_relative_rupture_volume': float(x_death_model_parameter.find('relative_rupture_volume').text)})
                                    d_mcds['setting']['units'].update({f'{s_celltype}_death_{s_model}_relative_rupture_volume': x_death_model_parameter.find('relative_rupture_volume').get('units')})
                                except AttributeError:
                                    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><death><model code={s_code} name={s_model}><parameters><relative_rupture_volume> node missing.')
                            else:
                                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><death><model code={s_code} name={s_model}><parameters> node missing.')
                    else:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><death> node missing.')

                    # find <phenotype><volume> node
                    x_volume = x_phenotype.find('volume')
                    if not (x_volume is None):
                        # <total>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_volume_total': float(x_volume.find('total').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_volume_total': x_volume.find('total').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><volume><total> node missing.')
                        # <fluid_fraction>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_volume_fluid_fraction': float(x_volume.find('fluid_fraction').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_volume_fluid_fraction': x_volume.find('fluid_fraction').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><volume><fluid_fraction> node missing.')
                        # <nuclear>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_volume_nuclear': float(x_volume.find('nuclear').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_volume_nuclear': x_volume.find('nuclear').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><volume><nuclear> node missing.')
                        # <fluid_change_rate>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_volume_fluid_change_rate': float(x_volume.find('fluid_change_rate').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_volume_fluid_change_rate': x_volume.find('fluid_change_rate').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><volume><fluid_change_rate> node missing.')
                        # <cytoplasmic_biomass_change_rate>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_volume_cytoplasmic_biomass_change_rate': float(x_volume.find('cytoplasmic_biomass_change_rate').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_volume_cytoplasmic_biomass_change_rate': x_volume.find('cytoplasmic_biomass_change_rate').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><volume><cytoplasmic_biomass_change_rate> node missing.')
                        # <nuclear_biomass_change_rate>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_volume_nuclear_biomass_change_rate': float(x_volume.find('nuclear_biomass_change_rate').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_volume_nuclear_biomass_change_rate': x_volume.find('nuclear_biomass_change_rate').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><volume><nuclear_biomass_change_rate> node missing.')
                        # <calcified_fraction>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_volume_calcified_fraction': float(x_volume.find('calcified_fraction').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_volume_calcified_fraction': x_volume.find('calcified_fraction').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><volume><calcified_fraction> node missing.')
                        # <calcification_rate>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_volume_calcification_rate': float(x_volume.find('calcification_rate').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_volume_calcification_rate': x_volume.find('calcification_rate').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><volume><calcification_rate> node missing.')
                        # <relative_rupture_volume>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_volume_relative_rupture_volume': float(x_volume.find('relative_rupture_volume').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_volume_relative_rupture_volume': x_volume.find('relative_rupture_volume').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><volume><relative_rupture_volume> node missing.')
                    else:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><volume> node missing.')

                    # find <phenotype><mechanics> node
                    x_mechanics = x_phenotype.find('mechanics')
                    if not (x_mechanics is None):
                        # <cell_cell_adhesion_strength>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_mechanics_cell_cell_adhesion_strength': float(x_mechanics.find('cell_cell_adhesion_strength').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_mechanics_cell_cell_adhesion_strength': x_mechanics.find('cell_cell_adhesion_strength').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><mechanics><cell_cell_adhesion_strength> node missing.')
                        # <cell_cell_repulsion_strength>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_mechanics_cell_cell_repulsion_strength': float(x_mechanics.find('cell_cell_repulsion_strength').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_mechanics_cell_cell_repulsion_strength': x_mechanics.find('cell_cell_repulsion_strength').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><mechanics><cell_cell_repulsion_strength> node missing.')
                        # cell_BM_adhesion_strength
                        #try:
                        #    d_mcds['setting']['parameters'].update({f'{s_celltype}_mechanics_cell_BM_adhesion_strength': float(x_mechanics.find('cell_BM_adhesion_strength').text)})
                        #    d_mcds['setting']['units'].update({f'{s_celltype}_mechanics_cell_BM_adhesion_strength': x_mechanics.find('cell_BM_adhesion_strength').get('units')})
                        #except AttributeError:
                        #    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><mechanics><cell_BM_adhesion_strength> node missing.')
                        # cell_BM_repulsion_strength
                        #try:
                        #    d_mcds['setting']['parameters'].update({f'{s_celltype}_mechanics_cell_BM_repulsion_strength': float(x_mechanics.find('cell_BM_repulsion_strength').text)})
                        #    d_mcds['setting']['units'].update({f'{s_celltype}_mechanics_cell_BM_repulsion_strength': x_mechanics.find('cell_BM_repulsion_strength').get('units')})
                        #except AttributeError:
                        #    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><mechanics><cell_BM_repulsion_strength> node missing.')
                        # <relative_maximum_adhesion_distance>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_mechanics_relative_maximum_adhesion_distance': float(x_mechanics.find('relative_maximum_adhesion_distance').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_mechanics_relative_maximum_adhesion_distance': x_mechanics.find('relative_maximum_adhesion_distance').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><mechanics><relative_maximum_adhesion_distance> node missing.')

                        # <cell_adhesion_affinities>
                        x_mechanics_affinities = x_mechanics.find('cell_adhesion_affinities')
                        if not (x_mechanics_affinities is None):
                            for x_affinity in x_mechanics_affinities.findall('cell_adhesion_affinity'):
                                s_name = x_affinity.get('name')
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_mechanics_relative_cell_adhesion_affinity_{s_name}': float(x_affinity.text)})
                        else:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><mechanics><cell_adhesion_affinities> node missing.')
                        # <options>
                        x_mechanics_options = x_mechanics.find('options')
                        if not (x_mechanics_options is None):
                            # <options><set_relative_equilibrium_distance>
                            try:
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_mechanics_relative_equilibrium_distance_enabled': str(x_mechanics_options.find('set_relative_equilibrium_distance').get('enabled').lower() == 'true')})
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_mechanics_relative_equilibrium_distance': float(x_mechanics_options.find('set_relative_equilibrium_distance').text)})
                                d_mcds['setting']['units'].update({f'{s_celltype}_mechanics_relative_equilibrium_distance': x_mechanics_options.find('set_relative_equilibrium_distance').get('units')})
                            except AttributeError:
                                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><mechanics><options><set_relative_equilibrium_distance> node missing.')
                            # <options><set_absolute_equilibrium_distance>
                            try:
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_mechanics_absolute_equilibrium_distance_enabled': str(x_mechanics_options.find('set_absolute_equilibrium_distance').get('enabled').lower() == 'true')})
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_mechanics_absolute_equilibrium_distance': float(x_mechanics_options.find('set_absolute_equilibrium_distance').text)})
                                d_mcds['setting']['units'].update({f'{s_celltype}_mechanics_absolute_equilibrium_distance': x_mechanics_options.find('set_absolute_equilibrium_distance').get('units')})
                            except AttributeError:
                                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><mechanics><options><set_absolute_equilibrium_distance> node missing.')
                        else:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><mechanics><options> node missing.')

                        # <attachment_elastic_constant>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_mechanics_elastic_constant': float(x_mechanics.find('attachment_elastic_constant').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_mechanics_elastic_constant': x_mechanics.find('attachment_elastic_constant').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><mechanics><attachment_elastic_constant> node missing.')
                        # <attachment_rate>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_mechanics_attachment_rate': float(x_mechanics.find('attachment_rate').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_mechanics_attachment_rate': x_mechanics.find('attachment_rate').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><mechanics><attachment_rate> node missing.')
                        # <detachment_rate>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_mechanics_detachment_rate': float(x_mechanics.find('detachment_rate').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_mechanics_detachment_rate': x_mechanics.find('detachment_rate').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><mechanics><detachment_rate> node missing.')

                    else:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><mechanics> node missing.')

                    # find <phenotype><motility> node
                    x_motility = x_phenotype.find('motility')
                    if not (x_mechanics is None):
                        # <speed>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_motility_speed': float(x_motility.find('speed').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_motility_speed': x_motility.find('speed').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><motility><speed> node missing.')
                        # <persistence_time>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_motility_persistence_time': float(x_motility.find('persistence_time').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_motility_persistence_time': x_motility.find('persistence_time').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><motility><persistence_time> node missing.')
                        # <migration_bias>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_motility_migration_bias': float(x_motility.find('migration_bias').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_motility_migration_bias': x_motility.find('migration_bias').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><motility><_motility_migration_bias> node missing.')
                        # <options>
                        x_motility_options = x_motility.find('options')
                        if not (x_motility_options is None):
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_motility_enabled': str(x_motility_options.find('enabled').text).lower() == 'true'})
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_motility_use_2D': str(x_motility_options.find('use_2D').text).lower() == 'true'})
                            # <options><chemotaxis>
                            x_motility_chemotaxis = x_motility_options.find('chemotaxis')
                            if not (x_motility_chemotaxis is None):
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_motility_chemotaxis_enabled': str(x_motility_chemotaxis.find('enabled').text).lower() == 'true'})
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_motility_chemotaxis_substrate': str(x_motility_chemotaxis.find('substrate').text)})
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_motility_chemotaxis_direction': float(x_motility_chemotaxis.find('direction').text)})
                            else:
                                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><motility><options><chemotaxis> node missing.')
                            # <options><advanced_chemotaxis>
                            x_motility_advancedchemotaxis = x_motility_options.find('advanced_chemotaxis')
                            if not (x_motility_advancedchemotaxis is None):
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_motility_advanced_chemotaxis_enabled': str(x_motility_advancedchemotaxis.find('enabled').text).lower() == 'true'})
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_motility_advanced_chemotaxis_normalize_each_gradient': str(x_motility_advancedchemotaxis.find('normalize_each_gradient').text).lower() == 'true'})
                                for x_chemotactic in x_motility_advancedchemotaxis.find('chemotactic_sensitivities').findall('chemotactic_sensitivity'):
                                    s_subs = str(x_chemotactic.get('substrate'))
                                    d_mcds['setting']['parameters'].update({f'{s_celltype}_motility_advanced_chemotaxis_chemotactic_sensitivity_{s_subs}': float(x_chemotactic.text)})
                            else:
                                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><motility><options><advanced_chemotaxis> node missing.')
                        else:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><motility><options> node missing.')
                    else:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><motility> node missing.')

                    # find <phenotype><secretion> node
                    x_secretion = x_phenotype.find('secretion')
                    if not (x_secretion is None):
                        for x_substrate in x_secretion.findall('substrate'):
                            s_subs = str(x_substrate.get('name'))
                            # <secretion_rate>
                            try:
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_{s_subs}_secretion_rate': float(x_substrate.find('secretion_rate').text)})
                                d_mcds['setting']['units'].update({f'{s_celltype}_{s_subs}_secretion_rate': x_substrate.find('secretion_rate').get('units')})
                            except AttributeError:
                                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><secretion><substrate><secretion_rate> node missing.')
                            # <secretion_target>
                            try:
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_{s_subs}_secretion_target': float(x_substrate.find('secretion_target').text)})
                                d_mcds['setting']['units'].update({f'{s_celltype}_{s_subs}_secretion_target': x_substrate.find('secretion_target').get('units')})
                            except AttributeError:
                                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><secretion><substrate><secretion_target> node missing.')
                            # <uptake_rate>
                            try:
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_{s_subs}_uptake_rate': float(x_substrate.find('uptake_rate').text)})
                                d_mcds['setting']['units'].update({f'{s_celltype}_{s_subs}_uptake_rate': x_substrate.find('uptake_rate').get('units')})
                            except AttributeError:
                                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><secretion><substrate><uptake_rate> node missing.')
                            # <net_export_rate>
                            try:
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_{s_subs}_net_export_rate': float(x_substrate.find('net_export_rate').text)})
                                d_mcds['setting']['units'].update({f'{s_celltype}_{s_subs}_net_export_rate': x_substrate.find('net_export_rate').get('units')})
                            except AttributeError:
                                print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><secretion><substrate><net_export_rate> node missing.')
                    else:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><secretion> node missing.')

                    # find <phenotype><cell_interactions> node
                    x_interaction = x_phenotype.find('cell_interactions')
                    if not (x_interaction is None):
                        # <dead_phagocytosis_rate>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_cell_interaction_dead_phagocytosis_rate': float(x_interaction.find('dead_phagocytosis_rate').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_cell_interaction_dead_phagocytosis_rate': x_interaction.find('dead_phagocytosis_rate').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><cell_interactions><dead_phagocytosis_rate> node missing.')
                        # <live_phagocytosis_rates>
                        x_interaction_phagocytosis = x_interaction.find('live_phagocytosis_rates')
                        if not (x_interaction_phagocytosis is None):
                            for x_rate in x_interaction_phagocytosis:
                                s_name = x_rate.get('name')
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_cell_interaction_{s_name}_live_phagocytosis_rate': float(x_rate.text)})
                                d_mcds['setting']['units'].update({f'{s_celltype}_cell_interaction_{s_name}_live_phagocytosis_rate': x_rate.get('units')})
                        else:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><cell_interactions><live_phagocytosis_rates> node missing.')
                        # <attack_rates>
                        x_interaction_attack = x_interaction.find('attack_rates')
                        if not (x_interaction_attack is None):
                            for x_rate in x_interaction_attack:
                                s_name = x_rate.get('name')
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_cell_interaction_{s_name}_attack_rate': float(x_rate.text)})
                                d_mcds['setting']['units'].update({f'{s_celltype}_cell_interaction_{s_name}_attack_rate': x_rate.get('units')})
                        else:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><cell_interactions><attack_rates> node missing.')

                        # <damage_rate>
                        try:
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_cell_interaction_damage_rate': float(x_interaction.find('damage_rate').text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_cell_interaction_damage_rate': x_interaction.find('damage_rate').get('units')})
                        except AttributeError:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><cell_interactions><damage_rate> node missing.')
                        # <fusion_rates>
                        x_interaction_fusion =  x_interaction.find('fusion_rates')
                        if not (x_interaction_fusion is None):
                            for x_rate in x_interaction_fusion:
                                s_name = x_rate.get('name')
                                d_mcds['setting']['parameters'].update({f'{s_celltype}_cell_interaction_{s_name}_fusion_rate': float(x_rate.text)})
                                d_mcds['setting']['units'].update({f'{s_celltype}_cell_interaction_{s_name}_fusion_rate': x_rate.get('units')})
                        else:
                            print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><cell_interactions><fusion_rates> node missing.')

                    else:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><cell_interactions> node missing.')

                    # find <phenotype><cell_transformations> node
                    x_transformation = x_phenotype.find('cell_transformations')
                    if not (x_transformation is None):
                        # <transformation_rates>
                        for x_rate  in x_transformation.find('transformation_rates').findall('transformation_rate'):
                            s_name = str(x_rate.get('name'))
                            d_mcds['setting']['parameters'].update({f'{s_celltype}_{s_name}_transformation_rate': float(x_rate.text)})
                            d_mcds['setting']['units'].update({f'{s_celltype}_{s_name}_transformation_rate': x_rate.get('units')})
                    else:
                        print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype><cell_transformations> node missing.')

                    # skip <phenotype><molecular> (extinct version)
                    # skip <phenotype><intracellular>

                # <cell_definition><phenotype>
                else:
                    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><phenotype> node missing.')

                # <cell_definition><custom_data>
                try:
                    for x_element in x_celltype.find('custom_data').iter():
                        if (x_element.tag != 'custom_data'):
                            try:
                                self.custom_type[x_element.tag]
                            except KeyError:
                                es_customdata.add(x_element.tag)
                except AttributeError:
                    print(f'Warning @ pyMCDS._read_setting_xml : <cell_definition name="{s_celltype}" ID="{s_id}"><custom_data> node missing.')

            # if <custom data> was found
            if (len(es_customdata) > 0):
                print(f'Warning @ pyMCDS._read_setting_xml : cell_definition custom_data without variable type setting detected. {sorted(es_customdata)}')

            # find <user_parameters> node
            for x_element in x_root.find('user_parameters').iter():
                if (x_element.tag != 'user_parameters'):
                    s_type = x_element.get('type')
                    if s_type in {'double', 'real', 'float'}:
                        o_type = float
                    elif s_type in {'int', 'integer'}:
                        o_type = int
                    elif s_type in {'str', 'string'}:
                        o_type = str
                    elif s_type in {'bool','boole'}:
                        print(f'Warning @  pyMCDS._read_setting.xml: user_parameter bool variable type detected. value {x_element.tag} might be ambigiuos. translate into string.')
                        o_type = str
                    else:
                        print(f'Warning @  pyMCDS._read_setting.xml: user_parameter unknown variable type {s_type} detected. will translate into string.')
                    o_type = str
                    d_mcds['setting']['parameters'].update({x_element.tag: o_type(x_element.text)})

            # skip <initial_conditions> node

            # find <cell_rules> node (protocol CBHG; version 2.0)
            x_rule = x_root.find('cell_rules')
            if not (x_rule is None):
                for x_ruleset in x_rule.find('rulesets'):
                    if (str(x_ruleset.get('enabled')).lower() == 'true'):
                        s_pathfile = self.path + '/' + x_ruleset.find('filename').text
                        try:
                            df_rule = pd.read_csv(s_pathfile, sep=',', header=None)
                            if (df_rule.shape[1] == 9):
                                df_rule.columns = ['cell_type','signal','direction','behavoir','base_value', 'saturation_value','half_max','hill_power','apply_to_dead']
                                df_rule.drop({'base_value'}, axis=1, inplace=True)
                            else:
                                df_rule.columns = ['cell_type','signal','direction','behavoir','saturation_value','half_max','hill_power','apply_to_dead']
                            df_rule = df_rule.astype({'apply_to_dead':bool})
                            df_rule.index.name = 'index'
                            if (d_mcds['setting']['rules'] is None):
                                d_mcds['setting']['rules'] = df_rule
                            else:
                                d_mcds['setting']['rules'] = pd.concate([d_mcds['setting']['rules'], df_rule], axis=0)
                        except pd._libs.parsers.EmptyDataError:
                            print(f'Warning @ pyMCDS._read_settings.xml : {s_pathfile} is empty.')
                if (d_mcds['setting']['rules'] is None):
                    print(f'Warning @ pyMCDS._read_setting_xml : no enabled cell_rules.csv detetcetd.')
            else:
                print(f'Warning @ pyMCDS._read_setting_xml : <cell_rules> node missing.')


        #######################################
        # read physicell output xml path/file #
        #######################################

        s_xmlpathfile = self.path + '/' + self.xmlfile
        x_tree = ET.parse(s_xmlpathfile)
        if self.verbose:
            print(f'reading: {s_xmlpathfile}')
        x_root = x_tree.getroot()


        ###################
        # handle metadata #
        ###################

        if self.verbose:
            print('working on metadata ...')

        ### find the metadata node ###
        x_metadata = x_root.find('metadata')

        # get multicellds xml version
        d_mcds['metadata']['multicellds_version'] = f"MultiCellDS_{x_root.get('version')}"

        # get physicell software version
        x_software = x_metadata.find('software')
        x_physicelln = x_software.find('name')
        x_physicellv = x_software.find('version')
        d_mcds['metadata']['physicell_version'] = f'{x_physicelln.text}_{x_physicellv.text}'

        # get timestamp
        x_time = x_metadata.find('created')
        d_mcds['metadata']['created'] = x_time.text

        # get current simulated time
        x_time = x_metadata.find('current_time')
        d_mcds['metadata']['current_time'] = float(x_time.text)
        d_mcds['metadata']['time_units'] = x_time.get('units')

        # get current runtime
        x_time = x_metadata.find('current_runtime')
        d_mcds['metadata']['current_runtime'] = float(x_time.text)
        d_mcds['metadata']['runtime_units'] = x_time.get('units')

        # update settings unit with metadata infromation
        d_mcds['setting']['units'].update({'time': d_mcds['metadata']['time_units']})
        d_mcds['setting']['units'].update({'runtime': d_mcds['metadata']['runtime_units']})


        ####################
        # handle mesh data #
        ####################

        if self.verbose:
            print('working on mesh data ...')

        ### find the mesh node ###
        x_microenv = x_root.find('microenvironment').find('domain')  # find the microenvironment node
        x_mesh = x_microenv.find('mesh')
        d_mcds['metadata']['spatial_units'] = x_mesh.get('units')

        # while we're at it, find the mesh
        s_x_coor = x_mesh.find('x_coordinates').text
        s_delim = x_mesh.find('x_coordinates').get('delimiter')
        ar_x_coor = np.array(s_x_coor.split(s_delim), dtype=np.float64)

        s_y_coor = x_mesh.find('y_coordinates').text
        s_delim = x_mesh.find('y_coordinates').get('delimiter')
        ar_y_coor = np.array(s_y_coor.split(s_delim), dtype=np.float64)

        s_z_coor = x_mesh.find('z_coordinates').text
        s_delim = x_mesh.find('z_coordinates').get('delimiter')
        ar_z_coor = np.array(s_z_coor.split(s_delim), dtype=np.float64)

        # reshape into a meshgrid
        d_mcds['mesh']['mnp_grid'] = np.array(np.meshgrid(ar_x_coor, ar_y_coor, ar_z_coor, indexing='xy'))

        # get mesh center axis
        d_mcds['mesh']['mnp_axis'] = [
            np.unique(ar_x_coor),
            np.unique(ar_y_coor),
            np.unique(ar_z_coor),
        ]

        # get mesh center range
        d_mcds['mesh']['mnp_range'] = [
           (d_mcds['mesh']['mnp_axis'][0].min(), d_mcds['mesh']['mnp_axis'][0].max()),
           (d_mcds['mesh']['mnp_axis'][1].min(), d_mcds['mesh']['mnp_axis'][1].max()),
           (d_mcds['mesh']['mnp_axis'][2].min(), d_mcds['mesh']['mnp_axis'][2].max()),
        ]

        # get voxel range
        d_mcds['mesh']['ijk_range'] = [
            (0, len(d_mcds['mesh']['mnp_axis'][0]) - 1),
            (0, len(d_mcds['mesh']['mnp_axis'][1]) - 1),
            (0, len(d_mcds['mesh']['mnp_axis'][2]) - 1),
        ]

        # get voxel axis
        d_mcds['mesh']['ijk_axis'] = [
            np.array(range(d_mcds['mesh']['ijk_range'][0][1] + 1)),
            np.array(range(d_mcds['mesh']['ijk_range'][1][1] + 1)),
            np.array(range(d_mcds['mesh']['ijk_range'][2][1] + 1)),
        ]

        # get mesh bounding box range [xmin, ymin, zmin, xmax, ymax, zmax]
        s_bboxcoor = x_mesh.find('bounding_box').text
        s_delim = x_mesh.find('bounding_box').get('delimiter')
        ar_bboxcoor = np.array(s_bboxcoor.split(s_delim), dtype=np.float64)

        d_mcds['mesh']['xyz_range'] = [
            (ar_bboxcoor[0], ar_bboxcoor[3]),
            (ar_bboxcoor[1], ar_bboxcoor[4]),
            (ar_bboxcoor[2], ar_bboxcoor[5]),
        ]

        # voxel data must be loaded from .mat file
        s_voxelpathfile = self.path + '/' + x_mesh.find('voxels').find('filename').text
        ar_mesh_initial = io.loadmat(s_voxelpathfile)['mesh']
        if self.verbose:
            print(f'reading: {s_voxelpathfile}')

        # center of voxel specified by first three rows [ x, y, z ]
        # volume specified by fourth row
        d_mcds['mesh']['mnp_coordinate'] = ar_mesh_initial[:3, :]
        d_mcds['mesh']['volumes'] = ar_mesh_initial[3, :]

        # update settings unit with mesh infromation
        d_mcds['setting']['units'].update({'spatial_unit': d_mcds['metadata']['spatial_units']})

        ################################
        # handle microenvironment data #
        ################################

        if self.microenv:
            if self.verbose:
                print('working on microenvironment data ...')

            # micro environment data is shape [4+n, len(voxels)] where n is the number
            # of species being tracked. the first 3 rows represent (x, y, z) of voxel
            # centers. The fourth row contains the voxel volume. The 5th row and up will
            # contain values for that species in that voxel.
            s_microenvpathfile = self.path + '/' +  x_microenv.find('data').find('filename').text
            ar_microenv = io.loadmat(s_microenvpathfile)['multiscale_microenvironment']
            if self.verbose:
                print(f'reading: {s_microenvpathfile}')

            # continuum_variables, unlike in the matlab version the individual chemical
            # species will be primarily accessed through their names e.g.
            # d_mcds['continuum_variables']['oxygen']['units']
            # d_mcds['continuum_variables']['glucose']['data']

            # substrate loop
            for i_s, x_substrate in enumerate(x_microenv.find('variables').findall('variable')):
                # i don't like spaces in species names!
                s_substrate = x_substrate.get('name').replace(' ', '_')

                d_mcds['continuum_variables'][s_substrate] = {}
                d_mcds['continuum_variables'][s_substrate]['units'] = x_substrate.get('units')

                if self.verbose:
                    print(f'parsing: {s_substrate} data')

                # update metadata substrate ID label dictionary
                d_mcds['metadata']['substrate'].update({str(i_s) : s_substrate})

                # initialize meshgrid shaped array for concentration data
                d_mcds['continuum_variables'][s_substrate]['data'] = np.zeros(d_mcds['mesh']['mnp_grid'][0].shape)

                # diffusion data for each species
                d_mcds['continuum_variables'][s_substrate]['diffusion_coefficient'] = {}
                d_mcds['continuum_variables'][s_substrate]['diffusion_coefficient']['value'] = float(x_substrate.find('physical_parameter_set').find('diffusion_coefficient').text)
                d_mcds['continuum_variables'][s_substrate]['diffusion_coefficient']['units'] = x_substrate.find('physical_parameter_set').find('diffusion_coefficient').get('units')

                # decay data for each species
                d_mcds['continuum_variables'][s_substrate]['decay_rate'] = {}
                d_mcds['continuum_variables'][s_substrate]['decay_rate']['value']  = float(x_substrate.find('physical_parameter_set').find('decay_rate').text)
                d_mcds['continuum_variables'][s_substrate]['decay_rate']['units']  = x_substrate.find('physical_parameter_set').find('decay_rate').get('units')

                # store data from microenvironment file as numpy array
                # iterate over each voxel
                # bue: i have a hunch this could be faster reimplemented.
                for vox_idx in range(d_mcds['mesh']['mnp_coordinate'].shape[1]):

                    # find the voxel coordinate
                    ar_center = d_mcds['mesh']['mnp_coordinate'][:, vox_idx]
                    i = np.where(np.abs(ar_center[0] - d_mcds['mesh']['mnp_axis'][0]) < 1e-10)[0][0]
                    j = np.where(np.abs(ar_center[1] - d_mcds['mesh']['mnp_axis'][1]) < 1e-10)[0][0]
                    k = np.where(np.abs(ar_center[2] - d_mcds['mesh']['mnp_axis'][2]) < 1e-10)[0][0]

                    # store value
                    d_mcds['continuum_variables'][s_substrate]['data'][j, i, k] = ar_microenv[4+i_s, vox_idx]

                # update settings unit wuth substrate parameters
                d_mcds['setting']['units'].update({s_substrate: d_mcds['continuum_variables'][s_substrate]['units']})
                # update settings unit with microenvironment parameters
                d_mcds['setting']['units'].update({f'{s_substrate}_diffusion_coefficient': d_mcds['continuum_variables'][s_substrate]['diffusion_coefficient']['units']})
                d_mcds['setting']['units'].update({f'{s_substrate}_decay_rate': d_mcds['continuum_variables'][s_substrate]['decay_rate']['units']})


        ####################
        # handle cell data #
        ####################

        if self.verbose:
            print('working on discrete cell data ...')

        # in order to get to the good stuff, we have to pass through a few different hierarchical levels
        x_cell = x_root.find('cellular_information').find('cell_populations').find('cell_population').find('custom')

        # we want the PhysiCell data, there is more of it
        for x_simplified_data in x_cell.findall('simplified_data'):
            if x_simplified_data.get('source') == 'PhysiCell':
                x_celldata = x_simplified_data
                break

        # iterate over labels which are children of labels these will be used to label data arrays
        ls_variable = []
        for label in x_celldata.find('labels').findall('label'):
            # I don't like spaces in my dictionary keys!
            s_variable = label.text.replace(' ', '_')
            i_variable = int(label.get('size'))
            s_unit = label.get('units')

            if s_variable in es_var_subs:
                if (len(d_mcds['metadata']['substrate']) > 0):
                    # continuum_variable id label sorting
                    ls_substrate = [s_substrate for _, s_substrate in sorted(d_mcds['metadata']['substrate'].items())]
                    for s_substrate in ls_substrate:
                        s_variable_subs = s_substrate + '_' + s_variable
                        ls_variable.append(s_variable_subs)
                        d_mcds['discrete_cells']['units'].update({s_variable_subs : s_unit})
                else:
                    ls_substrate = [str(i_substrate) for i_substrate in range(i_variable)]
                    for s_substrate in ls_substrate:
                        s_variable_subs = s_variable + '_' + s_substrate
                        ls_variable.append(s_variable_subs)
                        d_mcds['discrete_cells']['units'].update({s_variable_subs : s_unit})

            elif s_variable in es_var_death:
                for i_deathrate in range(i_variable):
                    s_variable_deathrate = s_variable + '_' + str(i_deathrate)
                    ls_variable.append(s_variable_deathrate)
                    d_mcds['discrete_cells']['units'].update({s_variable_deathrate : s_unit})

            elif s_variable in es_var_cell:
                if (len(d_mcds['metadata']['cell_type']) > 0):
                    # discrete_cells id label sorting
                    ls_celltype = [s_celltype for _, s_celltype in sorted(d_mcds['metadata']['cell_type'].items())]
                    for s_celltype in ls_celltype:
                        s_variable_celltype = s_celltype + '_' + s_variable
                        ls_variable.append(s_variable_celltype)
                        d_mcds['discrete_cells']['units'].update({s_variable_celltype : s_unit})
                else:
                    ls_celltype = [str(i_celltype) for i_celltype in range(i_variable)]
                    for s_celltype in ls_celltype:
                        s_variable_celltype = s_variable + '_' + s_celltype
                        ls_variable.append(s_variable_celltype)
                        d_mcds['discrete_cells']['units'].update({s_variable_celltype : s_unit})

            elif s_variable in es_var_spatial:
                for s_axis in ['_x','_y','_z']:
                    s_variable_spatial = s_variable + s_axis
                    ls_variable.append(s_variable_spatial)
                    d_mcds['discrete_cells']['units'].update({s_variable_spatial: s_unit})

            else:
                ls_variable.append(s_variable)
                d_mcds['discrete_cells']['units'].update({s_variable : s_unit})

        # update settings units from cell parameters
        d_mcds['setting']['units'].update(d_mcds['discrete_cells']['units'])
        del d_mcds['setting']['units']['ID']

        # load the file
        s_cellpathfile = self.path + '/' + x_celldata.find('filename').text
        try:
            ar_cell = io.loadmat(s_cellpathfile)['cells']
            if self.verbose:
                print(f'reading: {s_cellpathfile}')
        except ValueError:  # hack: some old PhysiCell versions generates a corrupt cells.mat file, if there are zero cells.
            print(f'Warning @ pyMCDS._read_xml : corrupt {s_cellpathfile} detected!\nassuming time step with zero cells because of a known bug in PhysiCell MultiCellDS version 0.5 output.')
            ar_cell = np.empty([len(ls_variable),0])

        # store data
        d_mcds['discrete_cells']['data'] = {}
        for col in range(len(ls_variable)):
            d_mcds['discrete_cells']['data'][ls_variable[col]] = ar_cell[col,:]


        #####################
        # handle graph data #
        #####################

        d_mcds['discrete_cells']['graph'] = {}
        d_mcds['discrete_cells']['graph'].update({'neighbor_cells': {}})
        d_mcds['discrete_cells']['graph'].update({'attached_cells': {}})

        if self.graph:
            if self.verbose:
                print('working on graph data ...')

            # neighborhood cell graph
            s_cellpathfile = self.path + '/' + x_cell.find('neighbor_graph').find('filename').text
            dei_graph = graphfile_parser(s_pathfile=s_cellpathfile)
            if self.verbose:
                print(f'reading: {s_cellpathfile}')

            # store data
            d_mcds['discrete_cells']['graph'].update({'neighbor_cells': dei_graph})

            # attached cell graph
            s_cellpathfile = self.path + '/' + x_cell.find('attached_cells_graph').find('filename').text
            dei_graph = graphfile_parser(s_pathfile=s_cellpathfile)
            if self.verbose:
                print(f'reading: {s_cellpathfile}')

            # store data
            d_mcds['discrete_cells']['graph'].update({'attached_cells': dei_graph})


        #########################
        # handle physiboss data #
        #########################

        d_mcds['discrete_cells']['physiboss'] = None

        if self.physiboss:
            if self.verbose:
                print('working on physiboss data ...')

            # intracellular file (hack because this is not yet in output.xml)
            df_physiboss = None
            s_intracellpathfile = self.path + '/' + f'states_{self.xmlfile.replace("output","").replace(".xml",".csv")}'
            if os.path.exists(s_intracellpathfile):
                if self.verbose:
                    print(f'reading: {s_intracellpathfile}')

                # load state
                df_physiboss = pd.read_csv(s_intracellpathfile, index_col=0)

                # add nodes
                df_physiboss[f'node_<nil>'] = df_physiboss.state.str.find('<nil>') > -1
                es_node = set()
                for s_state in df_physiboss.state.unique():
                    es_node = es_node.union(set(s_state.split(' -- ')))
                es_node.discard('<nil>')
                for s_node in sorted(es_node):
                    df_physiboss[f'node_{s_node}'] = df_physiboss.state.str.find(s_node) > -1

            else:
                print(f'Warning @ pyMCDS._read_xml : physiboss file missing {s_intracellpathfile}.')

            # store data
            d_mcds['discrete_cells']['physiboss'] = df_physiboss


        # output
        if self.verbose:
            print('done!')
        return d_mcds
