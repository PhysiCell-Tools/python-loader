#####
# title: test_snapshot_3d.py
#
# language: python3
# author: bue
# date: 2022-10-15
# license: BSD 3-Clause
#
# description:
#   pytest unit test library for heat.py
#   + https://docs.pytest.org/
#
#   note:
#   assert actual == expected, message
#   == value equality
#   is reference equality
#   pytest.approx for real values
#####

# load library
import pathlib
import pcDataLoader as pc

# const
s_path_3d = str(pathlib.Path(pc.__file__).parent.resolve()/'data_timeseries_3d')
s_file_3d = 'output00000024.xml'
s_pathfile_3d = f'{s_path_3d}/{s_file_3d}'


# load physicell data shortcut
class TestPyMcdsShortcut3D(object):
    ''' test for pc.pyMCDS data loader shortcut '''

    def test_pyMCDS(self):
        # load physicell data shortcut
        print(f'process: pc.pyMCDS(xml_file={s_pathfile_3d}) ...')
        mcds = pc.pyMCDS(xml_file=s_pathfile_3d)
        assert str(type(mcds)) == "<class 'pcDataLoader.pyMCDS.pyMCDS'>"


# load physicell data with microenvironment
class TestPyMcds3D(object):
    ''' test for pc.pyMCDS data loader microenvironment True'''
    mcds = pc.pyMCDS(xml_file=s_file_3d, output_path=s_path_3d, microenv=True)

    def test_pyMCDS(self, mcds=mcds):
        # load physicell data
        print(f'process: pc.pyMCDS(xml_file={s_file_3d}, output_path={s_path_3d}, microenv=True) ...')
        assert str(type(mcds)) == "<class 'pcDataLoader.pyMCDS.pyMCDS'>"

    ## metadata realted functions
    def test_mcds_get_physicel_version(self, mcds=mcds):
        s_pcversion = mcds.get_physicell_version()
        assert s_pcversion == '1.10.4'

    def test_mcds_get_timestamp(self, mcds=mcds):
        s_timestamp = mcds.get_timestamp()
        assert s_timestamp == '2022-10-18T15:28:37Z'

    def test_mcds_get_time(self, mcds=mcds):
        r_time = mcds.get_time()
        assert r_time == 1440.0

    def test_mcds_get_runtime(self, mcds=mcds):
        r_runtime = mcds.get_runtime()
        assert r_runtime == 313.35264

    ## graph related functions
    def test_mcds_get_attached_graph_dict(self, mcds=mcds):
        dei_graph = mcds.data['discrete_cells']['graph']['attached_cells']
        assert (str(type(dei_graph)) == "<class 'dict'>") and \
               (len(dei_graph) == 20460) and \
               (len(dei_graph[20459]) == 0)

    def test_mcds_get_neighbor_graph_dict(self, mcds=mcds):
        dei_graph = mcds.data['discrete_cells']['graph']['neighbor_cells']
        assert (str(type(dei_graph)) == "<class 'dict'>") and \
               (len(dei_graph) == 20460) and \
               (len(dei_graph[20459]) == 13)

    ## mesh related functions
    def test_mcds_get_x_range(self, mcds=mcds):
        tr_range = mcds.get_x_range()
        assert tr_range == (-30, 300)

    def test_mcds_get_y_range(self, mcds=mcds):
        tr_range = mcds.get_y_range()
        assert tr_range == (-20, 200)

    def test_mcds_get_z_range(self, mcds=mcds):
        tr_range = mcds.get_z_range()
        assert tr_range == (-10, 100)

    def test_mcds_get_mesh_m_range(self, mcds=mcds):
        tr_range = mcds.get_mesh_m_range()
        assert tr_range == (-15, 285)

    def test_mcds_get_mesh_n_range(self, mcds=mcds):
        tr_range = mcds.get_mesh_n_range()
        assert tr_range == (-10, 190)

    def test_mcds_get_mesh_p_range(self, mcds=mcds):
        tr_range = mcds.get_mesh_p_range()
        assert tr_range == (-5, 95)

    def test_mcds_get_voxel_i_range(self, mcds=mcds):
        tr_range = mcds.get_voxel_i_range()
        assert tr_range == (0, 11)

    def test_mcds_get_voxel_j_range(self, mcds=mcds):
        tr_range = mcds.get_voxel_j_range()
        assert tr_range == (0, 11)

    def test_mcds_get_voxel_k_range(self, mcds=mcds):
        tr_range = mcds.get_voxel_k_range()
        assert tr_range == (0, 11)

    def test_mcds_get_mesh_flat_false(self, mcds=mcds):
        lar_mesh = mcds.get_mesh(flat=False)
        assert (str(type(lar_mesh)) == "<class 'list'>") and \
               (len(lar_mesh) == 3) and \
               (str(type(lar_mesh[0])) == "<class 'numpy.ndarray'>") and \
               (str(lar_mesh[0].dtype) == "float64") and \
               (lar_mesh[0].shape == (11, 11, 11)) and \
               (str(type(lar_mesh[1])) == "<class 'numpy.ndarray'>") and \
               (str(lar_mesh[1].dtype) == "float64") and \
               (lar_mesh[1].shape == (11, 11, 11)) and \
               (str(type(lar_mesh[2])) == "<class 'numpy.ndarray'>") and \
               (str(lar_mesh[2].dtype) == "float64") and \
               (lar_mesh[2].shape == (11, 11, 11))

    def test_mcds_get_mesh_flat_true(self, mcds=mcds):
        lar_mesh = mcds.get_mesh(flat=True)
        assert (str(type(lar_mesh)) == "<class 'list'>") and \
               (len(lar_mesh) == 2) and \
               (str(type(lar_mesh[0])) == "<class 'numpy.ndarray'>") and \
               (str(lar_mesh[0].dtype) == "float64") and \
               (lar_mesh[0].shape == (11, 11)) and \
               (str(type(lar_mesh[1])) == "<class 'numpy.ndarray'>") and \
               (str(lar_mesh[1].dtype) == "float64") and \
               (lar_mesh[1].shape == (11, 11))

    def test_mcds_get_mesh_2d(self, mcds=mcds):
        lar_mesh_flat = mcds.get_mesh(flat=True)
        lar_mesh_2d = mcds.get_mesh_2D()
        assert (str(type(lar_mesh_2d)) == "<class 'list'>") and \
               (len(lar_mesh_2d) == 2) and \
               (lar_mesh_2d[0] == lar_mesh_flat[0]).all() and \
               (lar_mesh_2d[1] == lar_mesh_flat[1]).all()

    def test_mcds_get_mesh_linear(self, mcds=mcds):
        # cube coordinates
        ar_m_cube, ar_n_cube, ar_p_cube = mcds.get_mesh(flat=False)
        er_m_cube = set(ar_m_cube.flatten())
        er_n_cube = set(ar_n_cube.flatten())
        er_p_cube = set(ar_p_cube.flatten())
        # linear coordiantes
        lar_voxel = mcds.get_mesh_linear()
        assert (str(type(lar_voxel)) == "<class 'list'>") and \
               (len(lar_voxel) == 3) and \
               (str(type(lar_voxel[0])) == "<class 'numpy.ndarray'>") and \
               (str(lar_voxel[0].dtype) == "float64") and \
               (set(lar_voxel[0]) == er_m_cube) and \
               (lar_voxel[0].shape == (1331,)) and \
               (str(type(lar_voxel[1])) == "<class 'numpy.ndarray'>") and \
               (str(lar_voxel[1].dtype) == "float64") and \
               (set(lar_voxel[1]) == er_n_cube) and \
               (lar_voxel[1].shape == (1331,)) and \
               (str(type(lar_voxel[2])) == "<class 'numpy.ndarray'>") and \
               (str(lar_voxel[2].dtype) == "float64") and \
               (set(lar_voxel[2]) == er_p_cube) and \
               (lar_voxel[2].shape == (1331,))

    def test_mcds_get_mesh_spacing(self, mcds=mcds):
        lr_spacing = mcds.get_mesh_spacing()
        assert lr_spacing == [30.0, 20.0, 10.0]

    def test_mcds_get_voxel_volume(self, mcds=mcds):
        r_volume = mcds.get_voxel_volume()
        assert r_volume == 6000.0

    def test_mcds_get_voxel_ijk(self, mcds=mcds):
        li_voxel_0 = mcds.get_voxel_ijk(x=0, y=0, z=0)
        li_voxel_1 = mcds.get_voxel_ijk(x=15, y=10, z=5)
        li_voxel_2 = mcds.get_voxel_ijk(x=30, y=20, z=10)
        assert (li_voxel_0 == [0, 0, 0]) and \
               (li_voxel_1 == [1, 1, 1]) and \
               (li_voxel_2 == [2, 2, 2])

    def test_mcds_is_in_mesh(self, mcds=mcds):
        assert mcds.is_in_mesh(x=42, y=42, z=42, b_break=False) and \
               not mcds.is_in_mesh(x=-42, y=-42, z=-42, b_break=False)

    ## micro environment related functions
    def test_mcds_get_substrate_names(self, mcds=mcds):
        ls_substrate = mcds.get_substrate_names()
        assert ls_substrate == ['immunostimulatory_factor', 'oxygen']

    def test_mcds_get_substrate_df(self, mcds=mcds):
        df_substrate = mcds.get_substrate_df()
        assert (str(type(df_substrate)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_substrate.shape == (4, 1))

    def test_mcds_get_concentrations(self, mcds=mcds):
        ar_conc = mcds.get_concentrations(species_name='oxygen', z_slice=None)
        assert (str(type(ar_conc)) == "<class 'numpy.ndarray'>") and \
               (ar_conc.shape == (11, 11, 11))

    def test_mcds_get_concentrations_zslice(self, mcds=mcds):
        ar_conc = mcds.get_concentrations(species_name='oxygen', z_slice=-5)
        assert (str(type(ar_conc)) == "<class 'numpy.ndarray'>") and \
               (ar_conc.shape == (11, 11))

    def test_mcds_get_concentrations_at(self, mcds=mcds):
        ar_conc = mcds.get_concentrations_at(x=0, y=0, z=0)
        assert (str(type(ar_conc)) == "<class 'numpy.ndarray'>") and \
               (ar_conc.shape == (2,))

    def test_mcds_get_concentrations_df(self, mcds=mcds):
        df_conc = mcds.get_concentrations_df(z_slice=None)
        assert (str(type(df_conc)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_conc.shape == (1331, 8))

    def test_mcds_get_concentrations_df_zslice(self, mcds=mcds):
        df_conc = mcds.get_concentrations_df(z_slice=-5)
        assert (str(type(df_conc)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_conc.shape == (121, 8))

    ## cell realted functions
    def test_mcds_get_cell_variables(self, mcds=mcds):
        ls_variable = mcds.get_cell_variables()
        assert (str(type(ls_variable)) == "<class 'list'>") and \
               (len(ls_variable) == 97) and \
               (ls_variable[0] == 'ID')

    def test_mcds_get_cell_df(self, mcds=mcds):
        df_cell = mcds.get_cell_df()
        assert (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_cell.shape == (20460, 110))

    def test_mcds_get_cell_df_at(self, mcds=mcds):
        df_cell = mcds.get_cell_df_at(x=0, y=0, z=0)
        assert (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_cell.shape == (5, 110))

    ## unit related functions
    def test_mcds_get_unit_df(self, mcds=mcds):
        df_unit = mcds.get_unit_df()
        assert (str(type(df_unit)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_unit.shape == (105, 1))
