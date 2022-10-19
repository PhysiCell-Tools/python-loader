#####
# title: test_snapshot.py
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
class TestPyMcdsMicroenvTrue3D(object):
    ''' test for pc.pyMCDS data loader microenvironment True'''
    mcds = pc.pyMCDS(xml_file=s_file_3d, output_path=s_path_3d, microenv=True)

    def test_pyMCDS(self, mcds=mcds):
        # load physicell data
        print(f'process: pc.pyMCDS(xml_file={s_file_3d}, output_path={s_path_3d}, microenv=True) ...')
        assert str(type(mcds)) == "<class 'pcDataLoader.pyMCDS.pyMCDS'>"

    # metadata realted functions
    def test_mcds_get_time(self, mcds=mcds):
        print(f'process: mcds.get_time() ...')
        r_time = mcds.get_time()
        assert r_time == 1440.0

    # mesh related functions
    def test_mcds_get_mesh_flat_false(self, mcds=mcds):
        print(f'process: mcds.get_mesh(flat=False) ...')
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
        print(f'process: mcds.get_mesh(flat=True) ...')
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
        print(f'process: mcds.get_mesh_2d() ...')
        lar_mesh_flat = mcds.get_mesh(flat=True)
        lar_mesh_2d = mcds.get_mesh_2D()
        assert (str(type(lar_mesh_2d)) == "<class 'list'>") and \
               (len(lar_mesh_2d) == 2) and \
               (lar_mesh_2d[0] == lar_mesh_flat[0]).all() and \
               (lar_mesh_2d[1] == lar_mesh_flat[1]).all()

    def test_mcds_get_linear_voxels(self, mcds=mcds):
        print(f'process: mcds.get_linear_voxels() ...')
        # cube coordinates
        ar_m_cube, ar_n_cube, ar_p_cube = mcds.get_mesh(flat=False)
        er_m_cube = set(ar_m_cube.flatten())
        er_n_cube = set(ar_n_cube.flatten())
        er_p_cube = set(ar_p_cube.flatten())
        # linear coordiantes
        lar_voxel = mcds.get_linear_voxels()
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
        print(f'process: mcds.get_mesh_spacing() ...')
        lr_spacing = mcds.get_mesh_spacing()
        assert lr_spacing == [30.0, 20.0, 10.0]

    def test_mcds_get_containing_voxel_ijk(self, mcds=mcds):
        print(f'process: mcds.get_containing_voxel_ijk(x=0, y=0, z=0) ...')
        li_voxel_0 = mcds.get_containing_voxel_ijk(x=0, y=0, z=0)
        li_voxel_1 = mcds.get_containing_voxel_ijk(x=15, y=10, z=5)
        li_voxel_2 = mcds.get_containing_voxel_ijk(x=30, y=20, z=10)
        assert (li_voxel_0 == [0, 0, 0]) and \
               (li_voxel_1 == [1, 1, 1]) and \
               (li_voxel_2 == [2, 2, 2])

    # micro environment related functions
    def test_mcds_get_substrate_names(self, mcds=mcds):
        print(f'process: mcds.get_substrate_names() ...')
        ls_substrate = mcds.get_substrate_names()
        assert ls_substrate == ['immunostimulatory factor', 'oxygen']

    def test_mcds_get_concentrations(self, mcds=mcds):
        print(f'process: mcds.get_concentrations(species_name="oxygen", z_slice=None) ...')
        ar_conc = mcds.get_concentrations(species_name='oxygen', z_slice=None)
        assert (str(type(ar_conc)) == "<class 'numpy.ndarray'>") and \
               (ar_conc.shape == (11, 11, 11))

    def test_mcds_get_concentrations_zslice(self, mcds=mcds):
        print(f'process: mcds.get_concentrations(species_name="oxygen", z_slice=-5) ...')
        ar_conc = mcds.get_concentrations(species_name='oxygen', z_slice=-5)
        assert (str(type(ar_conc)) == "<class 'numpy.ndarray'>") and \
               (ar_conc.shape == (11, 11))

    def test_mcds_get_concentrations_df(self, mcds=mcds):
        print(f'process: mcds.get_concentrations_df(z_slice=None) ...')
        df_conc = mcds.get_concentrations_df(z_slice=None)
        assert (str(type(df_conc)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_conc.shape == (1331, 8))

    def test_mcds_get_concentrations_df_zslice(self, mcds=mcds):
        print(f'process: mcds.get_concentrations_df(z_slice=-5) ...')
        df_conc = mcds.get_concentrations_df(z_slice=-5)
        assert (str(type(df_conc)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_conc.shape == (121, 8))

    def test_mcds_get_concentrations_at(self, mcds=mcds):
        print(f'process: mcds.get_concentrations_at(x=0, y=0, z=0) ...')
        ar_conc = mcds.get_concentrations_at(x=0, y=0, z=0)
        assert (str(type(ar_conc)) == "<class 'numpy.ndarray'>") and \
               (ar_conc.shape == (2,))

    # cell realted functions
    def test_mcds_get_variables(self, mcds=mcds):
        print(f'process: mcds.get_cell_variables() ...')
        ls_variable = mcds.get_cell_variables()
        assert (str(type(ls_variable)) == "<class 'list'>") and \
               (len(ls_variable) == 97) and \
               (ls_variable[0] == 'ID')

    def test_mcds_get_cell_df(self, mcds=mcds):
        print(f'process: mcds.get_cell_df() ...')
        df_cell = mcds.get_cell_df()
        assert (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_cell.shape == (20460, 102))

    def test_mcds_get_cell_df_at(self, mcds=mcds):
        print(f'process: mcds.get_cell_df_at(x=0, y=0, z=0) ...')
        df_cell = mcds.get_cell_df_at(x=0, y=0, z=0)
        assert (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_cell.shape == (5, 102))


class TestPyMcdsMicroenvFalse3D(object):
    ''' test for pc.pyMCDS data loader microenvironment False '''
    mcds = pc.pyMCDS(xml_file=s_file_3d, output_path=s_path_3d, microenv=False)

    def test_pyMCDS(self, mcds=mcds):
        # load physicell data
        print(f'process: pc.pyMCDS(xml_file={s_file_3d}, output_path={s_path_3d}, microenv=False) ...')
        assert str(type(mcds)) == "<class 'pcDataLoader.pyMCDS.pyMCDS'>"

    # metadata realted functions
    def test_mcds_get_time(self, mcds=mcds):
        print(f'process: mcds.get_time() ...')
        r_time = mcds.get_time()
        assert r_time == 1440.0

    # mesh related functions
    def test_mcds_get_mesh_flat_false(self, mcds=mcds):
        print(f'process: mcds.get_mesh(flat=False) ...')
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
        print(f'process: mcds.get_mesh(flat=True) ...')
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
        print(f'process: mcds.get_mesh_2d() ...')
        lar_mesh_flat = mcds.get_mesh(flat=True)
        lar_mesh_2d = mcds.get_mesh_2D()
        assert (str(type(lar_mesh_2d)) == "<class 'list'>") and \
               (len(lar_mesh_2d) == 2) and \
               (lar_mesh_2d[0] == lar_mesh_flat[0]).all() and \
               (lar_mesh_2d[1] == lar_mesh_flat[1]).all()

    def test_mcds_get_linear_voxels(self, mcds=mcds):
        print(f'process: mcds.get_linear_voxels() ...')
        # cube coordinates
        ar_m_cube, ar_n_cube, ar_p_cube = mcds.get_mesh(flat=False)
        er_m_cube = set(ar_m_cube.flatten())
        er_n_cube = set(ar_n_cube.flatten())
        er_p_cube = set(ar_p_cube.flatten())
        # linear coordiantes
        lar_voxel = mcds.get_linear_voxels()
        assert (str(type(lar_voxel)) == "<class 'list'>") and \
               (len(lar_voxel) == 3) and \
               (str(type(lar_voxel[0])) == "<class 'numpy.ndarray'>") and \
               (str(lar_voxel[0].dtype) == "float64") and \
               (lar_voxel[0].shape == (1331,)) and \
               (set(lar_voxel[0]) == er_m_cube) and \
               (str(type(lar_voxel[1])) == "<class 'numpy.ndarray'>") and \
               (str(lar_voxel[1].dtype) == "float64") and \
               (lar_voxel[1].shape == (1331,)) and \
               (set(lar_voxel[1]) == er_n_cube) and \
               (str(type(lar_voxel[2])) == "<class 'numpy.ndarray'>") and \
               (str(lar_voxel[2].dtype) == "float64") and \
               (lar_voxel[2].shape == (1331,)) and \
               (set(lar_voxel[2]) == er_p_cube)

    def test_mcds_get_mesh_spacing(self, mcds=mcds):
        print(f'process: mcds.get_mesh_spacing() ...')
        lr_spacing = mcds.get_mesh_spacing()
        assert lr_spacing == [30.0, 20.0, 10.0]

    def test_mcds_get_containing_voxel_ijk(self, mcds=mcds):
        print(f'process: mcds.get_containing_voxel_ijk(x=0, y=0, z=0) ...')
        li_voxel_0 = mcds.get_containing_voxel_ijk(x=0, y=0, z=0)
        li_voxel_1 = mcds.get_containing_voxel_ijk(x=15, y=10, z=5)
        li_voxel_2 = mcds.get_containing_voxel_ijk(x=30, y=20, z=10)
        assert (li_voxel_0 == [0, 0, 0]) and \
               (li_voxel_1 == [1, 1, 1]) and \
               (li_voxel_2 == [2, 2, 2])

    # micro environment related functions
    # mcds.get_substrate_names # nop
    # mcds.get_concentrations # nop
    # mcds.get_concentrations_df # nop
    # mcds.get_concentrations_at # nop

    # cell realted functions
    def test_mcds_get_variables(self, mcds=mcds):
        print(f'process: mcds.get_cell_variables() ...')
        ls_variable = mcds.get_cell_variables()
        assert (str(type(ls_variable)) == "<class 'list'>") and \
               (len(ls_variable) == 97) and \
               (ls_variable[0] == 'ID')

    def test_mcds_get_cell_df(self, mcds=mcds):
        print(f'process: mcds.get_cell_df() ...')
        df_cell = mcds.get_cell_df()
        assert (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_cell.shape == (20460, 102))

    def test_mcds_get_cell_df_at(self, mcds=mcds):
        print(f'process: mcds.get_cell_df_at(x=0, y=0, z=0) ...')
        df_cell = mcds.get_cell_df_at(x=0, y=0, z=0)
        assert (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_cell.shape == (5, 102))
