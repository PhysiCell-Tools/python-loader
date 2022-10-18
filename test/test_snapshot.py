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
class TestPyMcdsShortcut(object):
    ''' test for pc.pyMCDS data loader shortcut '''

    def test_pyMCDS(self):
        # load physicell data shortcut
        print(f'process: pc.pyMCDS(xml_file={s_pathfile_3d}) ...')
        mcds = pc.pyMCDS(xml_file=s_pathfile_3d)
        assert str(type(mcds)) == "<class 'pcDataLoader.pyMCDS.pyMCDS'>"

# load physicell data with microenvironment
class TestPyMcdsMicroenvTrue(object):
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
               # x
               (str(type(lar_mesh[0])) == "<class 'numpy.ndarray'>") and \
               (str(lar_mesh[0].dtype) == "float64" and \
               (lar_mesh[0].shape == (75, 75, 75)) and \
               # y
               (str(type(lar_mesh[1])) == "<class 'numpy.ndarray'>") and \
               (str(lar_mesh[1].dtype) == "float64" and \
               (lar_mesh[1].shape == (75, 75, 75)) and \
               # z
               (str(type(lar_mesh[2])) == "<class 'numpy.ndarray'>") and \
               (str(lar_mesh[2].dtype) == "float64" and \
               (lar_mesh[2].shape == (75, 75, 75))

    def test_mcds_get_mesh_flat_true(self, mcds=mcds):
        print(f'process: mcds.get_mesh(flat=True) ...')
        lar_mesh = mcds.get_mesh(flat=True)
        assert (str(type(lar_mesh)) == "<class 'list'>") and \
               (len(lar_mesh) == 2) and \
               # x
               (str(type(lar_mesh[0])) == "<class 'numpy.ndarray'>") and \
               (str(lar_mesh[0].dtype) == "float64" and \
               (lar_mesh[0].shape == (75, 75)) and \
               # y
               (str(type(lar_mesh[1])) == "<class 'numpy.ndarray'>") and \
               (str(lar_mesh[1].dtype) == "float64" and \
               (lar_mesh[1].shape == (75, 75))

    def test_mcds_get_mesh_2d(self, mcds=mcds):
        print(f'process: mcds.get_mesh_2d() ...')
        lar_mesh_falt = mcds.get_mesh(flat=True)
        lar_mesh_2d = mcds.get_mesh_2D()
        assert (str(type(lar_mesh_2d)) == "<class 'list'>") and \
               (len(lar_mesh_2d) == 2) and \
               (lar_mesh_2d[0] == lar_mesh_flat[0]).all() and \
               (lar_mesh_2d[1] == lar_mesh_flat[1]).all()

    def test_mcds_get_linear_voxels(self, mcds=mcds):
        print(f'process: mcds.get_linear_voxels() ...')
        a_voxel = mcds.get_linear_voxels()
        assert (str(type(a_voxel)) == "<class 'numpy.ndarray'>") and \
               (str(a_voxel.dtype) == "float64" and \
               (a_voxel.shape == (3, 421875))

    def test_mcds_get_mesh_spacing(self, mcds=mcds):
        print(f'process: mcds.get_mesh_spacing() ...')
        i_spacing = mcds.get_mesh_spacing()
        assert (str(type(i_spacing)) == "<class 'int'>") and \
               (i_spacing == 20)

    def test_mcds_get_containing_voxel_ijk(self, mcds=mcds):
        print(f'process: mcds.get_containing_voxel_ijk(x=0, y=0, z=0) ...')
        li_voxel = mcds.get_containing_voxel_ijk(x=0, y=0, z=0)
        assert li_voxel == [37, 37, 37]

    # micro environment related functions
    def test_mcds_get_substrate_names(self, mcds=mcds):
        print(f'process: mcds.get_substrate_names() ...')
        ls_substrate = mcds.get_substrate_names()
        assert ls_substrate == ['immunostimulatory factor', 'oxygen']

    def test_mcds_get_concentrations(self, mcds=mcds):
        print(f'process: mcds.get_concentrations(species_name="oxygen", z_slice=None) ...')
        ar_conc = mcds.get_concentrations(species_name='oxygen', z_slice=None)
        assert (str(type(ar_conc)) == "<class 'numpy.ndarray'>") and \
               (ar_conc.shape == (75, 75, 75))

    def test_mcds_get_concentrations_zslice(self, mcds=mcds):
        print(f'process: mcds.get_concentrations(species_name="oxygen", z_slice=0) ...')
        ar_conc = mcds.get_concentrations(species_name='oxygen', z_slice=0)
        assert (str(type(ar_conc)) == "<class 'numpy.ndarray'>") and \
               (ar_conc.shape == (75, 75))

    def test_mcds_get_concentrations_df(self, mcds=mcds):
        print(f'process: mcds.get_concentrations_df(z_slice=None) ...')
        df_conc = mcds.get_concentrations_df(z_slice=None)
        assert (str(type(df_conc)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_conc.shape == (421875, 8))

    def test_mcds_get_concentrations_df_zslice(self, mcds=mcds):
        print(f'process: mcds.get_concentrations_df(z_slice=0) ...')
        df_conc = mcds.get_concentrations_df(z_slice=0)
        assert (str(type(df_conc)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_conc.shape == (5625, 8))

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
               (len(ls_variable) == 32) and \
               (ls_variable[0] == 'ID')

    def test_mcds_get_cell_df(self, mcds=mcds):
        print(f'process: mcds.get_cell_df() ...')
        df_cell = mcds.get_cell_df()
        assert (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_cell.shape == (66978, 37))

    def test_mcds_get_cell_df_at(self, mcds=mcds):
        print(f'process: mcds.get_cell_df_at(x=0, y=0, z=0) ...')
        df_cell = mcds.get_cell_df_at(x=0, y=0, z=0)
        assert (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_cell.shape == (0,37))


class TestPyMcdsMicroenvFalse(object):
    ''' test for pc.pyMCDS data loader microenvironment False '''

    def test_pyMCDS(self):
        # load physicell data
        print(f'process: pc.pyMCDS(xml_file={s_file_3d}, output_path={s_path_3d}, microenv=False) ...')
        mcds = pc.pyMCDS(xml_file=s_file_3d, output_path=s_path_3d, microenv=False)
        assert str(type(mcds)) == "<class 'pcDataLoader.pyMCDS.pyMCDS'>"

        # metadata realted functions
        #mcds.get_time()

        # mesh related functions
        #mcds.get_mesh(flat=False)
        #mcds.get_mesh_2D()
        #mcds.get_linear_voxels()
        #mcds.get_mesh_spacing()
        #mcds.get_containing_voxel_ijk(x, y, z)

        # micro environment related functions
        #mcds.get_substrate_names()
        #mcds.get_concentrations(species_name, z_slice=None)
        #mcds.get_concentrations_df(self, z_slice=None)
        #mcds.get_concentrations_at(self, x, y, z=0)

        # cell realted functions
        #mcds.get_cell_variables()
        #mcds.get_cell_df()
