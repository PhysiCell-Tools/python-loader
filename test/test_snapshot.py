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
s_path = str(pathlib.Path(pc.__file__).parent.resolve()/'data_snapshot')
s_file = 'output00003696.xml'
s_pathfile = f'{s_path}/{s_file}'

# load physicell data shortcut
class TestPyMcdsShortcut(object):
    ''' test for pc.pyMCDS data loader shortcut '''

    def test_pyMCDS(self):
        # load physicell data shortcut
        print(f'process: pc.pyMCDS(xml_file={s_pathfile}) ...')
        mcds = pc.pyMCDS(xml_file=s_pathfile)
        assert str(type(mcds)) == "<class 'pcDataLoader.pyMCDS.pyMCDS'>"

# load physicell data with microenvironment
class TestPyMcdsMicroenvTrue(object):
    ''' test for pc.pyMCDS data loader microenvironment True'''
    mcds = pc.pyMCDS(xml_file=s_file, output_path=s_path, microenv=True)

    def test_pyMCDS(self, mcds=mcds):
        # load physicell data
        print(f'process: pc.pyMCDS(xml_file=s_file, output_path=s_path, microenv=True) ...')
        assert str(type(mcds)) == "<class 'pcDataLoader.pyMCDS.pyMCDS'>"

        # metadata realted functions
    def test_mcds_get_time(self, mcds=mcds):
        print(f'process: mcds.get_time() ...')
        r_time = mcds.get_time()
        assert r_time == 30239.999998

        # mesh related functions
    def test_mcds_get_mesh_flat_false(self, mcds=mcds):
        print(f'process: mcds.get_mesh(flat=False) ...')
        o_mesh = mcds.get_mesh(flat=False)
        assert (str(type(o_mesh)) == "<class 'list'>") and \
               (len(o_mesh) == 3) and \
               (str(type(o_mesh[0])) == "<class 'numpy.ndarray'>") and \
               (o_mesh[0].shape == (75, 75, 75)) and \
               (str(type(o_mesh[1])) == "<class 'numpy.ndarray'>") and \
               (o_mesh[1].shape == (75, 75, 75)) and \
               (str(type(o_mesh[2])) == "<class 'numpy.ndarray'>") and \
               (o_mesh[2].shape == (75, 75, 75))

    def test_mcds_get_mesh_flat_true(self, mcds=mcds):
        print(f'process: mcds.get_mesh(flat=True) ...')
        o_mesh = mcds.get_mesh(flat=True)
        assert (str(type(o_mesh)) == "<class 'list'>") and \
               (len(o_mesh) == 2) and \
               (str(type(o_mesh[0])) == "<class 'numpy.ndarray'>") and \
               (o_mesh[0].shape == (75, 75)) and \
               (str(type(o_mesh[1])) == "<class 'numpy.ndarray'>") and \
               (o_mesh[1].shape == (75, 75))

    def test_mcds_get_mesh_2d(self, mcds=mcds):
        print(f'process: mcds.get_mesh_2d() ...')
        o_mesh = mcds.get_mesh_2D()
        assert (str(type(o_mesh)) == "<class 'list'>") and \
               (len(o_mesh) == 2) and \
               (str(type(o_mesh[0])) == "<class 'numpy.ndarray'>") and \
               (o_mesh[0].shape == (75, 75)) and \
               (str(type(o_mesh[1])) == "<class 'numpy.ndarray'>") and \
               (o_mesh[1].shape == (75, 75))

    def test_mcds_get_linear_voxels(self, mcds=mcds):
        print(f'process: mcds.get_linear_voxels() ...')
        a_voxel = mcds.get_linear_voxels()
        assert (str(type(a_voxel)) == "<class 'numpy.ndarray'>") and \
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
    #def test_mcds_get_(self, mcds=mcds):
    #    mcds.get_substrate_names()

    #def test_mcds_get_(self, mcds=mcds):
    #    mcds.get_concentrations(species_name, z_slice=None)

    #def test_mcds_get_(self, mcds=mcds):
    #    mcds.get_concentrations_df(self, z_slice=None)

    #def test_mcds_get_(self, mcds=mcds):
    #    mcds.get_concentrations_at(self, x, y, z=0)


        # cell realted functions
    #def test_mcds_get_(self, mcds=mcds):
    #    mcds.get_cell_variables()

    #def test_mcds_get_(self, mcds=mcds):
    #    mcds.get_cell_df()

    #def test_mcds_get_(self, mcds=mcds):
    #    mcds.get_cell_df_at()


class TestPyMcdsMicroenvFalse(object):
    ''' test for pc.pyMCDS data loader microenvironment False '''

    def test_pyMCDS(self):
        # load physicell data
        print(f'process: pc.pyMCDS(xml_file=s_file, output_path=s_path, microenv=False) ...')
        mcds = pc.pyMCDS(xml_file=s_file, output_path=s_path, microenv=False)
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
