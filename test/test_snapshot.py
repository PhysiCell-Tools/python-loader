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

# load physicell data
class TestPyMCDS(object):
    ''' test for pc.pyMCDS '''

    def test_pyMCDS(self):
        # load physicell data
        print(f'process: pc.pyMCDS(xml_file={s_pathfile}) ...')
        mcds = pc.pyMCDS(xml_file=s_pathfile)
        assert str(type(mcds)) == "<class 'pcDataLoader.pyMCDS.pyMCDS'>"

        # load physicell data
        print(f'process: pc.pyMCDS(xml_file=s_file, output_path=s_path, microenv=False) ...')
        mcds = pc.pyMCDS(xml_file=s_file, output_path=s_path, microenv=False)
        assert str(type(mcds)) == "<class 'pcDataLoader.pyMCDS.pyMCDS'>"

        # load physicell data
        print(f'process: pc.pyMCDS(xml_file=s_file, output_path=s_path, microenv=True) ...')
        mcds = pc.pyMCDS(xml_file=s_file, output_path=s_path, microenv=True)
        assert str(type(mcds)) == "<class 'pcDataLoader.pyMCDS.pyMCDS'>"

# commands to extract basic information
#print(mcds1.get_time())
#print(mcds1.get_cell_variables())
#print(mcds1.get_substrate_names())

# commnds to extract pandas compatible output
#print(mcds1.get_cell_df())
#print(mcds1.get_cell_df_at(x=39, y=83, z=0))

# commands to extract mesh related output
#print(mcds1.get_mesh_spacing())
#print(mcds1.get_mesh())
#print(mcds1.get_2D_mesh())
#print(mcds1.get_linear_voxels())
#print(mcds1.get_containing_voxel_ijk(x=0,y=0,z=0))
#print(mcds1.get_concentrations('quorum'))
