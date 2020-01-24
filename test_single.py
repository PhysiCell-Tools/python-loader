# load library
from pyMCDS import pyMCDS

# load physicell data
mcds1 = pyMCDS('output00000001.xml', 'timeseries_set')
mcds2 = pyMCDS('output00000008.xml', 'timeseries_set')

# commands to extract basic information
print(mcds1.get_time())
print(mcds1.get_cell_variables())
print(mcds1.get_substrate_names())

# commnds to extract pandas compatible output
print(mcds1.get_cell_df())
print(mcds1.get_cell_df_at(x=39, y=83, z=0))

# commands to extract mesh related output
print(mcds1.get_mesh_spacing())
print(mcds1.get_mesh())
print(mcds1.get_2D_mesh())
print(mcds1.get_linear_voxels())
print(mcds1.get_containing_voxel_ijk(x=0,y=0,z=0))
print(mcds1.get_concentrations('quorum'))
