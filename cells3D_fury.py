import numpy as np
from fury import window, actor, ui
import itertools
import math
from pyMCDS import pyMCDS

mcds = pyMCDS('output00003696.xml')  # replace this filename with yours
print(mcds.get_time())
print('\n------- cell_vars:')
cell_vars = mcds.get_cell_variables()
print(cell_vars)
#['ID', 'position_x', 'position_y', 'position_z', 'total_volume', 'cell_type', 'cycle_model', 'current_phase', 'elapsed_time_in_phase', 'nuclear_volume', 'cytoplasmic_volume', 'fluid_fraction', 'calcified_fraction', 'orientation_x', 'orientation_y', 'orientation_z', 'polarity', 'migration_speed', 'motility_vector_x', 'motility_vector_y', 'motility_vector_z', 'migration_bias', 'motility_bias_direction_x', 'motility_bias_direction_y', 'motility_bias_direction_z', 'persistence_time', 'motility_reserved', 'oncoprotein', 'elastic_coefficient', 'kill_rate', 'attachment_lifetime', 'attachment_rate']

print('\n------- cell_data:')
cell_data = mcds.get_cell_df()
print(cell_data)
print('cell_type min = ',cell_data['cell_type'].min())   # = 0.0
print('cell_type max = ',cell_data['cell_type'].max())   # = 1.0
#In [7]: cell_data['position_x'].min()
#Out[7]: -750.0005921113642
#  BEWARE - bogus max!
#In [6]: cell_data['position_x'].max()
#Out[6]: 2.8091642052341613e+229

# double cell_radius = cell_defaults.phenotype.geometry.radius; # int xc=0,yc=0,zc=0;
cell_radius = 8.412710547954228

xyz = np.empty((0,3))
#colors = np.ones((num_cells,4))  # white [0-2]; opaque [3]
colors = np.empty((0,4))
print('supposed len of cells = ',len(cell_data['position_x']))
#xyz = np.empty((len(cell_data['position_x']),3))

idx_cell = 0
for idx in range(len(cell_data['position_x'])):
	x = cell_data['position_x'][idx]
	y = cell_data['position_y'][idx]
	z = cell_data['position_z'][idx]
	if (abs(x) < 800. and abs(y) < 800. and abs(y) < 800.):
		if (x > 0.):
			xyz = np.append(xyz, np.array([[x,y,z]]), axis=0)
			if (cell_data['cell_type'][idx] > 0.0):
				colors = np.append(colors, np.array([[1,0,0,1]]), axis=0)  # red
			else:
				colors = np.append(colors, np.array([[0,1,1,1]]), axis=0)  # cyan
			idx_cell += 1
	else:
		print(x,y,z)

num_cells = idx_cell
print('actual len of cells = ',num_cells)
# default color
#colors = np.ones((num_cells,4))  # white [0-2]; opaque [3]
#colors[:,0] = 0  # make them cyan
#for idx in range(num_cells):

scene = window.Scene()
sphere_actor = actor.sphere(centers=xyz,
                            colors=colors,
                            radii=cell_radius)
scene.add(sphere_actor)
showm = window.ShowManager(scene,
                           size=(900, 768), reset_camera=False,
                           order_transparent=False)
showm.initialize()
showm.start()