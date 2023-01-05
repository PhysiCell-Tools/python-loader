# pcDataLoader Tutorial Man Page

Please install the latest version of pcDataLoader, as descripend in the [HowTo](https://github.com/elmbeech/pcDataLoader/blob/v3/man/HOWTO.md) section.\
If you are new to pcDataLoader, install from branch v3.\
Branch v2 exist mainly for backwards compatibility.


## Tutorial - branch v2
The original python-loader tutorial can be found here.
+ http://www.mathcancer.org/blog/python-loader/


## Tutorial - branch v3

### Introduction: PhysiCell's MultiCellDigital Snapshot (MCDS) Anatomy

Each time PhysiCell’s internal time tracker passes a time step where data is to be saved, it generates a number of files of various types.\
Each of these files will have a number at the end that indicates where it belongs in the sequence of outputs.\
All of the files from the first round of output will end in 00000000.* and the second round will be 00000001.* and so on.\
Have a look at this [PhysiCell data time series](https://github.com/elmbeech/pcDataLoader/tree/v3/pcDataLoader/data_timeseries_2d).

Let’s say we captured data every simmulation time hour, and we’re interested in the set of output, half a day through the run, the 13th set of output files.\
The files we care about most from this set consists of:

+ **output00000012.xml**: This file is the main organizer of the data.
    It contains an overview of the data stored in the MultiCellDS as well as some actual data including:\
    metadata (MultiCellDS version, PhysiCell or BioFVM version, simulation time, runtime, and processing time stamp),\
    coordinates for the computational domain (mesh),\
    parameters for diffusing substrates in the microenvironment (variables),\
    column labels and units for the cell data (cell_population),\
    file names for the files that contain microenvironment and cell data at this time step (mat and possibly graph.txt files),
+ **output00000012_cells.mat**: This is a MATLAB matrix file that contains all of the tracked information about the individual cells in the model.
    It tells us things like the cells’ position, volume, secretion, cell cycle status, and user-defined cell parameters.
+ **output00000012_microenvironment0.mat**: This is a MATLAB matrix file that contains all of the data about the microenvironment at this time step.


### Loading an MCDS into Python3
In this paragraph we will load the pcDataLoder library and data snapshot 00000012, described above, from [data timeseries 2d](https://github.com/elmbeech/pcDataLoader/tree/v3/pcDataLoader/data_timeseries_2d) from the test dataset.
(There is no need to extra install the test data set.
In fact, both test datasets are already installed.
They ship with the pip3 pcDataLoader installation.)

```python
import pathlib  # library to locate the test data
import pcDataLoader as pc  # the PhysiCell data loader library

s_path = str(pathlib.Path(pc.__file__).parent.joinpath('data_timeseries_2d'))  # local path to the installed test data set
s_file = 'output00000012.xml'  # the snapshot we want to analyse
s_pathfile = f'{s_path}/{s_file}'
print('mcds snapshot xml:', s_pathfile)

# load mcds - multi cell digital snapshot - object
mcds = pc.pyMCDS(s_pathfile)  # loads whole snapshot: the xml and all realated mat and graph files.
```

Side notes: in general, unix and windows path notation will work.\
Additionaly, legacy way of data loaing works too.
```python
# legacy way of loading an mcds object
mcds = pc.pyMCDS('output00000012.xml', s_path)
```

By default all data realted to the snapshot is loaded.\
For speed and less memory usage, it is however possible to only load the essential (xml and cell mat data), and exclude microenviroment and graph data loading.
```python
# fine tuned way of loading an mcds object
mcds = pc.pyMCDS(s_pathfile, graph=False, microenv=False)
```

### The Loaded Data Structure

All loaded  data lievs in `mcds.data` dictionary.\
As in **version 1**, we’ve tried to keep everything organized inside of this dictionary.\
Regarding **version 1**, the structure has slightly change.\
However, in **version 3**, all data is accessible by functions, thus there should be no need to fetch data directely form the `mcds.data` dictionary!\
Anyhow, let’s take a look at what we actually have in here.

![mcds.data dictionary blue print](img/pcdataloader_data_dictionary_v3.0.0.png)

```python
# main data branches
sorted(mcds.data.keys())  # metadata, mesh, substrate (continuum_variables), and agent (discrete_cells).

# metadata
sorted(mcds.data['metadata'].keys())  # multicellds version, physicell version, simulation time, runtime, time stamp, time unit, and spatial unit.

# mesh
sorted(mcds.data['mesh'].keys())  # voxel (ijk), mesh (nmp), and position (xyz) range, axis, coordinate, grid objects, and voxel volume.

# microenvironment
sorted(mcds.data['continuum_variables'].keys())  # list of all processed substrat, e.g. oxygen.
sorted(mcds.data['continuum_variables']['oxygen'].keys())  # substrat related data values, unit, diffusion coefficient, and decay rate.

# cell
sorted(mcds.data['discrete_cells'].keys())  # data, units, and graph dictionaries.
sorted(mcds.data['discrete_cells']['data'].keys()  # all cell realted, tracked data.
sorted(mcds.data['discrete_cells']['unit'].keys()  # all units from the cell realted, tracked data.
sorted(mcds.data['discrete_cells']['graph'].keys())  # neighbor_cells and attached_cells graph dictionaries.
```

### Accessing the Loaded Data by the pyMCDS Class Functions
Once agin, for the ones in the backrow, in version 3, all data is accessible by functions.\
There should be no need to fetch data directely form the `mcds.data` dictionary of dictionaries.\
We will explor this functions in this paragaraph.

#### Metadata

Fetch the data's MultiCellDS version, and the PhysiCell version the data was generated.
```python
mcds.get_multicellds_version()  # will retuns a string like MultiCellDS_2 or MultiCellDS_0.5
mcds.get_physicell_version()  # will returns a string like PhysiCell_1.10.4 or BioFVM_1.1.7
```

Fetch simulation time, runtime, and time stamp when the data was processed.
```python
mcds.pyMCDS.get_time()   # will retun a float value like 720.0
mcds.pyMCDS.get_runtime()  # will retun a float value like 15.596373
mcds.get_timestamp()  # will retun a sting like 2022-10-19T01:12:01Z
```

#### Mesh Data

Let's start by the voxel.\
It is common, but not necessary, that the voxels width, heigth, and depth is the same.\
However, you will find that in this test dataset this is not the case. \
In the related `PhysiCell_settings.xml` file the voxel was defined by 30[nm] high, 20[nm] wide, and 10[nm] deep.\
We can retrive the voxel spacing value and voxel volume. \
Additionaly, we can retrive the mesh spacing.\
You will notice that voxel and mesh spacing values differ.
This is because this data set is from a 2D simulation. Because of that, the mesh depth is set to 1.
In 3D simulation data voxel and mesh spacing will be the same.
```python
mcds.get_mesh_spacing()  # [30.0, 20.0, 1]
mcds.get_voxel_spacing() # [30.0, 20.0, 10.0]
mcds.get_voxel_volume()  # 30[nm] * 20[nm] * 10[nm] = 6000[nm\*\*3]
```

Since version 3, coordinate labels are distinct:
+ **x,y,z:** cell position coordinates are real values.
+ **m,n,p:** mesh center coordinates are real values.
+ **i,j,k:** voxel index coordinates are unsigned interger values.

For a better understanding, let's fetch the start and stop range values for each coordinate type.
Unlike in the python range function, the start and stop value is inclusive.
```python
mcds.get_voxel_ijk_range()  # [(0, 10), (0, 10), (0, 0)]
mcds.get_mesh_mnp_range()  # [(-15.0, 285.0), (-10.0, 190.0), (0.0, 0.0)]
mcds.get_xyz_range()  # [(-30.0, 300.0), (-20.0, 200.0), (-5.0, 5.0)]
```
The domain benchmarks are:
+ voxel index values stretch from 0 to 10 longitude (i), 0 to 10 latitude (j), the 2D domain is only 1 voxel deep (k).
+ mesh center values strech from -15[nm] to 285[nm] longitude (m), -10[nm] to 190[nm] latitude (n), and the depth mesh center (p) is at zero.
+ cell can hold a positioned between -30[nm] and 300[nm] longitude (x), -20[nm] and 200[nm] latitude (y), and -5[nm] and 5[nm] depth (z).

For voxel and  mesh centers we can fetch the axis tick lists.
Each of this function will return a list of 3 numpy arrays, orderd ijk or mnp, respective.
```python
mcds.get_voxel_ijk_axis()
mcds.get_mesh_mnp_axis()
```

We can even retive all mnp mesh center coordinate triplets, and the meshgrids in 2D or 3D form it self.
The coordinate triplets are returned as numpy array
The mesh grids are retunred as 3 and 4 way tensors numpy array, but showen below, it is easy to subseted these tensors.
```python
# all mnp coordinates
mnp_coordinats = mcds.get_mesh_coordinate()  # numpy array with shape (3, 121)
mnp_coordinats.shape  # (3, 121)

# 2D mesh grid
mn_meshgrid = mcds.get_mesh_2D()  # numpy array with shape (2, 11, 11)
mn_meshgrid.shape # (2, 11, 11)
m_meshgrid = mn_meshgrid[0]  # or explicite mn_meshgrid[0,:,:]
n_meshgrid = mn_meshgrid[1]  # or explicite mn_meshgrid[1,:,:]
m_meshgrid.shape  # (11, 11)
n_meshgrid.shape  # (11, 11)

# 3D mesh grid
mnp_meshgrid = mcds.get_mesh()  # numpy array with shape (3, 11, 11, 1)
mnp_meshgrid.shape # (3, 11, 11, 1)
m_meshgrid = mnp_meshgrid[0]  # or explicite mnp_meshgrid[0,:,:,:]
n_meshgrid = mnp_meshgrid[1]  # or explicite mnp_meshgrid[1,:,:,:]
p_meshgrid = mnp_meshgrid[2]  # or explicite mnp_meshgrid[2,:,:,:]
m_meshgrid.shape  # (11, 11, 1)
n_meshgrid.shape  # (11, 11, 1)
p_meshgrid.shape  # (11, 11, 1)
```

Furthmore, there are two helper function.
One to figure out of a particular x,y,z coordinat is still in side the mesh,
and other one to translate an x,y,z coordnate into i,j,k voxel indices, that by default as well checks, if the given coordinat is in the mesh.
```python
# is given x,y,z is in the mesh?
mcds.is_in_mesh(x=0, y=0, z=0)  # True
mcds.is_in_mesh(x=111, y=22, z=-5)  # True
mcds.is_in_mesh(x=111, y=22, z=-5.1)  # False and Warning @ pyMCDS.is_in_mesh : z = -5.1 out of bounds: z-range is (-5.0, 5.0).

# translate x,y,z into i,j,k!
mcds.get_voxel_ijk(x=0, y=0, z=0)  # [0,0,0]
mcds.get_voxel_ijk(x=111, y=22, z=-5)  # [4, 2, 0]
mcds.get_voxel_ijk(x=111, y=22, z=-5.1)  # None and Warning @ pyMCDS.is_in_mesh : z = -5.1 out of bounds: z-range is (-5.0, 5.0).
```

#### Microenvironment (Continuum Variables) Data

Let's have a look at the substrate.
We can retrieve an alphabetically ordered list of all the substrates processes in the simulation.
In the loaded data set only one substarte, oxygen, was part of the simulation.
```python
mcds.get_substrate_names()  # ['oxygen']
```

We can retrieve a pandas dataframe withe the parameters (decay\_rate, diffusion\_coefficient) that were set for in each substrate.
```python
df = mcds.get_substrate_df()
df.head()
```

Regaridng the concentrations, we can retieve:
+ a **numpy array** of all substratd concentrations at a particular xyz coordinate, ordered alphabetically by substrate name, like the list retieved by the get\_substrate\_names function.
+ substrate specific 3D or 2D **meshgrid numpy arrays**.
  To get a 2D meshgrids you can slice though any z stack value, the fuction will always pick the closest mesh coordinate, the smaller coordinate, if you hit the saddlepont between two voxels.
+ a **pandas dataframe** with voxel ijk coordinate values, mesh center mnp coordinate values, and concentrations values for all substrates.
```python
# all contration values at particular coordinate
mcds.get_concentration_at(x=0, y=0, z=0)  # array([34.4166271])
mcds.get_concentration_at(x=111, y=22, z=-5)  # array([18.80652216])
mcds.get_concentration_at(x=111, y=22, z=-5.1)  # None and Warning @ pyMCDS.is_in_mesh : z = -5.1 out of bounds: z-range is (-5.0, 5.0).

# concentration meshgrid for a particular substrate
oxygen_3d = mcds.get_concentration('oxygen')
oxygen_2d = mcds.get_concentration('oxygen', z_slice=0)
oxygen_3d.shape  # (11, 11, 1)
oxygen_2d.shape  # (11, 11)

# dataframe with complete voxel cordinate (ijk), mesh coordiante (mnp), and substarte concentration mapping.
df = mcds.get_concentration_df()
df.head()
```

#### Cell (Discrete Cells) Data

Now, lets have a look at all variables we track from each agent.\
It is quite a list.\
This list is ordered alphabetically, except the for first entry, which always is ID.
```python
mcds.get_cell_variables()  # ['ID', 'attachment_elastic_constant', ..., 'velocity_z']
len(mcds.get_cell_variables())  # 77
```

The data itslef we can retireve as dataframe.\
There exit two function to do so.\
One that retriefs data from all agents in the whole domain, the oder retrievs only data fromm the agenetis in a particular voxel, specified by the xyz coordinate.\
Please note, this dataframes not only hold the exact x,y,z coordinate and all dicret variables mentioned above,
they list additionally all agent surrounding substrate concentartion, the containg i,j,k voxel coordinate and m,n,p mesh center coordinate.

```python
# data from all agent in the domain.
df = mcds.get_cell_df()
df.shape  # (992, 87)  this mean: 992 agent in the whole domain, 87 tracked variables.
df.info()
df.head()

# data from all agents in the xyz specified voxel.
df = mcds.get_cell_df_at(x=0,y=0,z=0)
df.shape  # (4, 87)

df = mcds.get_cell_df_at(x=111,y=22,z=-5)
df.shape  # (3, 87)

df = mcds.get_cell_df_at(x=111,y=22,z=-5.1)  # None and Warning @ pyMCDS.is_in_mesh : z = -5.1 out of bounds: z-range is (-5.0, 5.0).
```

Since lately, PysiCell tracks for each cell, if it touches other cells.\
This data is stored in two dictionaries of sets of integers which link the corresponding cell ids.\
We have here a closer look at the neighbour graph because the attached graph is in this particular study not really intresting.
```python
# attached graph
graph = mcds.get_attached_graph_dict()
len(graph)  # 992

# neighbor graph
graph = mcds.get_neighbor_graph_dict()
len(graph)  # 992
graph.keys()  # dict_keys([0, 1, ..., 993])
graph[0]  # {1, 31, 33, 929, 935, 950}
```

#### Units

Finally, it is possible to retrieve a dataframe that lists all units from all tracked variables, from metadata, mesh, continuum\_variables, and discrete\_cells.
```python
df = mcds.get_unit_df()
df.shape  # (82, 1)
df.columns  # Index(['unit'], dtype='object')
df.index  # Index(['attachment_elastic_constant', 'attachment_rate', ..., 'velocity_z'], dtype='object', name='parameter')
df.head()
```


### Working With MCDS Time Series in Python3


################
## HERE I AM ##
###############

from first implementation:
In this tutorial we will load the 24[h] time step data snapshot from the 3D cancer-immune-sample project, provided as test data with this libarary.


One of the big advantages of working with PhysiCell data in python is that we have access to its plotting tools. For the sake of example let’s plot the partial pressure of oxygen throughout the computational domain along the z=0 plane. Once we’ve loaded our data by initializing a pyMCDS object, we can work entirely within python to produce the plot.
?

from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt

# load data
mcds = pyMCDS('output00003696.xml', '../output')

# Set our z plane and get our substrate values along it
z_val = 0.00
plane_oxy = mcds.get_concentrations('oxygen', z_slice=z_val)

# Get the 2D mesh for contour plotting
xx, yy = mcds.get_mesh()

# We want to be able to control the number of contour levels so we
# need to do a little set up
num_levels = 21
min_conc = plane_oxy.min()
max_conc = plane_oxy.max()
my_levels = np.linspace(min_conc, max_conc, num_levels)

# set up the figure area and add data layers
fig, ax = plt.subplot()
cs = ax.contourf(xx, yy, plane_oxy, levels=my_levels)
ax.contour(xx, yy, plane_oxy, color='black', levels = my_levels,
           linewidths=0.5)

# Now we need to add our color bar
cbar1 = fig.colorbar(cs, shrink=0.75)
cbar1.set_label('mmHg')

# Let's put the time in to make these look nice
ax.set_aspect('equal')
ax.set_xlabel('x (micron)')
ax.set_ylabel('y (micron)')
ax.set_title('oxygen (mmHg) at t = {:.1f} {:s}, z = {:.2f} {:s}'.format(
                                        mcds.get_time(),
                                        mcds.data['metadata']['time_units'],
                                        z_val,
                                        mcds.data['metadata']['spatial_units'])

plt.show()
oxygen partial pressures over z=0

Adding a cells layer

We can also use pandas to do fairly complex selections of cells to add to our plots. Below we use pandas and the previous plot to add a cells layer.
?

from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt

# load data
mcds = pyMCDS('output00003696.xml', '../output')

# Set our z plane and get our substrate values along it
z_val = 0.00
plane_oxy = mcds.get_concentrations('oxygen', z_slice=z_val)

# Get the 2D mesh for contour plotting
xx, yy = mcds.get_mesh()

# We want to be able to control the number of contour levels so we
# need to do a little set up
num_levels = 21
min_conc = plane_oxy.min()
max_conc = plane_oxy.max()
my_levels = np.linspace(min_conc, max_conc, num_levels)

# get our cells data and figure out which cells are in the plane
cell_df = mcds.get_cell_df()
ds = mcds.get_mesh_spacing()
inside_plane = (cell_df['position_z'] < z_val + ds) \ & (cell_df['position_z'] > z_val - ds)
plane_cells = cell_df[inside_plane]

# We're going to plot two types of cells and we want it to look nice
colors = ['black', 'grey']
sizes = [20, 8]
labels = ['Alive', 'Dead']

# set up the figure area and add microenvironment layer
fig, ax = plt.subplot()
cs = ax.contourf(xx, yy, plane_oxy, levels=my_levels)

# get our cells of interest
# alive_cells = plane_cells[plane_cells['cycle_model'] < 6]
# dead_cells = plane_cells[plane_cells['cycle_model'] > 6]
# -- for newer versions of PhysiCell
alive_cells = plane_cells[plane_cells['cycle_model'] < 100]
dead_cells = plane_cells[plane_cells['cycle_model'] >= 100]

# plot the cell layer
for i, plot_cells in enumerate((alive_cells, dead_cells)):
    ax.scatter(plot_cells['position_x'].values,
            plot_cells['position_y'].values,
            facecolor='none',
            edgecolors=colors[i],
            alpha=0.6,
            s=sizes[i],
            label=labels[i])

# Now we need to add our color bar
cbar1 = fig.colorbar(cs, shrink=0.75)
cbar1.set_label('mmHg')

# Let's put the time in to make these look nice
ax.set_aspect('equal')
ax.set_xlabel('x (micron)')
ax.set_ylabel('y (micron)')
ax.set_title('oxygen (mmHg) at t = {:.1f} {:s}, z = {:.2f} {:s}'.format(
                                        mcds.get_time(),
                                        mcds.data['metadata']['time_units'],
                                        z_val,
                                        mcds.data['metadata']['spatial_units'])
ax.legend(loc='upper right')

plt.show()

adding a cell layer to the oxygen plot

Future Direction

The first extension of this project will be timeseries functionality. This will provide similar data loading functionality but for a time series of MultiCell Digital Snapshots instead of simply one point in time.

