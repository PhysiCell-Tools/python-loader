# PhysiCell Data Loader Python Tutorial

In this chapter, we will load the pcdl library and use its TimeStep class to load the data snapshot 00000012, described above, from [data\_timeseries\_ 2d](https://github.com/elmbeech/physicelldataloader/tree/master/pcdl/data_timeseries_2d) from the 2D time series test dataset.

First, please install the latest version of physicelldataloader (pcdl),
as described in the [HowTo](https://github.com/elmbeech/physicelldataloader/blob/master/man/HOWTO.md) chapter.

And, if not already doen so, have a quick read through the pcdl [background](https://github.com/elmbeech/physicelldataloader/tree/master/man/TUTORIAL_introduction.md) infromation.


## Preparation: install the test data set and extract path where the dataset is stored.

```python
import pcdl  # the physicell data loader library
import pathlib  # the library to locate the test data

# install the test datasets
pcdl.install_data()

# local path to the installed test data set
s_path = str(pathlib.Path(pcdl.__file__).parent.joinpath('data_timeseries_2d'))
s_file = 'output00000012.xml'  # the snapshot we want to analyze
s_pathfile = f'{s_path}/{s_file}'

# output
print('mcds xml path:', s_path)
print('mcds xml file:', s_file)
print('mcds xml path/file:', s_pathfile)
```

For paths, in general, unix (slash) and windows (backslash) notation will work.


## Loading an MCDS Time Step

By default, all data related to the snapshot is loaded.\
For speed and less memory usage, it is however possible to only load the essential (output xml and cell mat data), and exclude microenvironment, graph data, and ID label mapping loading.\
Additionally, it is possible to specify data types for the custom\_data variables, apart from the generic float type, namely: int, bool, and str.

The basic way to load a mcds time step object:
```python
import pcdl  # the physicell data loader library
print('pcdl version:', pcdl.__version__)  # it is easy to figure out which pcdl version you run

mcds = pcdl.TimeStep(s_pathfile)  # loads the whole snapshot: the xml and all related mat and graph files
```

The fine tuned way of loading a mcds time step object:
```python
import pcdl  # the physicell data loader library
print('pcdl version:', pcdl.__version__)  # it is easy to figure out which pcdl version you run

mcds = pcdl.TimeStep(s_pathfile, custom_data_type={}, microenv=True, graph=True, physiboss=True, settingxml='PhysiCell_settings.xml', verbose=True)
```

Verbosity you can tune, evnen after loading the data.
```python3
mcds.set_verbose_false()
mcds.set_verbose_true()
```


## The MCDS Data Structure

All loaded data lives in `mcds.data` dictionary.\
As in original python-loader, we tried to keep everything organized inside this dictionary.\
Regarding the origial python-loader, the structure has slightly changed.\

**In pcdl, all data is accessible by functions. There should be no need to fetch data directly from the `mcds.data` dictionary!**

Anyhow, let's take a look at what we actually have in here.

![mcds.data dictionary blueprint](img/physicelldataloader_data_dictionary_v3.2.2.png)


```python
# main data branches
sorted(mcds.data.keys())  # metadata, mesh, substrate mircoenviroment (continuum_variables), and cell agent (discrete_cells)

# metadata
sorted(mcds.data['metadata'].keys())  # multicellds version, physicell version, simulation time, runtime, time stamp, time unit, spatial unit, and substrate and cell type ID label mappings

# mesh
sorted(mcds.data['mesh'].keys())  # voxel (ijk), mesh (nmp), and position (xyz) range, axis, coordinate, grid objects, and voxel volume

# microenvironment
sorted(mcds.data['continuum_variables'].keys())  # list of all processed substrates, e.g. oxygen
sorted(mcds.data['continuum_variables']['oxygen'].keys())  # substrate related data values, unit, diffusion coefficient, and decay rate

# cell
sorted(mcds.data['discrete_cells'].keys())  # data, units, and graph dictionaries
sorted(mcds.data['discrete_cells']['data'].keys())  # all cell related, tracked data
sorted(mcds.data['discrete_cells']['units'].keys())  # all units from the cell related, tracked data
sorted(mcds.data['discrete_cells']['graph'].keys())  # neighbor_cells and attached_cells graph dictionaries
```

**Once again, loud, for the ones in the back, in pcdl >= version 3, all data is accessible by functions.
There should be no need to fetch data directly from the `mcds.data` dictionary of dictionaries.**

We will explore these functions in the upcomming sections.


## Metadata Related Functions

Fetch the data's MultiCellDS version, and the PhysiCell version the data was generated.

```python
mcds.get_multicellds_version()  # will return a string like MultiCellDS_2 or MultiCellDS_0.5
mcds.get_physicell_version()  # will return a string like PhysiCell_1.10.4 or BioFVM_1.1.7
```

Fetch simulation time, runtime, and time stamp when the data was processed.

```python
mcds.get_time()   # will return a float value like 720.0
mcds.get_runtime()  # will return a float value like 15.596373
mcds.get_timestamp()  # will return a sting like 2022-10-19T01:12:01Z
```

Finally, it is possible to retrieve a dictionary that lists all units from all tracked variables, from metadata, mesh, continuum\_variables, and discrete\_cells.

```python
mcds.get_unit_dict()  # will return a dictionary, which maps the variables to the units specified in the setting.xml.
sorted(mcds.get_unit_dict().items())
```


## Microenvironment Data Related Functions (Continuum Variables)

In the loaded dataset, only one substrate, oxygen, was part of the simulation.\
However, let's have a look at the substrate microenviroment related functions.

We can retrieve an alphabetically ordered list of all substrates processed in the simulation.
(This function might become deprected in a future version of pcdl.)

```python
mcds.get_substrate_list()  # ['oxygen']
```

We can retrieve a dictionary that maps substrate indexes to labels.

```python
mcds.get_substrate_dict()  # {'0': 'oxygen'}

# substrates sorted alpahbetically
sorted(mcds.get_substrate_dict().values())  # ['oxygen']

# substrates sorted by index
ds_substrate = mcds.get_substrate_dict())
[ds_substrate[s_key] for s_key in sorted(ds_substrate, key=int)]  # ['oxygen']
```

We can also retrieve a pandas dataframe with the parameters (decay\_rate, diffusion\_coefficient) that were set for each substrate.

```python
df = mcds.get_substrate_df()
df.head()
```

Regarding the concentrations, we can retrieve a **pandas dataframe** with voxel ijk coordinate values, mesh center mnp coordinate values, and concentrations values for all substrates.

```python
# dataframe with complete voxel coordinate (ijk), mesh coordinate (mnp), and substrate concentration mapping
df = mcds.get_conc_df()
df.head()

# detect substrates that not over whole domain have the same concentration
df = mcds.get_conc_df(values=2)
df.head()  # oxygen concentration varies over the domain
```

Additionally, there is a less often used function to retrieve substrate specific 3D or 2D **meshgrid numpy arrays**.
To get a 2D meshgrids you can slice though any z stack value, the function will always pick the closest mesh center coordinate, the smaller coordinate, if you hit the saddle point between two voxels.
(This function might become deprected in a future version of pcdl.)

```python
# concentration meshgrid for a particular substrate
oxygen_3d = mcds.get_concentration('oxygen')
oxygen_2d = mcds.get_concentration('oxygen', z_slice=0)
oxygen_3d.shape  # (11, 11, 1)
oxygen_2d.shape  # (11, 11)
```

Additionally, there is a less often used functions to retrieve a **numpy array** of all substrate concentrations at a particular xyz coordinate, ordered alphabetically by substrate name, like the list retrieved by the get\_substrate\_names function.
(This function might become deprected in a future version of pcdl as well.)

```python
# all concentration values at a particular coordinate
mcds.get_concentration_at(x=0, y=0, z=0)  # array([34.4166271])
mcds.get_concentration_at(x=111, y=22, z=-5)  # array([18.80652216])
mcds.get_concentration_at(x=111, y=22, z=-5.1)  # None and Warning @ pyMCDS.is_in_mesh : z = -5.1 out of bounds: z-range is (-5.0, 5.0)
```

For substrate concentration visualization **matplotlib contour and contourf plots**,
for any substrate, through any z\_slice can be retrived.\
The mcds.plot_contour function has many parameters to fine tune the plot.
Please have a look at it's docstring in the [REFERENCE.md]()manual to learn more.\
The mcds.plot_contour function output can combine with the mcds.plot_scatter output.
Please have a look at the [TUTORIAL_python3_matplotlib.md]() to learn more.

```python
# contour plot
fig = mcds.plot_contour('oxygen')
fig = mcds.plot_contour('oxygen', z_slice=3.333)
fig.show()
```

For substrate concentration visualization **rectilinear grid vtk files**
for any substrate can be retrived.\
This files can be analysied, for example with the [Paraview](https://en.wikipedia.org/wiki/ParaView) software.
Please have a look at the [TUTORIAL_paraview.md]() to learn more.

```python
mcds.make_conc_vtk()
```


## Cell Data Related Functions (Discrete Cells)

We can retrieve a dictionary that maps cell type indexes to labels.

```python
mcds.get_celltype_dict()  # {'0': 'cancer_cell'}
```


Now, let's have a look at all variables we track from each agent.\
It is quite a list.\
This list is ordered alphabetically, except for the first entry, which always is ID.

```python
mcds.get_celltype_list()  # ['ID', 'attachment_elastic_constant', ..., 'velocity_z']
len(mcds.get_celltype_list())  # 77
```


The data itself we can retrieve as dataframe.\
There exist two functions to do so.\
One that retrieves data from all agents in the whole domain, the other retrieves only data from the agents in a particular voxel, specified by the xyz coordinate.\
Please note, this dataframes not only hold the exact xyz coordinate, and all discrete variables mentioned above, they list additionally all agent surrounding substrate concentrations, the containing ijk voxel coordinate, and the mnp mesh center coordinate.

```python
# data from all agents in the domain
df = mcds.get_cell_df()
df.shape  # (992, 94)  this means: 992 agents in the whole domain, 94 tracked variables
df.info()
df.head()

# data variables that are in all agents the same carry no information
# let's filter for variables that carry at least 2 values
df = mcds.get_cell_df(values=2)
df.shape # (992, 38) this means: 992 agents in the whole domain, 38 tracked variables have more than 2 values
```
BUE 20240808: Data Triage

Cell variables that have no variance, zero entropy, that are in all agent overall time steps in the same state, have always exacted the same value, carry no information.
Similarly, substrates variables that over the whole domain overall time steps have the same concentration are not interesting.

```python
# fetch data
mcdsts = pcdl.TimeSeries(s_path)
```

There are functions to help triage over the entier time series for attributes that more likely might carry information, by checking for variables with variation.
```
# cell data min max values
dl_cell = mcdsts.get_cell_attribute()  # returns a dictionary with all attributes, listing all accessed values
len(dl_cell)  # 84 attributes
dl_cell.keys()  # list attribute names
dl_cell['oxygen']  # list min and max oxygen values found, surrounding a cell, over the whole series

# cell data number of values
di_state = {}
[di_state.update({s_attribute: len(li_state)}) for s_attribute, li_state in mcdsts.get_cell_attribute(allvalues=True).items()]
di_state['oxygen']  # cell surrounding oxygen was found occupying 2388 different values (states) over the whole time series

# substrate data
dl_conc = mcdsts.get_conc_attribute()
dl_conc.keys()  # list attribute names
dl_conc['oxygen']  # list min and max oxygen values found in the domain over the whole series
```
BUE 20240808: Data Triage


BUE 20240808: HELPER FUNCTIONS

```python
# data from all agents in the xyz specified voxel
df = mcds.get_cell_df_at(x=0,y=0,z=0)
df.shape  # (4, 94)

df = mcds.get_cell_df_at(x=111,y=22,z=-5)
df.shape  # (3, 94)

df = mcds.get_cell_df_at(x=111,y=22,z=-5.1)  # None and Warning @ pyMCDS.is_in_mesh : z = -5.1 out of bounds: z-range is (-5.0, 5.0)
```
BUE 20240808: HELPER FUNCTIONS


Additionally, there is a scatter plot function available, shaped for df\_cell dataframe content.
```
# scatter plot
fig = mcds.plot_scatter()  # default focus is cell_type and z_slice=0
fig.show()

fig = mcds.plot_scatter('oxygen', z_slice=3.333)
fig.show()
```

BUE 20240808: DANG.
**make_cell_vtk**


BUE 20240808: SINGLECELL DATA!
#### Transform an MCDS Time Step into an AnnData Object

The [AnnData](https://anndata.readthedocs.io/en/latest/) format is the de facto standard for single cell data analysis in python.\
It is a one-liner to transform a mcds object onto an anndata object.\
This is the gate to a whole new universe.

+ https://scverse.org/

```python
ann = mcds.get_anndata()
print(ann)  # AnnData object with n_obs × n_vars = 1099 × 80
            #     obs: 'ID', 'current_phase', 'cycle_model'
            #     obsm: 'spatial'
```

The output tells us that we have loaded a time step  with 1099 cells (agents) and 80 attributes.
And that we have spatial coordinate annotation (position\_x, position\_y, position\_z, time) of the loaded data.

Whatever you d'like to do with your physicell data, it most probably was already done with single cell wet lab data.
That's being said: PhysiCell data is different scdata than scRNA seq, for example.
scRNA seq data is higher dimensional (e.g. for the human genome, over 20000 genes each time step) than PhysiCell data (tens, maybe hundreds of attributes).
scRNA seq data is always single time step data because the measurement consumes the sample.
PhysiCell data is always time series data, even we look at this moment only at one time step.
This all means, the wet lab bioinformatics will partially try to solve problems (e.g. trajectory inference), that simply are no problems for us and the other way around.

Anyhow, the gate is open.
For the shake of appearances, let's do a cluster analysis on PhysiCell output using scanpy.

```python
import scanpy as sc

# loads only attribute that have not the same value in all cells.
# max absolute scales the attributes into a range between -1 and 1.
ann = mcds.get_anndata(values=2, scale='maxabs')

# principal component analysis
sc.tl.pca(ann)  # process anndata object with the pca tool.
sc.pl.pca(ann)  # plot pca result.
ann.var_names  # list the numerical attributes we have at hand (alternative way: ann.var.index).
ann.obs_keys()  # list the categories attributes we have at hand (alternative way: ann.obs.columns).
sc.pl.pca(ann, color=['current_phase','oxygen'])  # plot the pca results colored by some attributes.
sc.pl.pca(ann, color=list(ann.var_names)+list(ann.obs_keys()))  # gotta catch 'em all!
sc.pl.pca_variance_ratio(ann)  # plot how much of the variation each principal component captures.

# neighborhood graph clustering
sc.pp.neighbors(ann, n_neighbors=15)  # compute the neighborhood graph with the neighbors preprocess step.
sc.tl.leiden(ann, resolution=0.01)  # cluster the neighborhood graph with the leiden tool.
sc.pl.pca(ann, color='leiden')  # plot the pca results colored by leiden clusters.

# umap dimensional reduction embedding
sc.tl.umap(ann)  # process anndata object with the umap tool.
sc.pl.umap(ann, color=['current_phase','oxygen','leiden'])  # plot the umap result colored by some attributes.
```

## Graph data Functions

Since lately, PhysiCell tracks for each cell, if this cell touches other cells.\
This data is stored in two dictionaries of sets of integers which link the corresponding cell IDs.\
We have here a closer look at the neighbor graph because the attached graph is in this particular study not really interesting.

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


## Mesh Data Functions

For data analysis use, the functions realted to the mesh are most probably the least one you have to deal with.
However, lets have a look these functions anyway:\
It is common, but not necessary, that the voxel's width, height, and depth is the same.\
In fact, in this test dataset you will find that this is not the case. \
In the related `PhysiCell_settings.xml` file the voxel was specified as 30[nm] high, 20[nm] wide, and 10[nm] deep.\

We pcdl we can retrieve anytime the voxel spacing value and voxel volume.\
Additionally, we can retrieve the mesh spacing.\
You will notice that voxel and mesh spacing values differ.\
This is because this data set is from a 2D simulation. For that reason, the mesh depth is set to 1.\
In 3D simulation data, voxel and mesh spacing will be the same because voxel and mesh depth is the same.

```python
mcds.get_mesh_spacing()  # [30.0, 20.0, 1]
mcds.get_voxel_spacing() # [30.0, 20.0, 10.0]
mcds.get_voxel_volume()  # 30[nm] * 20[nm] * 10[nm] = 6000[nm\*\*3]
```

Since version 3, for clarity, coordinate labels are distinct for cell position, mesh center, and voxel index:
+ **x,y,z:** stand for cell position coordinates, and are real values.
+ **m,n,p:** stand for mesh center coordinates, and are real values.
+ **i,j,k:** stand for voxel index coordinates, and are unsigned integer values.

For a better understanding, let's fetch start and stop range values for each coordinate type.\
Unlike in the python range function, not only the start value, but also the stop value is inclusive.

```python
mcds.get_voxel_ijk_range()  # [(0, 10), (0, 10), (0, 0)]
mcds.get_mesh_mnp_range()  # [(-15.0, 285.0), (-10.0, 190.0), (0.0, 0.0)]
mcds.get_xyz_range()  # [(-30.0, 300.0), (-20.0, 200.0), (-5.0, 5.0)]
```

The domain benchmarks are:
+ voxel index values stretch from 0 to 10 longitude (i), 0 to 10 latitude (j), the 2D domain is only 1 voxel deep (k).
+ mesh center values stretch from -15[nm] to 285[nm] longitude (m), -10[nm] to 190[nm] latitude (n), and the depth mesh center (p) is at zero.
+ cells can hold a positioned between -30[nm] and 300[nm] longitude (x), -20[nm] and 200[nm] latitude (y), and -5[nm] and 5[nm] depth (z).

For voxel and mesh centers, we can fetch the axis tick lists.\
Each of this function will return a list of 3 numpy arrays, ordered by ijk or mnp, respective.

```python
mcds.get_voxel_ijk_axis()
mcds.get_mesh_mnp_axis()
```

We can even retrieve all mnp mesh center coordinate triplets, and the meshgrids in 2D or 3D shape itself.\
The coordinate triplets are returned as numpy array.\
The meshgrids are returned as a 3 way tensor (2D) or 4 way tensor (3D) numpy array. \
As shown below, these tensors can be subset in the common numpy way.

```python
# all mnp coordinates
mnp_coordinats = mcds.get_mesh_coordinate()  # numpy array with shape (3, 121)
mnp_coordinats.shape  # (3, 121)

# 2D mesh grid
mn_meshgrid = mcds.get_mesh_2D()  # numpy array with shape (2, 11, 11)
mn_meshgrid.shape # (2, 11, 11)
m_meshgrid = mn_meshgrid[0]  # or explicit mn_meshgrid[0,:,:]
n_meshgrid = mn_meshgrid[1]  # or explicit mn_meshgrid[1,:,:]
m_meshgrid.shape  # (11, 11)
n_meshgrid.shape  # (11, 11)

# 3D mesh grid
mnp_meshgrid = mcds.get_mesh()  # numpy array with shape (3, 11, 11, 1)
mnp_meshgrid.shape # (3, 11, 11, 1)
m_meshgrid = mnp_meshgrid[0]  # or explicit mnp_meshgrid[0,:,:,:]
n_meshgrid = mnp_meshgrid[1]  # or explicit mnp_meshgrid[1,:,:,:]
p_meshgrid = mnp_meshgrid[2]  # or explicit mnp_meshgrid[2,:,:,:]
m_meshgrid.shape  # (11, 11, 1)
n_meshgrid.shape  # (11, 11, 1)
p_meshgrid.shape  # (11, 11, 1)
```

Furthermore, there are two helper function.\
One to figure out if a particular xyz coordinate is still in side the mesh, another one to translate a xyz coordinate into ijk voxel indices.\
The translation function will by default checks if the given coordinate is in the mesh.

```python
# is given xyz is in the mesh?
mcds.is_in_mesh(x=0, y=0, z=0)  # True
mcds.is_in_mesh(x=111, y=22, z=-5)  # True
mcds.is_in_mesh(x=111, y=22, z=-5.1)  # False and Warning @ pyMCDS.is_in_mesh : z = -5.1 out of bounds: z-range is (-5.0, 5.0)

# translate xyz into ijk!
mcds.get_voxel_ijk(x=0, y=0, z=0)  # [0,0,0]
mcds.get_voxel_ijk(x=111, y=22, z=-5)  # [4, 2, 0]
mcds.get_voxel_ijk(x=111, y=22, z=-5.1)  # None and Warning @ pyMCDS.is_in_mesh : z = -5.1 out of bounds: z-range is (-5.0, 5.0)
```


### Data Clean Up

After you are done checking out the results, you can uninstall the test datasets and all files stored within its folders.

```
pcdl.uninstall_data()
```
