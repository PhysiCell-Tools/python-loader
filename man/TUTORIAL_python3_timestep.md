# PhysiCell Data Loader Python Tutorial

In this chapter, we will load the pcdl library and use its TimeStep class to load the data snapshot 00000012, from [data\_timeseries\_ 2d](https://github.com/elmbeech/physicelldataloader/tree/master/pcdl/data_timeseries_2d) from the 2D time series test dataset.

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
For speed and less memory usage, it is however possible to only load the essential (output xml and cell mat data), and exclude microenvironment, graph data, PhysiBoss data, and the PhysiCell\_settings.xml cell type ID label mapping.\
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
```
```python3
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
There should be no need to fetch data directly from the `mcds.data` dictionaries.**

We will explore these functions in the upcomming sections.


## Metadata Related Functions

Fetch the data's MultiCellDS version, and the PhysiCell version the data was generated.

```python
mcds.get_multicellds_version()  # will return a string like MultiCellDS_2 or MultiCellDS_0.5
```
```python
mcds.get_physicell_version()  # will return a string like PhysiCell_1.10.4 or BioFVM_1.1.7
```

Fetch simulation time, runtime, and time stamp when the data was processed.

```python
mcds.get_time()   # will return a float value like 720.0
```
```python
mcds.get_runtime()  # will return a float value like 15.596373
```
```python
mcds.get_timestamp()  # will return a sting like 2022-10-19T01:12:01Z
```

Finally, it is possible to retrieve a dictionary that lists all units from all tracked variables, from metadata, mesh, continuum\_variables, and discrete\_cells.

```python
mcds.get_unit_dict()  # will return a dictionary, which maps the variables to the units specified in the setting.xml.
```
```python
# sorted alphabetically
sorted(mcds.get_unit_dict().items())
```


## Microenvironment Data Related Functions (Continuum Variables)

In the loaded dataset, only one substrate, oxygen, was part of the simulation.\
However, let's have a look at the substrate microenviroment related functions.

We can retrieve list of all substrates processed in the simulation,
ordered by substrate ID.

```python
mcds.get_substrate_list()  # ['oxygen']
```

We can retrieve a dictionary that maps substrate indexes to labels.

```python
mcds.get_substrate_dict()  # {'0': 'oxygen'}
```
```python
# substrates sorted alpahbetically
sorted(mcds.get_substrate_dict().values())  # ['oxygen']
```

We can also retrieve a pandas dataframe with the parameters (decay\_rate, diffusion\_coefficient) that were set for each substrate.

```python
df = mcds.get_substrate_df()
df.head()
```

Further, we can retrieve a **pandas dataframe** with voxel ijk coordinate values, mesh center mnp coordinate values, and concentrations values for all substrates.

```python
# dataframe with complete voxel coordinate (ijk), mesh coordinate (mnp), and substrate concentration mapping
df_conc = mcds.get_conc_df()
df_conc.head()
```

Substartes that have in the whole domain the same concentration, carry for that time step no information.
Let's filter for substrates that have at least 2 different values over the whole domain.
```python
df_conc = mcds.get_conc_df(values=2)
df_conc.head()  # oxygen concentration varies over the domain
```

Let's filter for substarte concentrations voxel i2 j1 k0.
```python
df_conc.loc[(df_conc.voxel_i == 2) & (df_conc.voxel_j == 1) & (df_conc.voxel_k == 0), :]
```

Additionally, there is a less often used function to retrieve substrate specific 3D or 2D **meshgrid numpy arrays**.
To get a 2D meshgrids you can slice though any z stack value, the function will always pick the closest mesh center coordinate, the smaller coordinate, if you hit the saddle point between two voxels.
(This function might become deprected in a future version of pcdl.)

```python
# concentration meshgrid for a particular substrate
oxygen_2d = mcds.get_concentration('oxygen', z_slice=0)
oxygen_2d.shape  # (11, 11)
```
```python
# concentration meshgrid for a particular substrate
oxygen_3d = mcds.get_concentration('oxygen')
oxygen_3d.shape  # (11, 11, 1)
```

Additionally, there is a less often used functions to retrieve a **numpy array** of all substrate concentrations at a particular xyz coordinate, ordered alphabetically by substrate name, like the list retrieved by the get\_substrate\_names function.
(This function might become deprected in a future version of pcdl as well.)

```python
# all concentration values at a particular coordinate
mcds.get_concentration_at(x=0, y=0, z=0)  # array([34.4166271])
```
```python
# all concentration values at a particular coordinate
mcds.get_concentration_at(x=111, y=22, z=-5)  # array([18.80652216])
```
```python
# all concentration values at a particular coordinate
mcds.get_concentration_at(x=111, y=22, z=-5.1)  # None and Warning @ pyMCDS.is_in_mesh : z = -5.1 out of bounds: z-range is (-5.0, 5.0)
```

### &#x2728; Microenvironment Data Analysis with [Pandas](https://pandas.pydata.org/)

Pandas is the library you will probably spend a lot of the time when you are analysis you data.
Pandas mimics the computer language R.
Pandas provides us with the Series and Dataframe data types.

For getting started with pandas,
and get a clearer idea what pandas is all about and capable of,
I recommend you to work through this pandas cookbook from Julia Evens:
+ https://jvns.ca/blog/2013/12/22/cooking-with-pandas/

The Pandas library as such is very well documented.
Please spend some time to familarize your self with the homepage, the users guides, and the API reference.
+ https://pandas.pydata.org/
+ http://pandas.pydata.org/pandas-docs/stable/user_guide/index.html
+ http://pandas.pydata.org/pandas-docs/stable/reference/index.html


### &#x2728; Microenvironment Data Analysis with [Mathplotlib](https://matplotlib.org/)

For substrate concentration visualization **matplotlib contour and contourf plots**,
for any substrate, through any z\_slice can be retrived.\
The mcds.plot\_contour function has many parameters to fine tune the plot.
Please have a look at it's docstring to learn more.\
The mcds.plot\_contour function output can combine with the mcds.plot\_scatter output.

+ https://matplotlib.org/

Please have a look at the [TUTORIAL_python3_matplotlib.md]() to learn more.

```python
# contour plot
fig = mcds.plot_contour('oxygen')
fig = mcds.plot_contour('oxygen', z_slice=3.333)
fig.show()
```
```python
help(mcds.plot_contour)
```


### &#x2728; Microenvironment Data Analysis with [Vtk](https://vtk.org/)

For substrate concentration visualization **rectilinear grid vtk files**
for any substrate can be retrived.\
This files can be analysied, for example with the [Paraview](https://en.wikipedia.org/wiki/ParaView) software.

+ https://vtk.org/

Please have a look at the [TUTORIAL_paraview.md]() to learn more.

```python
mcds.make_conc_vtk()
```

## Cell Data Related Functions (Discrete Cells)

In the loaded dataset, only one cell type, cancer\_cell, was part of the simulation.\
However, let's have a look at the cell agenet related functions.

We can retrieve list of all cell types processed in the simulation,
ordered by cell type ID.

```python
mcds.get_celltype_list()  # ['cancer_cell']
```

We can retrieve a dictionary that maps cell type IDs to labels.

```python
mcds.get_celltype_dict()  # {'0': 'cancer_cell'}
```

Now, let's have a look at all variables we track from each agent.\
This data we can retrieve as dataframe. \
By the way, this dataframes not only hold the exact xyz coordinate and all variables.
The datafram list additionally ijk voxel coordinate, mnp mesh center coordinate,
the cell surrounding substrate concentrations, and, if available, PhysiBosss output.

```python
# data from all agents in the domain
df_cell = mcds.get_cell_df()
df_cell.shape  # (992, 95)  this means: 992 agents in the whole domain, 96 tracked variables
```
```python
df_cell.info()
```
```python
df_cell.head()
```
```python
sorted(df_cell.columns)
```

Cell attributes that carry in all agents the same value, carry no information.
Let's filter for variables that carry at least 2 different values.

```python
df_cell = mcds.get_cell_df(values=2)
df_cell.shape # (992, 40) this means: 992 agents in the whole domain, 40 tracked variables have more than 2 different values
```

Let's filter for cells in voxel i2 j1 k0.

```python
df_cell.loc[(df_cell.voxel_i == 2) & (df_cell.voxel_j == 1) & (df_cell.voxel_k == 0), :]  # cells: 5, 7, 39
```

There exist an additional, less often used function,
to filter for cells in xyz position plus minus (voxel spacing / 2).
(This function might become deprected in a future version of pcdl.)

```python
mcds.get_cell_df_at(x=45, y=10, z=0)  # cells: 5, 7, 39
```

### &#x2728; Cell Data Analysis with [Pandas](https://pandas.pydata.org/)


BUE PANDAS!
+ https://pandas.pydata.org/


### &#x2728; Cell Data Analysis with in the [Scverse](https://scverse.org/)

To be able to analyse cell agent data the same way single cell RNA seq data is analysed,
pcdl has a function to translate cell agent data into [AnnData](https://anndata.readthedocs.io/en/latest/) format.
Anndata is the de facto standard for sc RNA seq  data analysis in python.
Anndata is the backbone of the sc verse project.

+ https://scverse.org/

````python
ann = mcds.get_anndata(values=2)
print(ann)  # AnnData object with n_obs × n_vars = 992 × 26
            #     obs: 'z_layer', 'time', 'current_phase', 'cycle_model'
            #     uns: 'neighbor'
            #     obsm: 'spatial'
            #     obsp: 'physicell_neighbor_conectivities', 'physicell_neighbor_distances'
```

variables
```python
ann.var_name  # numerical cell attributes:  Index(['cell_BM_repulsion_strength', ... , 'total_volume'], dtype='object')
```

observation

```python
ann.obs_names  # cell IDs: Index(['0', ..., '994'], dtype='object', name='ID', length=992)
```
```python
ann.obs_keys()  # categorical cell attributes:  ['z_layer', 'time', 'current_phase', 'cycle_model']
```
```python
ann.obsm_keys()  # cell coordinates: ['spatial']
```
```python
ann.obsm['spatial']  # the coordinate values (position\_x, position\_y) of the loaded data.
```
```python
ann.obsp  # cell neigborhood graph data.
```

unstructured data

```python
ann.uns_keys()  # ['neighbor']
```
```python
ann.uns['neighbor']  # metadata about the neighborhood graph.
```

The output tells us that we have loaded a time step  with 992 cell agents and 26 numerical attributes (vars).
Further, we have 4 categorical cell agent attributes (obs).
We have each cell agnet's spatial coordinate information (obsm).
And we have cell neighbor graph infromation (obsp, uns).

Whatever you d'like to do with your physicell data, it most probably was already done with single cell wet lab data.
That's being said: PhysiCell data is different scdata than scRNA seq!
For example, scRNA seq data is higher dimensional (e.g. the human genome has over 20000 genes each time step) than PhysiCell data (tens, maybe hundreds of cell attributes).
For example, scRNA seq data is always single time step data because the measurement consumes the sample. PhysiCell data is always time series data, even we look at this moment only at one time step.
This means, the wet lab bioinformatics will partially try to solve problems (for example trajectory inference), that simply are no problems for us and the other way around.
Anyhow, there are a lot of scRNA seq data analysis methodes around, whoch make sense to apply to both of this data types.

Please have a look at the [TUTORIAL_python3_scverse.md]() to learn more.


### &#x2728; Cell Data Analysis with in the [Networkx](https://networkx.org/) and [Igraph](https://igraph.org/)

+ https://networkx.org/
+ https://igraph.org/

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

### &#x2728; Cell Data Analysis with [Mathplotlib](https://matplotlib.org/)

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


## Microenvironment and Cell Data Related Functions

### &#x2728; PhysiCell Data Analysis with [Napari](https://napari.org/stable/) and [Fiji Imagej](https://fiji.sc/)

The open microscopy's [ome.tiff](https://www.openmicroscopy.org/ome-files/) file fromat.

+ https://napari.org/stable/
+ https://fiji.sc/


## Mesh Data Related Functions

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
