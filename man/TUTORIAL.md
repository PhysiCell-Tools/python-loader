# PhysiCell Data Loader Tutorial Man Page

Please install the latest version of physicelldataloader (pcdl) branch v3,
as described in the [HowTo](https://github.com/elmbeech/physicelldataloader/blob/master/man/HOWTO.md) section.\
Branch v1 and v2 exists mainly for backwards compatibility.


## Tutorial - branch v1 and v2
The original python-loader tutorial can be found here.
+ http://www.mathcancer.org/blog/python-loader/


## Tutorial - branch v3

### Introduction: PhysiCell's MultiCellDigital Snapshot (MCDS) Anatomy

Each time PhysiCell's internal time tracker passes a time step where data is to be saved, it generates a number of files of various types.\
Each of these files will have a number at the end that indicates where it belongs in the sequence of outputs.\
All files from the first round of output will end in 00000000.\*, and the second round will be 00000001.\*, and so on.\
Have a look at this [PhysiCell data time series](https://github.com/elmbeech/physicelldataloader/tree/master/pcdl/data_timeseries_2d).

Let's assume we captured data every simulation time hour, and we're interested in the set of output half a day through the run, the 13th set of output files.\
The files we care about most from this set consists of:

+ **output00000012.xml**: This file is the main organizer of the data.
    It contains an overview of the data stored in the MultiCellDS as well as some actual data, including:\
    metadata (MultiCellDS version, PhysiCell or BioFVM version, simulation time, runtime, and processing time stamp),\
    coordinates for the computational domain (mesh),\
    parameters for diffusing substrates in the microenvironment (continuum\_variables),\
    column labels and units for the cell data (cell\_population),\
    file names for the files that contain microenvironment and cell data at this time step (mat and possibly graph.txt files),
+ **output00000012_cells.mat**: This is a MATLAB matrix file that contains tracked information about the individual cells in the model.
    It tells us things like the cells' position, volume, secretion, cell cycle status, and user defined cell parameters.
+ **output00000012_microenvironment0.mat**: This is a MATLAB matrix file that contains data about the microenvironment at this time step.


### Loading an MCDS into Python3
In this section, we will load the pcDataLoder library, and use its pyMCDS class to load the data snapshot 00000012, described above, from [data\_timeseries\_ 2d](https://github.com/elmbeech/physicelldataloader/tree/master/pcdl/data_timeseries_2d) from the 2D timeseries test dataset.

```python
import pathlib  # library to locate the test data
import pcdl  # the physicell data loader library

print('pcdl version:', pcdl.__version__)  # it is easy to figure out which pcdl version you run

pcdl.install_data()  # to install the test datasets

s_path = str(pathlib.Path(pcdl.__file__).parent.joinpath('data_timeseries_2d'))  # local path to the installed test data set
s_file = 'output00000012.xml'  # the snapshot we want to analyze
s_pathfile = f'{s_path}/{s_file}'
print('mcds xml:', s_pathfile)

# load mcds - multi cell digital snapshot - object
mcds = pcdl.pyMCDS(s_pathfile)  # loads the whole snapshot: the xml and all related mat and graph files.
```

Side note: in general, unix and windows path notation will work.\
The legacy way of loading data works too.
```python
# legacy way of loading a mcds object
mcds = pcdl.pyMCDS('output00000012.xml', s_path)
```

By default, all data related to the snapshot is loaded.\
For speed and less memory usage, it is however possible to only load the essential (output xml and cell mat data),
and exclude microenvironment, graph data, and ID label mapping loading.\
Additionally, it is possible to specify for custom_data variable types other than the generic float type, namely: int, bool, and str.
```python
# fine tuned way of loading a mcds object
mcds = pcdl.pyMCDS(s_pathfile, custom_type={}, microenv=False, graph=False, settingxml=False)
```

### The Loaded Data Structure

All loaded data lives in `mcds.data` dictionary.\
As in **version 1**, we have tried to keep everything organized inside this dictionary.\
Regarding **version 1**, the structure has slightly changed.\
However, in **version 3**, all data is accessible by functions. There should be no need to fetch data directly from the `mcds.data` dictionary!\
Anyhow, let's take a look at what we actually have in here.

![mcds.data dictionary blueprint](img/physicelldataloader_data_dictionary_v3.2.2.png)

```python
# main data branches
sorted(mcds.data.keys())  # metadata, mesh, substrate (continuum_variables), and agent (discrete_cells).

# metadata
sorted(mcds.data['metadata'].keys())  # multicellds version, physicell version, simulation time, runtime, time stamp, time unit, spatial unit, and substrate and cell type ID label mappings.

# mesh
sorted(mcds.data['mesh'].keys())  # voxel (ijk), mesh (nmp), and position (xyz) range, axis, coordinate, grid objects, and voxel volume.

# microenvironment
sorted(mcds.data['continuum_variables'].keys())  # list of all processed substrates, e.g. oxygen.
sorted(mcds.data['continuum_variables']['oxygen'].keys())  # substrate related data values, unit, diffusion coefficient, and decay rate.

# cell
sorted(mcds.data['discrete_cells'].keys())  # data, units, and graph dictionaries.
sorted(mcds.data['discrete_cells']['data'].keys())  # all cell related, tracked data.
sorted(mcds.data['discrete_cells']['units'].keys())  # all units from the cell related, tracked data.
sorted(mcds.data['discrete_cells']['graph'].keys())  # neighbor_cells and attached_cells graph dictionaries.
```

### Accessing the Loaded Data by the pyMCDS Class Functions
Once again, loud, for the ones in the back, in version 3, all data is accessible by functions.\
There should be no need to fetch data directly from the `mcds.data` dictionary of dictionaries.\
We will explore these functions in this section.

#### Metadata

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

Fetch substrate and cell type ID label mappings, read out from the PhysiCell_settings.xml file.
```python
mcds.get_substarte_dict()  # will return a dictionary, which maps substrate IDs to labels.
mcds.get_celltype_dict()   # will return a dictionary, which maps cell type IDs to labels.
```

#### Mesh Data

Let's start by the voxel.\
It is common, but not necessary, that the voxel's width, height, and depth is the same.\
However, you will find that in this test dataset this is not the case. \
In the related `PhysiCell_settings.xml` file the voxel was specified as 30[nm] high, 20[nm] wide, and 10[nm] deep.\
We can retrieve this voxel spacing value and voxel volume. \
Additionally, we can retrieve the mesh spacing.\
You will notice that voxel and mesh spacing values differ.\
This is because this data set is from a 2D simulation. For that reason, the mesh depth is set to 1.\
In 3D simulation data, voxel and mesh spacing will be the same because voxel and mesh depth is the same.
```python
mcds.get_mesh_spacing()  # [30.0, 20.0, 1]
mcds.get_voxel_spacing() # [30.0, 20.0, 10.0]
mcds.get_voxel_volume()  # 30[nm] * 20[nm] * 10[nm] = 6000[nm\*\*3]
```

Since version 3, coordinate labels are distinct:
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
As shown below, it is easy to subset these tensors.
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
mcds.is_in_mesh(x=111, y=22, z=-5.1)  # False and Warning @ pyMCDS.is_in_mesh : z = -5.1 out of bounds: z-range is (-5.0, 5.0).

# translate xyz into ijk!
mcds.get_voxel_ijk(x=0, y=0, z=0)  # [0,0,0]
mcds.get_voxel_ijk(x=111, y=22, z=-5)  # [4, 2, 0]
mcds.get_voxel_ijk(x=111, y=22, z=-5.1)  # None and Warning @ pyMCDS.is_in_mesh : z = -5.1 out of bounds: z-range is (-5.0, 5.0).
```

#### Microenvironment (Continuum Variables) Data

Let's have a look at the substrate.\
We can retrieve an alphabetically ordered list of all substrates processes in the simulation.\
In the loaded dataset, only one substrate, oxygen, was part of the simulation.
```python
mcds.get_substrate_names()  # ['oxygen']
```

We can retrieve a pandas dataframe with the parameters (decay\_rate, diffusion\_coefficient) that were set for each substrate.
```python
df = mcds.get_substrate_df()
df.head()
```

Regarding the concentrations, we can retrieve:
+ a **numpy array** of all substrate concentrations at a particular xyz coordinate, ordered alphabetically by substrate name, like the list retrieved by the get\_substrate\_names function.
+ substrate specific 3D or 2D **meshgrid numpy arrays**.
  To get a 2D meshgrids you can slice though any z stack value, the function will always pick the closest mesh center coordinate, the smaller coordinate, if you hit the saddle point between two voxels.
+ **matplotlib contour and contourf plots**, for any substrate, through any z\_slice.
+ a **pandas dataframe** with voxel ijk coordinate values, mesh center mnp coordinate values, and concentrations values for all substrates.
```python
# all concentration values at a particular coordinate
mcds.get_concentration_at(x=0, y=0, z=0)  # array([34.4166271])
mcds.get_concentration_at(x=111, y=22, z=-5)  # array([18.80652216])
mcds.get_concentration_at(x=111, y=22, z=-5.1)  # None and Warning @ pyMCDS.is_in_mesh : z = -5.1 out of bounds: z-range is (-5.0, 5.0).

# concentration meshgrid for a particular substrate
oxygen_3d = mcds.get_concentration('oxygen')
oxygen_2d = mcds.get_concentration('oxygen', z_slice=0)
oxygen_3d.shape  # (11, 11, 1)
oxygen_2d.shape  # (11, 11)

# contour plot
fig = mcds.get_contour('oxygen')
fig = mcds.get_contour('oxygen', z_slice=3.333)
fig.show()

# dataframe with complete voxel coordinate (ijk), mesh coordinate (mnp), and substrate concentration mapping.
df = mcds.get_concentration_df()
df.head()
```

#### Cell (Discrete Cells) Data

Now, let's have a look at all variables we track from each agent.\
It is quite a list.\
This list is ordered alphabetically, except for the first entry, which always is ID.
```python
mcds.get_cell_variables()  # ['ID', 'attachment_elastic_constant', ..., 'velocity_z']
len(mcds.get_cell_variables())  # 77
```

The data itself we can retrieve as dataframe.\
There exist two functions to do so.\
One that retrieves data from all agents in the whole domain, the other retrieves only data from the agents in a particular voxel, specified by the xyz coordinate.\
Please note, this dataframes not only hold the exact xyz coordinate, and all discrete variables mentioned above,
they list additionally all agent surrounding substrate concentrations, the containing ijk voxel coordinate, and the mnp mesh center coordinate.

```python
# data from all agents in the domain.
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

#### Unites

Finally, it is possible to retrieve a dataframe that lists all units from all tracked variables, from metadata, mesh, continuum\_variables, and discrete\_cells.
```python
df = mcds.get_unit_df()
df.shape  # (82, 1)
df.columns  # Index(['unit'], dtype='object')
df.index  # Index(['attachment_elastic_constant', 'attachment_rate', ..., 'velocity_z'], dtype='object', name='parameter')
df.head()
```

### MCDS Plotting

Since microenvironment data and cell data can be retrieved as pandas datafarme,
basic plotting
(line plot, bar plot, histogram, boxplot, kernel density estimation plot, area plot, pie plot, scatter plot, hexbin plot)
can easily be generated with the **pandas plot function**.
As mentioned above, for microenvironment data the pyMCDS class has a mcds.get\_contour function because pandas has no contour and contourf plots implementation.
All these plots are mathplotlib plots, hence fine tuning can always be done using the matplotlib library.

+ https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.plot.html
+ https://matplotlib.org/

```python
# load library
import matplotlib.pyplot as plt

# set constantes
z_slice = 0

# generate plot
fig, ax = plt.subplots(figsize=(7,4))

# plot microenvironment
mcds.get_contour('oxygen', z_slice=z_slice, ax=ax)

# plot cells
df = mcds.get_cell_df()
df = df.loc[(df.mesh_center_p == z_slice),:]
df.plot(
    kind = 'scatter',
    x = 'position_x',
    y = 'position_y',
    s = 10,
    c = 'oncoprotein',
    cmap = 'magma',
    vmin = 0,
    vmax = 2,
    grid = True,
    title = 'cells and microenvironment',
    ax = ax,
)
ax.axis('equal')
plt.tight_layout()

# save plot to file
# bue 20230101: note, in this test dataset, cells were seeded outside the actual domain.
fig.savefig('pymcds_2d_cell_and_microenvironment.png', facecolor='white')
plt.close()
```


### Working With MCDS Time Series in Python3

An exciting thing about modeling is to have time series data.\
Pcdl's pyMCDSts class is here to make the handling of a time series of MCD snapshots easy.
```python
import pathlib  # library to locate the test data
import pcdl  # the physicell data loader library

s_path = str(pathlib.Path(pcdl.__file__).parent.joinpath('data_timeseries_2d'))  # local path to the installed test data set
print('mcds time series:', s_path)

# generate a mcds time series instance
mcdsts = pcdl.pyMCDSts(s_path)
```

Like in the pyMCDs class, for memory consumption and processing speed control, we can specify if we want to load microenvironment data and graph data from the snapshots we later on analyze.\
By default, all data related to the snapshot is loaded.\
Additionally, we can exclude to read the PhysiCell_settings.xml file, if it is not available, and fine tune variable typing.
```python
# fine tuned the way mcds objects will be loaded
mcdsts = pcdl.pyMCDSts(s_path, custom_type={}, microenv=False, graph=False, settingxml=False)
```

#### Times Series MCDS Data

Now, let's have a look at the pyMCDSts instances data analysis functions.\

+ **get_xmlfile_list** will simply return a list of absolute path and mcds xml file name strings.\
This list can be manipulated, and later be fed into the objects read\_xml function, to only read specific snapshots into memory.
For example, if you want only to analyze hour 11, 12, and 13 from your run, or if data from every other hour will be enough.
```python
mcdsts = pcdl.pyMCDSts(s_path)
ls_xml = mcdsts.get_xmlfile_list()
ls_xml   # ['/path/to/output00000000.xml', '/path/to/output00000001.xml', ..., '/path/to/output00000024.xml']

# filter for snapshot 11,12, and 13
ls_xml_11_12_13 = ls_xml[11:14]
ls_xml_11_12_13  # ['/path/to/output00000011.xml', '/path/to/output00000012.xml', '/path/to/output00000013.xml']

# filter for every other snapshot
ls_xml_even = [s_xml for i, s_xml in enumerate(ls_xml) if (i%2 == 0)]
ls_xml_even  # ['/path/to/output00000000.xml', '/path/to/output00000002.xml', ..., /path/to/output00000024.xml']
```

+ **read_xml** will actually read the snapshots into RAM.\
The default setting will read all snapshots from a time series. \
However, you can always use a filtered ls\_xml list to only load a subset of snapshots.
```python
# load all snapshots
mcdsts.read_mcds()
mcdsts.l_mcds  # YMMV! [<pcdl.pyMCDS.pyMCDS at 0x7fa660996b00>, <pcdl.pyMCDS.pyMCDS at 0x7fa67673e1d0>, ..., <pcdl.pyMCDS.pyMCDS at 0x7fa660950f10>]
len(mcdsts.l_mcds)  # 25

# load snapshot 11, 12, and 13
mcdsts.l_mcds = mcdsts.read_mcds(ls_xml_11_12_13)
len(mcdsts.l_mcds)  # 3

# load all even snapshots
mcdsts.ls_mcds = mcdsts.read_mcds(ls_xml_even)
len(mcdsts.ls_mcds)  # 13
```

Single snapshots can now be accessed by indexing.\
With a single snapshot, you work exactly in the same way as with an object loaded by pyMCDS.
```python
# get the simulation time
mcdsts.l_mcds[12].get_time()  # 720.0

# get the cell data frame
df = mcdsts.l_mcds[12].get_cell_df()
df.shape  # (992, 87)
df.head()
```

It is very easy to loop over the whole time series, to analyze time series data, and generate time series plots.
```python
# load library
import matplotlib.pyplot as plt
import pandas as pd

# fetch data
mcdsts = pcdl.pyMCDSts(s_path)
mcdsts.read_mcds()
lr_time = [mcds.get_time() for mcds in mcdsts.l_mcds]  # [0.0, 60.0, ..., 1440.0]
li_cellcount = [mcds.get_cell_df().shape[0] for mcds in mcdsts.l_mcds]  # [889, 898, ..., 1099]

# pack data into a pandas datafarm
df = pd.DataFrame([lr_time,li_cellcount], index=['time_min','cell_count']).T
df.head()

# plot data
df.plot(
    kind = 'scatter',
    x = 'time_min',
    y = 'cell_count',
    s = 36,
    ylim = (800,1200),
    grid=True,
    title='2D time series test data'
)

# save plot to file
fig, ax = plt.subplots()
df.plot(
    kind = 'scatter',
    x = 'time_min',
    y = 'cell_count',
    s = 36,
    ylim = (800,1200),
    grid = True,
    title = '2D time series test data',
    ax = ax,
)
plt.tight_layout()
fig.savefig('pymcdsts_2d_cellcount.png', facecolor='white')
plt.close()
```

#### Times Series Scatter Plot Images, Contour Plot Images, and Movies

With PhysiCell it is not only possible to take data snapshots, but as well [svg](https://en.wikipedia.org/wiki/SVG) vector graphics images snapshots.\
PhysiCell's [Makefile](https://en.wikipedia.org/wiki/Make_(software)) has code to translate those svg images into [gif](https://en.wikipedia.org/wiki/GIF), [jpeg](https://en.wikipedia.org/wiki/JPEG), [png](https://en.wikipedia.org/wiki/Portable_Network_Graphics), or [tiff](https://en.wikipedia.org/wiki/TIFF) format, making use of the [image magick](https://en.wikipedia.org/wiki/ImageMagick) library, and to translate the jpeg, png, or tiff images into a [mp4](https://en.wikipedia.org/wiki/MP4_file_format) movie, therefore making use from the [ffmpeg](https://en.wikipedia.org/wiki/FFmpeg) library.\
pyMCDSts instances provides a similar functionality, although the images are generated straight from data and not for the svg files.\
However, mp4 movies and gif images are generated in the same way.\
This means the mcdsts.make_gif code will only run if image magick and mcdsts.make_movie code will only run if ffmpeg is installed on your computer.

```python
# fetch data
mcdsts = pcdl.pyMCDSts(s_path)
```

Translate physicell output data into static raster graphic images:

```python
## cell data ##

# cell images colored by cell_type
mcdsts.make_imgcell()  # jpeg images colored by cell_type
mcdsts.make_imgcell(ext='png')  # png images colored by cell_type
mcdsts.make_imgcell(ext='tiff')  # tiff images colored by cell_type

# cell images can be colored by any cell dataframe column,
# like, for example, by pressure.
mcdsts.make_imgcell(focus='pressure')  # jpeg images colored by pressure values

# there are parameters to adjust cell size and more.


## substrate data ##

# there is an equivalent function for substrate data,
# being able to display any substrate dataframe column.
mcdsts.make_imgsubs(focus='oxygen')  # jpeg images colored by oxygen values
```

Translate raster graphic images into a dynamic gif image:\
Gif images can only be generated for already existing jpeg, png, or tiff images!\
By default, jpeg files will be used, to generate the movie.\
If png or tiff files should be used as source, then this has to be explicitly stated.

```python
# using jpeg images for input
s_path = mcdsts.make_imgcell()
mcdsts.make_gif(s_path)

# using jpeg images one-liner
mcdsts.make_gif(mcdsts.make_imgcell())

# using tiff images for input
s_path = mcdsts.make_imgcell(ext='tiff')
mcdsts.make_gif(s_path, interface='tiff')
```

Translate physicell svg images into an mp4 movie:\
Movies can only be generated for already existing jpeg, png, or tiff images!

```python
# using jpeg images for input
s_path = mcdsts.make_imgcell()
mcdsts.make_movie(s_path)

# using jpeg images one-liner
mcdsts.make_movie(mcdsts.make_imgcell())

# using tiff images for input
s_path = mcdsts.make_imgcell(ext='tiff')
mcdsts.make_movie(s_path, interface='tiff')
```


### Data Clean Up ###
After you are done checking out the results,
you can uninstall the test datasets and all files stored within its folders.
```
pcdl.uninstall_data()
```


**That's all Folks!**
