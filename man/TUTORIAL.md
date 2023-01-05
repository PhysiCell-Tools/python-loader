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
All files from the first round of output will end in 00000000.* and the second round will be 00000001.* and so on.\
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
+ **output00000012_cells.mat**: This is a MATLAB matrix file that contains all tracked information about the individual cells in the model.
    It tells us things like the cells’ position, volume, secretion, cell cycle status, and user-defined cell parameters.
+ **output00000012_microenvironment0.mat**: This is a MATLAB matrix file that contains all data about the microenvironment at this time step.


### Loading an MCDS into Python3
In this paragraph we will load the pcDataLoder library and use its pyMCDS class to load the data snapshot 00000012, described above, from [data timeseries 2d](https://github.com/elmbeech/pcDataLoader/tree/v3/pcDataLoader/data_timeseries_2d) from the test dataset.
(There is no need to extra install the test data set.
In fact, both test datasets are already installed.
They ship with the pip3 pcDataLoader installation.)

```python
import pathlib  # library to locate the test data
import pcDataLoader as pc  # the PhysiCell data loader library

s_path = str(pathlib.Path(pc.__file__).parent.joinpath('data_timeseries_2d'))  # local path to the installed test data set
s_file = 'output00000012.xml'  # the snapshot we want to analyse
s_pathfile = f'{s_path}/{s_file}'
print('mcds xml:', s_pathfile)

# load mcds - multi cell digital snapshot - object
mcds = pc.pyMCDS(s_pathfile)  # loads whole snapshot: the xml and all realated mat and graph files.
```

Side note: in general, unix and windows path notation will work.\
The legacy way of data loaing works too.
```python
# legacy way of loading a mcds object
mcds = pc.pyMCDS('output00000012.xml', s_path)
```

By default, all data realted to the snapshot is loaded.\
For speed and less memory usage, it is however possible to only load the essential (xml and cell mat data), and exclude microenviroment and graph data loading.
```python
# fine tuned way of loading a mcds object
mcds = pc.pyMCDS(s_pathfile, graph=False, microenv=False)
```

This is the easiest way to check, which version of pcDataLoader you run.
```python
pc.__version__
```

### The Loaded Data Structure

All loaded data lievs in `mcds.data` dictionary.\
As in **version 1**, we’ve tried to keep everything organized inside this dictionary.\
Regarding **version 1**, the structure has slightly changed.\
However, in **version 3**, all data is accessible by functions, thus there should be no need to fetch data directely form the `mcds.data` dictionary!\
Anyhow, let’s take a look at what we actually have in here.

![mcds.data dictionary blueprint](img/pcdataloader_data_dictionary_v3.0.0.png)

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
mcds.get_physicell_version()  # will return a string like PhysiCell_1.10.4 or BioFVM_1.1.7
```

Fetch simulation time, runtime, and time stamp when the data was processed.
```python
mcds.pyMCDS.get_time()   # will retun a float value like 720.0
mcds.pyMCDS.get_runtime()  # will retun a float value like 15.596373
mcds.get_timestamp()  # will retun a sting like 2022-10-19T01:12:01Z
```

#### Mesh Data

Let's start by the voxel.\
It is common, but not necessary, that the voxel's width, heigth, and depth is the same.\
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

For voxel and mesh centers we can fetch the axis tick lists.
Each of this function will return a list of 3 numpy arrays, orderd ijk or mnp, respective.
```python
mcds.get_voxel_ijk_axis()
mcds.get_mesh_mnp_axis()
```

We can even retive all mnp mesh center coordinate triplets, and the meshgrids in 2D or 3D form itself.
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

We can retrieve a pandas dataframe with the parameters (decay\_rate, diffusion\_coefficient) that were set for in each substrate.
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

Now, let's have a look at all variables we track from each agent.\
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

An exciting thing about modeling is to have time series data.\
pcDataLoader's pyMCDSts class is here to make handle of a time series of MCD snapshots easy.\
```python
import pathlib  # library to locate the test data
import pcDataLoader as pc  # the PhysiCell data loader library

s_path = str(pathlib.Path(pc.__file__).parent.joinpath('data_timeseries_2d'))  # local path to the installed test data set
print('mcds time series:', s_path)

# load mcds - multi cell digital snapshot - object
mcdsts = pc.pyMCDSts(s_path)  
```

Like in the pyMCDs class, for memory consumption and processing speed control, we can specify if we later on want to load microenviroment data and graph data from the snapshots we later on analyse.\
By default all data realted to the snapshot is loaded.\
```python
# fine tuned way of loading a mcds object
mcdsts = pc.pyMCDSts(s_pathfile, graph=False, microenv=False)
```

#### Times Series MCDS Data

Now, let's have a look at the pyMCDSts instances data analysis functions.\

+ **get_xmlfile_list** will simply return a list of absolute path of the mcds xml files.\
This list can be maipulated a later be fed into the objetcs read\_xml function to only read sudden snapshots into memory.\
For example, you want only to analys the last hour 11, 12, 13 from your run, or data form every other hour will be enough.
```python
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
The default setting will read all snapshots from a timeseris, however, you can alway use a filtered ls\_xml list to only load a subset of snapshots.
```python
# load all snapshots
l_mcds = mcdsts.read_mcds()
l_mcds  # YMMV! [<pcDataLoader.pyMCDS.pyMCDS at 0x7fa660996b00>, <pcDataLoader.pyMCDS.pyMCDS at 0x7fa67673e1d0>, ..., <pcDataLoader.pyMCDS.pyMCDS at 0x7fa660950f10>]
len(l_mcds)  # 25

# load snapshot 11, 12, and 13
l_mcds_11_12_13 = mcdsts.read_mcds(ls_xml_11_12_13)
len(l_mcds_11_12_13)  # 3

# load all even snapshots
ls_xml_even = mcdsts.read_mcds(ls_xml_even)
len(ls_xml_even)  # 13
```

Single now snapshots can be accessed by indexing.\
With a single snapshot you work in exactely same way as with an object loaded by pyMCDS.\
```python
# get the simmulation time
l_mcds[12].get_time()  # 720.0 

# get the cell data frame
df = l_mcds[12].get_cell_df()
df.shape  # (992, 87)
df.head()
```

It is very easy to loop over the whole time series to analyse timeseris data and generate time series plots.
```python
# load library
import pandas as pd

# fetch data
lr_time = [mcds.get_time() for mcds in l_mcds]  # [0.0, 60.0, ..., 1440.0]
li_cellcount = li_cellcount = [mcds.get_cell_df().shape[0] for mcds in l_mcds]  # [889, 898, ..., 1099]

# plot data
df = pd.DataFrame([lr_time,li_cellcount], index=['time_min','cell_count']).T
df.head()
df.plot(kind='scatter', x='time_min', y='cell_count', grid=True)
```

#### Times Series Scatterplot Images and Movies

With PhysiCell is not only possible to take data snapshots, but as well [svg](https://en.wikipedia.org/wiki/SVG) vector graphics images snapshots.\
PhysiCell's [Makefile](https://en.wikipedia.org/wiki/Make_(software)) has code to translate those svg images into [gif](https://en.wikipedia.org/wiki/GIF), [jpeg](https://en.wikipedia.org/wiki/JPEG), [png](https://en.wikipedia.org/wiki/Portable_Network_Graphics), or [tiff](https://en.wikipedia.org/wiki/TIFF) fomat, making use of the [image magick](https://en.wikipedia.org/wiki/ImageMagick) library, and to translate the jpeg, png, or tiff images into a [mp4](https://en.wikipedia.org/wiki/MP4_file_format) movie, therefore making use from the [ffmpeg](https://en.wikipedia.org/wiki/FFmpeg) library.\
pyMCDSts instances provides the similar functionallity.\
This means this following code will only run if image magick and ffmpeg are installed on your operating system.\
The svg images might be quite huge. You can always use the `resize_factor` parameter to scale donw the image size for the resulting images. Resizing will lead to processing time speed up and saves disk space.\

Translate physicells svg images into static raster graphic images:
```python
# resize factor 1 will leave the image size as it is.
mcds.make_jpeg()  # resize factor 1 
mcds.make_png()  # resize factor 1
mcds.make_tiff()  # esize factor 1

# resize factor 0.2 will down scale to 20% width and length of the image.
mcds.make_jpeg(0.2) 
mcds.make_png(0.2)
mcds.make_tiff(0.2)
```

Translate physicells svg images into a dynamic gif image:\
The default file name for the resulting gif image is timeseries.gif.
```python
# resize factor 1 
make_gif()

# resize factor 0.2
make_gif(0.2)
```

Translate physicells svg images into a mp4 movie:\
Movies can only be generated for already existing jpeg, png, or tiff moves!\
By default jpegs files will be used to generat the movie.\
If png or tiff files should be used as source, then this have to be explicite stated.\
The default file name from the resulting movie is movie.mp4.
```python
make_movie()  # generate move from jpeg files
make_movie('jpeg')  # generate move from jpeg files
make_movie('png')  # generate move from png files
make_movie('tiff')  # generate move from tiff files
```

**That's all Folks!**
