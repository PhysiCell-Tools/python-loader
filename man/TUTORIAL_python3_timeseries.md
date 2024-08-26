# PhysiCell Data Loader Python Tutorial

Please install the latest version of physicelldataloader (pcdl),
as described in the [HowTo](https://github.com/elmbeech/physicelldataloader/blob/master/man/HOWTO.md) section.

And have a quick read of the pcdl [background](https://github.com/elmbeech/physicelldataloader/tree/master/man/TUTORIAL_introduction.md) infromation.


## Loading an MCDS (initialize)

In this section, we will load the pcdl library and use its TimeStep class to load the data snapshot 00000012, described above, from [data\_timeseries\_ 2d](https://github.com/elmbeech/physicelldataloader/tree/master/pcdl/data_timeseries_2d) from the 2D time series test dataset.

```python
import pathlib  # library to locate the test data
import pcdl  # the physicell data loader library

print('pcdl version:', pcdl.__version__)  # it is easy to figure out which pcdl version you run

pcdl.install_data()  # to install the test datasets

s_path = str(pathlib.Path(pcdl.__file__).parent.joinpath('data_timeseries_2d'))  # local path to the installed test data set
s_file = 'output00000012.xml'  # the snapshot we want to analyze
s_pathfile = f'{s_path}/{s_file}'
print('mcds xml:', s_pathfile)

# load mcds - multi cell data standard - object
mcds = pcdl.TimeStep(s_pathfile)  # loads the whole snapshot: the xml and all related mat and graph files
```

Side note: for path, in general, unix (slash) and windows (backslash) notation will work.

The legacy way of loading data, where filename and path had to be separated, works too.

```python
# legacy way of loading a mcds object
mcds = pcdl.TimeStep('output00000012.xml', s_path)
```

By default, all data related to the snapshot is loaded.\
For speed and less memory usage, it is however possible to only load the essential (output xml and cell mat data), and exclude microenvironment, graph data, and ID label mapping loading.\
Additionally, it is possible to specify for custom\_data variable types other than the generic float type, namely: int, bool, and str.

```python
# fine tuned way of loading a mcds object
mcds = pcdl.TimeStep(s_pathfile, custom_data_type={}, microenv=False, graph=False, settingxml=None)
```


### Working With MCDS Time Series in Python3

An exciting thing about modeling is to have time series data.\
Pcdl's TimeSeries class is here to make the handling of a time series of MCD snapshots easy.

```python
import pathlib  # library to locate the test data
import pcdl  # the physicell data loader library

s_path = str(pathlib.Path(pcdl.__file__).parent.joinpath('data_timeseries_2d'))  # local path to the installed test data set
print('mcds time series:', s_path)

# generate a mcds time series instance
mcdsts = pcdl.TimeSeries(s_path)
```

Like in the pyMCDs class, for memory consumption and processing speed control, we can specify if we want to load microenvironment data and graph data from the snapshots we later on analyze.\
By default, all data related to the snapshot is loaded, if needed, we can set load to False.\
Additionally, we can exclude to read the PhysiCell\_settings.xml file, if it is not available, and fine tune variable typing.

```python
# fine tuned the way mcds objects will be loaded
mcdsts = pcdl.TimeSeries(s_path, custom_data_type={}, load=True, microenv=False, graph=False, settingxml=None)
```


#### MCDS Times Series Data Analysis - the Basics

Now, let's have a look at the TimeSeries instances data analysis functions.

```python
# fetch data
mcdsts = pcdl.TimeSeries(s_path, load=False)
```

+ **get_xmlfile_list** will simply return a list of absolute path and mcds xml file name strings.\
This list can be manipulated, and later be fed into the objects read\_xml function, to only read specific snapshots into memory.
For example, if you want only to analyze hour 11, 12, and 13 from your run, or if data from every other hour will be enough.

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
The default setting will read all snapshots from a time series. \
However, you can always use a filtered ls\_xml list to only load a subset of snapshots.

```python
# load all snapshots
mcdsts.read_mcds()
mcdsts.get_mcds_list()  # YMMV! [<pcdl.pyMCDS.pyMCDS at 0x7fa660996b00>, <pcdl.pyMCDS.pyMCDS at 0x7fa67673e1d0>, ..., <pcdl.pyMCDS.pyMCDS at 0x7fa660950f10>]
len(mcdsts.get_mcds_list())  # 25

# load snapshot 11, 12, and 13
mcdsts.read_mcds(ls_xml_11_12_13)
len(mcdsts.get_mcds_list())  # 3

# load all even snapshots
mcdsts.read_mcds(ls_xml_even)
len(mcdsts.get_mcds_list())  # 13
```

Single snapshots can now be accessed by indexing.\
With a single snapshot, you work exactly in the same way as with an object loaded by TimeStep.

```python
# get the simulation time
mcdsts.get_mcds_list()[12].get_time()  # 720.0

# get the cell data frame
df = mcdsts.get_mcds_list()[12].get_cell_df()
df.shape  # (992, 87)
df.head()
```

It is very easy to loop over the whole time series, to analyze time series data, and generate time series plots.

```python
# load library
import matplotlib.pyplot as plt
import pandas as pd

# fetch data
mcdsts = pcdl.TimeSeries(s_path)

# loop over the time series to gather temporal information like, for example, data for a growth curve
lr_time = [mcds.get_time() for mcds in mcdsts.get_mcds_list()]  # [0.0, 60.0, ..., 1440.0]
li_cellcount = [mcds.get_cell_df().shape[0] for mcds in mcdsts.get_mcds_list()]  # [889, 898, ..., 1099]

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


#### Transform MCDS Time Series into AnnData Objects

```python
# fetch data
mcdsts = pcdl.TimeSeries(s_path)
```

A mcds time series can be translated into one single anndata (default).
```python
ann = mcdsts.get_anndata(values=2, scale='maxabs', collapse=True)
print(ann)  # AnnData object with n_obs × n_vars = 889 × 26
            #     obs: 'z_layer', 'time', 'current_phase', 'cycle_model'
            #     obsm: 'spatial'
```
The output tells us that we have loaded a time series with 24758 cell (agent) snapshots and 27 attributes.
And that we have spatial coordinate annotation (position\_x, position\_y, position\_z, time) of the loaded data.

A mcds time series can be translated into a chronological list of anndata objects, where each entry is a single time step.
After running get\_anndata, you can access the objects by the get\_annmcds\_list function.
```python
l_ann = mcdsts.get_anndata(values=2, scale='maxabs', collapse=False)
len(l_ann)  # 25
l_ann[24]  # AnnData object with n_obs × n_vars = 1099 × 79
           #     obs: 'z_layer', 'time', 'cell_type', 'current_death_model', 'current_phase', 'cycle_model'
           #     obsm: 'spatial'

# load all snapshots
mcdsts.get_annmcds_list()  # YMMV! [AnnData object with n_obs × n_vars ..., ..., ... obsm: 'spatial']
len(mcdsts.get_annmcds_list())  # 25
```


#### Times Series Data Scatter Plot Images, Contour Plot Images, and Movies

With PhysiCell it is not only possible to take data snapshots, but as well [svg](https://en.wikipedia.org/wiki/SVG) vector graphics images snapshots.\
PhysiCell's [Makefile](https://en.wikipedia.org/wiki/Make_(software)) has code to translate those svg images into [gif](https://en.wikipedia.org/wiki/GIF), [jpeg](https://en.wikipedia.org/wiki/JPEG), [png](https://en.wikipedia.org/wiki/Portable_Network_Graphics), or [tiff](https://en.wikipedia.org/wiki/TIFF) format, making use of the [image magick](https://en.wikipedia.org/wiki/ImageMagick) library.
The Makefile also has  code to translate the jpeg, png, or tiff images into a [mp4](https://en.wikipedia.org/wiki/MP4_file_format) movie, therefore utilizing the [ffmpeg](https://en.wikipedia.org/wiki/FFmpeg) library.\
TimeSeries instances provide similar functionality, although the jpeg, png, and tiff images are generated straight from data and not from the svg files.
However, mp4 movies and gif images are generated in the same way.
This means the mcdsts.make\_gif code will only run if image magick and mcdsts.make\_movie code will only run if ffmpeg is installed on your computer.

```python
# fetch data
mcdsts = pcdl.TimeSeries(s_path)
```

Translate physicell output data into static raster graphic images:

```python
# cell images colored by cell_type
mcdsts.plot_scatter()  # jpeg images colored by cell_type
mcdsts.plot_scatter(ext='png')  # png images colored by cell_type
mcdsts.plot_scatter(ext='tiff')  # tiff images colored by cell_type

# cell images can be colored by any cell dataframe column, like, for example, by pressure
# there are parameters to adjust cell size and more
mcdsts.plot_scatter(focus='pressure')  # jpeg images colored by pressure values

# substrate data owns an equivalent function,
# being able to display any substrate dataframe column
mcdsts.plot_contour(focus='oxygen')  # jpeg images colored by oxygen values
```

Translate raster graphic images into a dynamic gif image:\
Gif images can only be generated for already existing jpeg, png, or tiff images!\
By default, jpeg files will be used, to generate the gif.\
If png or tiff files should be used as source, then this has to be explicitly stated.

```python
# using jpeg images for input
s_path = mcdsts.plot_scatter()
mcdsts.make_gif(s_path)

# using jpeg images one-liner
mcdsts.make_gif(mcdsts.plot_scatter())

# using tiff images for input
s_path = mcdsts.plot_scatter(ext='tiff')
mcdsts.make_gif(s_path, interface='tiff')
```

Translate raster graphic images into an mp4 movie:\
Movies can only be generated for already existing jpeg, png, or tiff images!

```python
# using jpeg images for input
s_path = mcdsts.plot_scatter()
mcdsts.make_movie(s_path)

# using jpeg images one-liner
mcdsts.make_movie(mcdsts.plot_scatter())

# using tiff images for input
s_path = mcdsts.plot_scatter(ext='tiff')
mcdsts.make_movie(s_path, interface='tiff')
```


#### MCDS Times Series Data Triage

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

Cell variables that have no variance, zero entropy, that are in all agent overall time steps in the same state, have always exacted the same value, carry no information.
Similarly, substrates variables that over the whole domain overall time steps have the same concentration are not interesting.

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






### Data Clean Up

After you are done checking out the results, you can uninstall the test datasets and all files stored within its folders.

```
pcdl.uninstall_data()
```
