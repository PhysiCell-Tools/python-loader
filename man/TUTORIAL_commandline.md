# PhysiCell Data Loader Tutorial: pcdl from the Commandline

<!-- bue 2024-08-22: have to check if this works from dos and power shell. it will work somehow, because i can run the cli unit tests -->

The most important commands for down stream data analysis,
available in the pcdl TimeStep and TimeSeries class,
can be run straight from a command line shell like: [bash](https://en.wikipedia.org/wiki/Bash_(Unix_shell)), [csh](https://en.wikipedia.org/wiki/C_shell), [dos](https://en.wikipedia.org/wiki/DOS), [fish](https://en.wikipedia.org/wiki/Fish_(Unix_shell)), [ksh](https://en.wikipedia.org/wiki/KornShell), [powershell](https://en.wikipedia.org/wiki/PowerShell), [sh](https://en.wikipedia.org/wiki/Bourne_shell), [tsh](https://en.wikipedia.org/wiki/Tcsh), or [zsh](https://en.wikipedia.org/wiki/Z_shell), to name a view.

The command names are drived from the related python3 function. \
The command parameter mimic the related python3 function arguments as closely as possible. \
You can always call the [help](https://en.wikipedia.org/wiki/Help!) parameter ( pcdl_command -h),
to access the man page for a pcdl command!

Please spend some time to learn about each of the 18 commands, by studing it's man page.
this will truely make you a power user!


## Preparation

To runs this tutorial,
you can install the 2D unittest dataset into your PhysiCell output folder,
by executing the following command sequence.

**Waring: all data currentlty in your PhysiCell/output folder will be overwritten!**

```bash
cd path/to/PhysiCell
```
```bash
make data-cleanup
python3 -c"import pathlib, pcdl, shutil; s_ipath=str(pathlib.Path(pcdl.__file__).parent.resolve()/'data_timeseries_2d'); shutil.copytree(s_ipath, 'output', dirs_exist_ok=True)"
```


## Metadata related commands

### &#x2728 pcdl\_get\_version

Output PhysiCell, MCDS, and pcdl version on the screen.

```bash
pcdl_get_version output
```
```bash
pcdl_get_version output/output00000000.xml
```
```bash
pcdl_get_version -h
```

### pcdl\_get\_unit dict

Generate a [csv](https://en.wikipedia.org/wiki/Comma-separated_values) file that maps attribute and units, as specified in the settings.xml.

```bash
pcdl_get_unit_dict output
```
```bash
pcdl_get_unit_dict output/output00000000.xml
```
```bash
pcdl_get_unit_dict -h
```


## Microenvironment realted commands


### pcdl\_get\_substrate\_list

Output all substrated modeled in the microenviroment.

```bash
pcdl_get_substrate_list output
```
```bash
pcdl_get_substrate_list output/output00000000.xml
```
```bash
pcdl_get_substrate_list -h
```

### pcdl\_get\_conc\_attribute

Generate a [json](https://en.wikipedia.org/wiki/JSON) file, that lists all substrate attributes.
For each such attribute the min and the max value are listed.

In the example below:
+ all substrates attributes are listed, that over the whole timeseris have at least 2 diffrent values.
+ all substrates attributes are listed, that in this particular time step have at least 2 different values.

```bash
pcdl_get_conc_attribute output 2
```
```bash
pcdl_get_conc_attribute output/output00000000.xml 2
```
```bash
pcdl_get_conc_attribute -h
```

<!-- bue 20240822: link extra tutorials -->


### pcdl\_get\_conc\_df

Generate a [csv](https://en.wikipedia.org/wiki/Comma-separated_values) file that lists one voxel per row,
all substarte concentrations.

In the example below the generated csv contains:
+ all substarte concentration values over the whole timeseries.
+ all substarte that within this particular time step have more than 2 different concentration values.

```bash
pcdl_get_conc_df output
```
```bash
pcdl_get_conc_df output/output00000000.xml 2
```
```bash
pcdl_get_conc_df -h
```

<!-- bue 20240822: link extra tutorials -->


### Substrate concentration visualization

For oxygen generate a [jpeg](https://en.wikipedia.org/wiki/JPEG) file
for the the whole time series or for a single time step.

```bash
pcdl_plot_contour output oxygen
pcdl_plot_contour output/output00000000.xml oxygen
pcdl_plot_contour -h
```

Generate a rectilinear grid [vtk](https://en.wikipedia.org/wiki/VTK) file from a single time step,
or rectilinear grid vtk files from the whole time series,
containing all the substrates from the model.

These vtk files can be further analysed,
for example with the [paraview](https://www.paraview.org/) or [blender](https://www.blender.org/) software,
as described in the extra tutorials.

```bash
pcdl_make_conc_vtk output
pcdl_make_conc_vtk output/output00000000.xml
pcdl_make_conc_vtk -h
```

<!-- bue 20240822: link extra tutorials -->


## Cell agent related commands


### Cell agent attributes

Generate a [json](https://en.wikipedia.org/wiki/JSON) file, that lists all cell attributes.
For each such attribute the min and the max value are listed.

In the example below:
+ all cell attributes are listed, that over the whole timeseris have at least 2 diffrent values.
+ all cell attributes are listed, that in this particular time step have at least 2 different values.

```bash
pcdl_get_cell_attribute output 2
pcdl_get_cell_attribute output/output00000000.xml 2
pcdl_get_cell_attribute -h
```

<!-- bue 20240822: link extra tutorials -->


### Cell agent dataframe

Generate a [csv](https://en.wikipedia.org/wiki/Comma-separated_values) file that lists one cell per row,
all attributes.

In the example below the generated csv contains:
+ all cell attributes, that over the whole time series have more than 2 different values.
+ from that particular time step all cell attributes available.

```bash
pcdl_get_cell_df output 2
pcdl_get_cell_df output/output00000000.xml
pcdl_get_cell_df -h
```

<!-- bue 20240822: link extra tutorials -->


### Cell agent anndata hd5 file

Frome the whole timeseries or from a single time step generate h5ad [anndata](https://anndata.readthedocs.io/en/latest/) [hd5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) files.

Anndata is the standard data format in the python single cell community.
Data stored in this fromat can be analyzed the same way as usually sc RNA seq data is analysed.

```bash
pcdl_get_anndata output
pcdl_get_anndata output/output00000000.xml
pcdl_get_anndata -h
```

<!-- bue 20240822: link extra tutorials -->


### Cell agent neighborhood graph file

Generate [gml](https://github.com/elmbeech/physicelldataloader/blob/master/man/publication/himsolt1996gml_a_portable_graph_file_format.pdf) files.
One gml file per time step.

Gml files can be read by graph analysis libraries like [networkx](https://networkx.org/) and [igraph](https://igraph.org/).

```bash
pcdl_make_graph_gml output
pcdl_make_graph_gml output/output00000000.xml --node_attribute cell_type dead oxygen pressure
pcdl_make_graph_gml -h
```

<!-- bue 20240822: link extra tutorials -->


### Cell agent visualization

Generate a [jpeg](https://en.wikipedia.org/wiki/JPEG) file that displaying all cells.

In the example below:
+ generate plots for the the whole time series, color the cells by cell\_type.
+ generate a plot for a single timestep, color the cells by pressure.

```bash
pcdl_plot_scatter output
pcdl_plot_scatter output/output00000000.xml pressure
pcdl_plot_scatter -h
```

Generate a 3D glyph [vtk](https://en.wikipedia.org/wiki/VTK) file from a single time step,
or rectilinear grid vtk files from the whole time series,
with information for the attributes listed.
The deafalt atribute listed is cell_type.

These vtk files can be further analysed,
for example with the [paraview](https://www.paraview.org/) or [blender](https://www.blender.org/) software,
as described in the extra tutorials.

```bash
pcdl_make_cell_vtk output
pcdl_make_cell_vtk output/output00000000.xml cell_type dead oxygen pressure
pcdl_make_cell_vtk -h
```

<!-- bue 20240822: link extra tutorials -->


## Microenvironment and cell agent related commands


### Microenviroment and cell agent visualization

Generate a timeseries plot and save it as a [jpeg](https://en.wikipedia.org/wiki/JPEG) file.

The default plots outputs total cell count over time.
However, this is a very powerfull command!
Below we generate a plot for:
+ total cell count
+ cell count per cell_type
+ mean oxygen concentration detected per cell_type
+ max oxygen concentration detected per cell_type
+ mean oxygen concentration detected in the cells.
+ mean oxygen concentration detected in the domain.

```bash
pcdl_plot_timeseries output
pcdl_plot_timeseries output none
pcdl_plot_timeseries output cell_type
pcdl_plot_timeseries output cell_type oxygen
pcdl_plot_timeseries output cell_type oxygen max
pcdl_plot_timeseries output none oxygen
pcdl_plot_timeseries output none oxygen --frame conc
pcdl_plot_timeseries -h
```

Generate a [ome.tiff](https://ome-model.readthedocs.io/en/stable/index.html) file,
to analyse a single time step or the whole time series,
the same way as usually fluorescent microscopy data is analysed.

By default, the cell_attribute outputted is the cell ID + 1.
However, any numerical (bool, int, float) cell_attribute can be outputted.
For example: dead, cells_per_voxel, or pressure.

These ome.tiff files can be further analysed,
for example with the [napari](https://napari.org/stable/) or [fiji imagej](https://fiji.sc/) software,
as described in the extra tutorials.

```bash
pcdl_make_ome_tiff output
pcdl_make_ome_tiff output/output00000000.xml pressure
pcdl_make_ome_tiff -h
```

<!-- bue 20240822: link extra tutorials -->


## Making movies

Make an [mp4](https://en.wikipedia.org/wiki/MP4_file_format) movie from the jpeg plots from a time series.

```bash
pcdl_plot_scatter output
pcdl_make_movie output/cell_cell_type_z0.0/
pcdl_make_movie -h
```

Make a [gif](https://en.wikipedia.org/wiki/GIF) image from the jpeg plots from a time series.

```bash
pcdl_plot_scatter output
pcdl_make_gif output/cell_cell_type_z0.0/
pcdl_make_gif -h
```
