# PhysiCell Data Loader Reference Man Pages

This is the technical descriptions of the machinery and how to operate it.\
References are maintained in straight on the [source code](https://github.com/elmbeech/physicelldataloader/tree/master/pcdl), \
in each function's [docstring](https://en.wikipedia.org/wiki/Docstring), \
or in case of the command Line interface functions in the [argparse](https://docs.python.org/3/library/argparse.html) help strings.

Docstings you can access by clicking on the commands in the listing below, \
or within the python shell with the `help(command)` , \
or in the ipython shell with the `help(command)` or `command?` .

For the command line interface functions you can access help by clickin on the commands in the listing below, \
or straight from the command line with `command --help` or `command -h` [help](https://en.wikipedia.org/wiki/Help!) parameter.


Within python3, you can load the PhysiCell data loader module like this:
```python3
import pcdl

mcds = pcdl.TimeStep('path/to/outputnnnnnnnn.xml')
mcdsts = pcdl.TimeStep('path/to/output')
```

The links from the most important data analysis functions are highlited as **#! workhorse function** . \
Familiarize yourself well with their parameters!

There are four types of function:
+ set_ : set a python3 varaiale.
+ get_ : recal a python3 variable.
+ make_ : make functions generates file output.
+ plot_ : plot functions enerate a matplotlib figure or axis object or file output, depending on your parameter settings.


# TimeStep initialize
+ help([pcdl.TimeStep.__init__](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.__init__.md))  #! workhorse function
+ help([mcds.set_verbose_false](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.set_verbose_false.md)) #! workhorse function
+ help([mcds.set_verbose_true](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.set_verbose_true.md)) #! workhorse function


## TimeStep medata
**version**
+ help([mcds.get_multicellds_version](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_multicellds_version.md))
+ help([mcds.get_physicell_version](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_physicell_version.md))
**time**
+ help([mcds.get_timestamp](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_timestamp.md))
+ help([mcds.get_time](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_time.md))
+ help([mcds.get_runtime](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_runtime.md))


## TimeStep mesh

**range and axis**
+ help([mcds.get_voxel_ijk_range](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_voxel_ijk_range.md))
+ help([mcds.get_mesh_mnp_range](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_mesh_spacing.md))
+ help([mcds.get_xyz_range](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_xyz_range.md))
+ help([mcds.get_voxel_ijk_axis](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_voxel_ijk_axis.md))
+ help([mcds.get_mesh_mnp_axis](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_mesh_mnp_axis.md))

**mesh mnp**
+ help([mcds.get_mesh](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_mesh.md))
+ help([mcds.get_mesh_2D](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_mesh_2D.md))
+ help([mcds.get_mesh_coordinate](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_mesh_coordinate.md))
+ help([mcds.get_mesh_spacing](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_mesh_spacing.md))
+ help([mcds.is_in_mesh](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.is_in_mesh.md))

**voxel ijk**
+ help([mcds.get_voxel_spacing](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_voxel_spacing.md))
+ help([mcds.get_voxel_volume](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_voxel_volume.md))
+ help([mcds.get_voxel_ijk](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_voxel_ijk.md))


## TimeStep microenvironment
+ help([mcds.get_substrate_names](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_substrate_names.md))
+ help([mcds.get_substrate_dict](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_substrate_dict.md))
+ help([mcds.get_substrate_df](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_substrate_df.md))
+ help([mcds.get_concentration](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_concentration.md))
+ help([mcds.get_concentration_at](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_concentration_at.md))
+ help([mcds.get_conc_df](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_conc_df.md))  #! workhorse function
+ help([mcds.plot_contour](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.plot_contour.md))  #! workhorse function

## TimeStep cells and other agents
+ help([mcds.get_cell_variables](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_cell_variables.md))
+ help([mcds.get_celltype_dict](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_celltype_dict.md))
+ help([mcds.get_cell_df_at](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_cell_df_at.md))
+ help([mcds.get_cell_df](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_cell_df.md))  #! workhorse function
+ help([mcds.plot_scatter](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.plot_scatter.md))  #! workhorse function

## TimeStep graphs
+ help([mcds.get_attached_graph_dict](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_attached_graph_dict.md))
+ help([mcds.get_neighbor_graph_dict](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_neighbor_graph_dict.md))
+ help([mcds.make_graph_gml](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.make_graph_gml.md))  #! workhose function

## TimeStep unit
+ help([mcds.get_unit_se](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_unit_se.md))  #! workhorse function

## TimeStep anndata
+ help([mcds.get_anndata](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_anndata.md))  #! workhorse function


## TimeStep internal functions
```python
help(pcdl.pyMCDS.\_read\_xml)
help(pcdl.pyMCDS.graphfile\_parser)
help(pcdl.pyAnnData.\_anndextract)
help(pcdl.pyAnnData.scaler)
```

# TimeSeries Class
```python3
help(pcdl.TimeSeries)  # ! make class instance
help(pcdl.TimeSeries.__init__)

# TimeSeries load data
help(pcdl.TimeSeries.get_xmlfile_list)
help(pcdl.TimeSeries.read_mcds)
help(pcdl.TimeSeries.get_mcds_list)  # ! workhorse function
help(pcdl.TimeSeries.get_annmcds_list)  # ! workhorse function

# TimeSeries triage data
help(pcdl.TimeSeries.get_cell_df_features)  # ! workhorse function; cli
help(pcdl.TimeSeries.get_conc_df_features)  # ! workhorse function; cli

# TimeSeries images and movies
help(pcdl.TimeSeries.plot_timeseries)  # ! workhorse function; cli
help(pcdl.TimeSeries.plot_scatter)  # ! workhorse function; cli
help(pcdl.TimeSeries.plot_contour)  # ! workhorse function; cli
help(pcdl.TimeSeries.make_gif)  # ! workhorse function; cli (no-need for mcds?)
help(pcdl.TimeSeries.make_movie)  # ! workhorse function; cli (no-need for mcds?)

# TimeSeries graphs
help(pcdl.TimeSeries.get_graph_gml)  # ! workhose function; cli

# TimeSeries anndata
help(pcdl.TimeSeries.get_anndata)  # ! workhorse function
help(pcdl.TimeSeries.get_annmcds_list)

# TimeSeries internal functions
help(pcdl.TimeSeries._handle_magick)
```


# test data sets
```python3
help(pcdl.install_data)
help(pcdl.uninstall_data)
```
