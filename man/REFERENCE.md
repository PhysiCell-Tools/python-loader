# PhysiCell Data Loader Reference Man Pages

This is the technical descriptions of the machinery and how to operate it.\
References are maintained straight in the [source code](https://github.com/elmbeech/physicelldataloader/tree/master/pcdl),
in each function's [docstring](https://en.wikipedia.org/wiki/Docstring),
or in the case of the command line interface functions in the [argparse](https://docs.python.org/3/library/argparse.html) help strings.

Docstings can be accessed by clicking on the commands in the listing below,
or within the python3 shell with the `help(command)` ,
or in the ipython3 shell with the `help(command)` or `command?` .

For the command line interface functions, you can access help by clicking on the commands in the listing below,
or straight from the command line with the `command --help` or `command -h` [help](https://en.wikipedia.org/wiki/Help!) parameter.

Within python3, you can load the PhysiCell data loader module like this:
```python3
import pcdl

mcds = pcdl.TimeStep('path/to/outputnnnnnnnn.xml')
mcdsts = pcdl.TimeStep('path/to/output')
```
The links to the most important data analysis functions are highlighted as **#! workhorse function** . \
Familiarize yourself well with their parameters!


# Test datasets
+ help([pcdl.install_data](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl.install_data.md))
+ help([pcdl.uninstall_data](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl.uninstall_data.md))


# TimeStep

There are basically four types of functions:
+ set_ : set a python3 variable.
+ get_ : recall a python3 variable.
+ make_ : make functions generate file output.
+ plot_ : plot functions generate a matplotlib figure, or axis object, or file output, depending on your parameter settings.

## TimeStep initialize
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
+ help([mcds.get_anndata](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_anndata.md))  #! workhorse function

## TimeStep graphs
+ help([mcds.get_attached_graph_dict](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_attached_graph_dict.md))
+ help([mcds.get_neighbor_graph_dict](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_neighbor_graph_dict.md))
+ help([mcds.make_graph_gml](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.make_graph_gml.md))  #! workhose function

## TimeStep settings
+ help([mcds.get_parameter_dict](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_parameter_dict.md))  #! workhorse function
+ help([mcds.get_rule_df](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_rule_df.md))  #! workhorse function
+ help([mcds.get_unit_se](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_unit_se.md))  #! workhorse function

## TimeStep internal functions
+ help([mcds.scaler](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl.scaler.md))  # anndata
+ help([mcds.graphfile_parser](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl.graphfile_parser.md))  # mcds
```python3
help(pcdl.pyMCDS._read_xml)
help(pcdl.pyAnnData._anndextract)
```


# TimeSeries

There are basically four types of functions:
+ set_ : set a python3 variable.
+ get_ : recall a python3 variable.
+ make_ : make functions generate file output.
+ plot_ : plot functions generate a matplotlib figure, or axis object, or file output, depending on your parameter settings.

## TimeSeries initialization
+ help([pcdl.TimeSeries.__init__](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcdsts.__init__.md)) #! workhosefunction
+ help([mcdsts.get_xmlfile_list](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcdsts.get_xmlfile_list.md))  #! workhosefunction
+ help([mcdsts.read_mcds](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcdsts.read_mcds.md))
+ help([mcdsts.get_mcds_list](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcdsts.get_mcds_list.md))  #! workhose function
+ help([mcdsts.get_annmcds_list](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcdsts.get_annmcds_list.md))  #! workhose function
+ help([mcdsts.set_verbose_false](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcdsts.set_verbose_false.md)) #! workhorse function
+ help([mcdsts.set_verbose_true](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcdsts.set_verbose_true.md)) #! workhorse function

## TimeSeries microenvironment
+ help([mcdsts.get_conc_df_features](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcdsts.get_conc_df_features.md))  #! workhorse function
+ help([mcdsts.plot_contour](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.plot_contour.md))  # !workhorse function

## TimeSeries cells and other agents
+ help([mcdsts.get_cell_df_features](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcdsts.get_cell_df.md))  #! workhorse function
+ help([mcdsts.plot_scatter](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcdsts.plot_scatter.md))  # !workhorse function
+ help([mcdsts.get_anndata](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcdsts.get_anndata.md))  #! workhorse function

## TimeSeries graphs
+ help([mcdsts.get_graph_gml](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.make_graph_gml.md))  #! workhose function

## Timeseries timeseries
+ help([mcdsts.plot_timeseries](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcdsts.plot_timeseries.md))  #! workhorse function

## TimeSeries making movies
+ help([mcdsts.make_gif](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl.make_gif.md))  # ! workhorse function
+ help([mcdsts.make_movie](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl.make_movie.md))  # ! workhorse function

## TimeSeries internal functions
+ help([mcdsts.scaler](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl.scaler.md))  # anndata
```
help(pcdl.TimeSeries._handle_magick)
help(pcdl.pyAnnData._anndextract)
```


# Command line

The command line interface functions mimic the name and parameter arguments as closely as possible to the related python3 functions.

## Command line metadata
+ [pcdl_get_version --help](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl_get_version.md)  # ! workhorse function

## Command line microenvironment
+ [pcdl_get_conc_df_features --help](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl_get_conc_df_features.md)  # ! workhorse function
+ [pcdl_get_conc_df --help](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl_get_conc_df.md)  # ! workhorse function
+ [pcdl_plot_contour --help](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl_plot_contour.md)  # ! workhorse function

## Command line cells and other agents
+ [pcdl_get_cell_df_features --help](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl_get_cell_df_features.md)  # ! workhorse function
+ [pcdl_get_cell_df --help](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl_get_cell_df.md)  # ! workhorse function
+ [pcdl_plot_scatter --help](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl_plot_scatter.md)  # ! workhorse function
+ [pcdl_get_anndata --help](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl_get_anndata.md)  # ! workhorse function

## Command line graph
+ [pcdl_make_graph_gml --help](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl_make_graph_gml.md)  # ! workhorse function

## Command line setting
+ [pcdl_get_parameter_dict --help](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl_get_parameter_dict.md)  # ! workhorse function
+ [pcdl_get_rule_df --help](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl_get_rule_df.md)  # ! workhorse function
+ [pcdl_get_unit_se --help](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl_get_unit_se.md)  # ! workhorse function

## Command line timeseries
+ [pcdl_plot_timeseries --help](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl_plot_timeseries.md)  # ! workhorse function

## Command line making movies
+ [pcdl_make_gif --help](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl_make_gif.md)  # ! workhorse function
+ [pcdl_make_movie --help](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/pcdl_make_movie.md)  # ! workhorse function

