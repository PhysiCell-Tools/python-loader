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

The links from the most important data analysis functions are highlited as **#! workhorse functions** . \
Familiarize yourself well with their parameters!


# TimeStep Class
+ help([pcdl.TimeStep.__init__](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.__init__.md))  #! workhorse functions

## TimeStep medata
+ help([mcds.get_multicellds_version](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_multicellds_version.md))
+ help([mcds.get_physicell_version](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_physicell_version.md))
+ help([mcds.get_timestamp](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_timestamp.md))
+ help([mcds.get_time](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_time.md))
+ help([mcds.get_runtime](https://github.com/elmbeech/physicelldataloader/tree/master/man/docstring/mcds.get_runtime.md))

## TimeStep mesh
help(pcdl.TimeStep.get_voxel_ijk_range)
help(pcdl.TimeStep.get_mesh_mnp_range)
help(pcdl.TimeStep.get_xyz_range)
help(pcdl.TimeStep.get_voxel_ijk_axis)
help(pcdl.TimeStep.get_mesh_mnp_axis)
help(pcdl.TimeStep.get_mesh)
help(pcdl.TimeStep.get_mesh_2D)
help(pcdl.TimeStep.get_mesh_coordinate)
help(pcdl.TimeStep.get_mesh_spacing)
help(pcdl.TimeStep.get_voxel_spacing)
help(pcdl.TimeStep.get_voxel_volume)
help(pcdl.TimeStep.get_voxel_ijk)
help(pcdl.TimeStep.is_in_mesh)

## TimeStep microenvironment
help(pcdl.TimeStep.get_substrate_names)
help(pcdl.TimeStep.get_substrate_dict) cli ??? nop.
help(pcdl.TimeStep.get_substrate_df); cli ??? nop this is input.
help(pcdl.TimeStep.get_concentration) cli ??? nop this is bs.
help(pcdl.TimeStep.get_concentration_df)  # ! workhorse function
help(pcdl.TimeStep.get_conc_df)  # ! shorthand; cli
help(pcdl.TimeStep.get_concentration_at)
help(pcdl.TimeStep.plot_contour)  # ! workhorse function; (cli)

## TimeStep cells and other agents
help(pcdl.TimeStep.get_celltype_dict)
help(pcdl.TimeStep.get_cell_variables)
help(pcdl.TimeStep.get_cell_df)  # ! workhorse function; cli
help(pcdl.TimeStep.get_cell_df_at)
help(pcdl.TimeStep.plot_scatter)  # ! workhorse function; (cli)

## TimeStep graphs
help(pcdl.TimeStep.get_attached_graph_dict)
help(pcdl.TimeStep.get_neighbor_graph_dict)
help(pcdl.TimeStep.get_graph_gml)  # ! workhose function; (cli)

## TimeStep unit
help(pcdl.TimeStep.get_unit_se)  # ! workhorse function; cli

## TimeStep anndata
help(pcdl.TimeStep.get_anndata)  # ! workhorse function; nop!

## TimeStep internal functions
help(pcdl.TimeStep._anndextract)
help(pcdl.TimeStep._read_xml)
help(pcdl.TimeStep.graphfile_parser)
help(pcdl.TimeStep.scaler)


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

# TimeSeries internal functions
help(pcdl.TimeSeries._handle_magick)
```


# test data sets
```python3
help(pcdl.install_data)
help(pcdl.uninstall_data)
```
