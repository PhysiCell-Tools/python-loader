# PhysiCell Data Loader Reference Man Page

This is the technical descriptions of the machinery and how to operate it.\
References are maintained in each module's [docstring](https://en.wikipedia.org/wiki/Docstring).\
You can access them through the [source code](https://github.com/elmbeech/physicelldataloader/tree/master/pcdl), or by first load PhysiCell Data Loader.

```python3
import pcdl
```

Then, for each pcdl module, get on the fly reference information with the [help](https://en.wikipedia.org/wiki/Help!) command.\
The **workhorse functions** are the ones most important for data analysis.
Familiarize yourself with all of their parameters!

# TimeStep
```python3
help(pcdl.TimeStep)  # ! make class instance
help(pcdl.TimeStep.__init__)

# TimeStep medata
help(pcdl.TimeStep.get_multicellds_version)
help(pcdl.TimeStep.get_physicell_version)
help(pcdl.TimeStep.get_timestamp)
help(pcdl.TimeStep.get_time)
help(pcdl.TimeStep.get_runtime)

# TimeStep mesh
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

# TimeStep microenvironment
help(pcdl.TimeStep.get_substrate_names)
help(pcdl.TimeStep.get_substrate_dict)
help(pcdl.TimeStep.get_substrate_df)
help(pcdl.TimeStep.get_concentration)
help(pcdl.TimeStep.get_concentration_df)  # ! workhorse function
help(pcdl.TimeStep.get_conc_df)  # ! shorthand
help(pcdl.TimeStep.get_concentration_at)
help(pcdl.TimeStep.get_contour)  # ! workhorse function

# TimeStep cells and other agents
help(pcdl.TimeStep.get_celltype_dict)
help(pcdl.TimeStep.get_cell_variables)
help(pcdl.TimeStep.get_cell_df)  # ! workhorse function
help(pcdl.TimeStep.get_cell_df_at)

# TimeStep graphs
help(pcdl.TimeStep.get_attached_graph_dict)  # !
help(pcdl.TimeStep.get_neighbor_graph_dict)  # !

# TimeStep unit
help(pcdl.TimeStep.get_unit_df)  # ! workhorse function

# TimeStep anndata
help(pcdl.TimeStep.get_anndata)  # ! workhorse function

# TimeStep internal functions
help(pcdl.TimeStep._read_xml)
help(pcdl.graphfile_parser)
help(pcdl.extract)
help(pcdl.scaler)
```


# TimeStep time series
```
help(pcdl.TimeSeries)  # ! make class instance
help(pcdl.TimeSeries.__init__)

# TimeSeries load data
help(pcdl.TimeSeries.get_xmlfile_list)
help(pcdl.TimeSeries.read_mcds)

# TimeSeries triage data
help(pcdl.TimeSeries.get_cell_df_columns_min_states)
help(pcdl.TimeSeries.get_conc_df_columns_min_states)

# TimeSeries images and movies
help(pcdl.TimeSeries.make_imgcell)  # ! workhorse function
help(pcdl.TimeSeries.make_imgsubs)  # ! workhorse function
help(pcdl.TimeSeries.make_gif)  # ! workhorse function
help(pcdl.TimeSeries.make_movie)  # ! workhorse function

# TimeSeries anndata
help(pcdl.TimeSeries.get_anndata)  # ! workhorse function

# TimeSeries internal functions
help(pcdl.TimeSeries._handle_magick)
```


# test data sets
```
help(pcdl.install_data)
help(pcdl.uninstall_data)
```
