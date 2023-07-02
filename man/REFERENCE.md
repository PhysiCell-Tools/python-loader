# PhysiCell Data Loader Reference Man Page

This is the technical descriptions of the machinery and how to operate it.\
References are maintained in each module's [docstring](https://en.wikipedia.org/wiki/Docstring).\
You can access them through the [source code](https://github.com/elmbeech/physicelldataloader/tree/master/pcdl), or by first load PhysiCell Data Loader.

```python3
import pcdl
```

Then, for each pcdl module, get on the fly reference information with the [help](https://en.wikipedia.org/wiki/Help!) command.

# pyMCDS
```
help(pcdl.pyMCDS)
help(pcdl.pyMCDS.__init__)

# pyMCDS medata
help(pcdl.pyMCDS.get_multicellds_version)
help(pcdl.pyMCDS.get_physicell_version)
help(pcdl.pyMCDS.get_timestamp)
help(pcdl.pyMCDS.get_time)
help(pcdl.pyMCDS.get_runtime)

# pyMCDS mesh
help(pcdl.pyMCDS.get_voxel_ijk_range)
help(pcdl.pyMCDS.get_mesh_mnp_range)
help(pcdl.pyMCDS.get_xyz_range)
help(pcdl.pyMCDS.get_voxel_ijk_axis)
help(pcdl.pyMCDS.get_mesh_mnp_axis)
help(pcdl.pyMCDS.get_mesh)
help(pcdl.pyMCDS.get_mesh_2D)
help(pcdl.pyMCDS.get_mesh_coordinate)
help(pcdl.pyMCDS.get_mesh_spacing)
help(pcdl.pyMCDS.get_voxel_spacing)
help(pcdl.pyMCDS.get_voxel_volume)
help(pcdl.pyMCDS.get_voxel_ijk)
help(pcdl.pyMCDS.is_in_mesh)

# pyMCDS microenvironment
help(pcdl.pyMCDS.get_substrate_names)
help(pcdl.pyMCDS.get_substrate_dict)
help(pcdl.pyMCDS.get_substrate_df)
help(pcdl.pyMCDS.get_concentration)
help(pcdl.pyMCDS.get_concentration_df)
help(pcdl.pyMCDS.get_concentration_at)
help(pcdl.pyMCDS.get_contour)

# pyMCDS cells and other agents
help(pcdl.pyMCDS.get_celltype_dict)
help(pcdl.pyMCDS.get_cell_variables)
help(pcdl.pyMCDS.get_cell_df)
help(pcdl.pyMCDS.get_cell_df_at)

# pyMCDS graphs
help(pcdl.pyMCDS.get_attached_graph_dict)
help(pcdl.pyMCDS.get_neighbor_graph_dict)

# pyMCDS unit
help(pcdl.pyMCDS.get_unit_df)

# pyMCDS internal functions
help(pcdl.pyMCDS._read_xml)
help(pcdl.graphfile_parser)
```

# pyMCDS time series
```
help(pcdl.pyMCDSts)
help(pcdl.pyMCDSts.__init__)

# pyMCDSts load data
help(pcdl.pyMCDSts.get_xmlfile_list)
help(pcdl.pyMCDSts.read_mcds)

# pyMCDSts triage data
help(pcdl.pyMCDSts.get_cell_minstate_col)
help(pcdl.pyMCDSts.get_concentration_minstate_col)

# pyMCDSts images and movies
help(pcdl.pyMCDSts.make_imgcell)
help(pcdl.pyMCDSts.make_imgsubs)
help(pcdl.pyMCDSts.make_gif)
help(pcdl.pyMCDSts.make_movie)

# pyMCDSts internal functions
help(pcdl.pyMCDSts._handle_magick)
```

# test data sets
```
help(pcdl.install_data)
help(pcdl.uninstall_data)
```
