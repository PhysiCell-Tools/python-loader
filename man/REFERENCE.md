# pcDataLoader Reference Man Page

This is the technical descriptions of the machinery and how to operate it.\
References are maintained in each module's [docstring](https://en.wikipedia.org/wiki/Docstring).\
You can access them through the [source code](https://github.com/elmbeech/pcDataLoader/tree/master/pcDataLoader), or by first load PhysiCell Data Loader.

```python3
import pcDataLoader as pc
```

Then, for each pcDataLoader module, get on the fly reference information with the [help](https://en.wikipedia.org/wiki/Help!) command.

# pyMCDS
```
help(pc.pyMCDS)
help(pc.pyMCDS.__init__)

# pyMCDS medata
help(pc.pyMCDS.get_multicellds_version)
help(pc.pyMCDS.get_physicell_version)
help(pc.pyMCDS.get_timestamp)
help(pc.pyMCDS.get_time)
help(pc.pyMCDS.get_runtime)

# pyMCDS mesh
help(pc.pyMCDS.get_voxel_ijk_range)
help(pc.pyMCDS.get_mesh_mnp_range)
help(pc.pyMCDS.get_xyz_range)
help(pc.pyMCDS.get_voxel_ijk_axis)
help(pc.pyMCDS.get_mesh_mnp_axis)
help(pc.pyMCDS.get_mesh)
help(pc.pyMCDS.get_mesh_2D)
help(pc.pyMCDS.get_mesh_coordinate)
help(pc.pyMCDS.get_mesh_spacing)
help(pc.pyMCDS.get_voxel_spacing)
help(pc.pyMCDS.get_voxel_volume)
help(pc.pyMCDS.get_voxel_ijk)
help(pc.pyMCDS.is_in_mesh)

# pyMCDS microenvironment
help(pc.pyMCDS.get_substrate_names)
help(pc.pyMCDS.get_substrate_df)
help(pc.pyMCDS.get_concentration)
help(pc.pyMCDS.get_concentration_df)
help(pc.pyMCDS.get_concentration_at)
help(pc.pyMCDS.get_contour)

# pyMCDS cells and other agents
help(pc.pyMCDS.get_cell_variables)
help(pc.pyMCDS.get_cell_df)
help(pc.pyMCDS.get_cell_df_at)

# pyMCDS graphs
help(pc.pyMCDS.get_attached_graph_dict)
help(pc.pyMCDS.get_neighbor_graph_dict)

# pyMCDS unit
help(pc.pyMCDS.get_unit_df)

# pyMCDS internal functions
help(pc.pyMCDS._read_xml)
help(pc.graphfile_parser)
```

# pyMCDS time series
```
help(pc.pyMCDSts)
help(pc.pyMCDSts.__init__)

# pyMCDSts load data
help(pc.pyMCDSts.get_xmlfile_list)
help(pc.pyMCDSts.read_mcds)

# pyMCDSts images and movies
help(pc.pyMCDSts.make_jpeg)
help(pc.pyMCDSts.make_png)
help(pc.pyMCDSts.make_tiff)
help(pc.pyMCDSts.make_gif)
help(pc.pyMCDSts.make_movie)

# pyMCDSts internal functions
help(pc.pyMCDSts._handle_magick)
help(pc.pyMCDSts._handle_resize)
```

