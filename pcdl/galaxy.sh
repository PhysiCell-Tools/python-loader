## pcdl_get_anndata
#planemo tool_init --force \
#--id pcdl_get_anndata \
#--name pcdl_get_anndata \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_get_anndata output_3d 1 --custom_data_type --microenv true --graph true --physiboss true --settingxml PhysiCell_settings.xml --verbose false --drop --keep --scale maxabs --collapse true' \
#--example_input output_3d/ \
#--example_output output_3d/timeseries_cell_maxabs.h5ad \
#--test_case \
#--help_from_command 'pcdl_get_anndata --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

## pcdl_get_cell_attribute
#planemo tool_init --force \
#--id pcdl_get_cell_attribute \
#--name pcdl_get_cell_attribute \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_get_cell_attribute output_3d 1 --custom_data_type --microenv true --physiboss true --settingxml PhysiCell_settings.xml --verbose false --drop --keep --allvalues false' \
#--example_input output_3d/ \
#--example_output output_3d/timeseries_cell_attribute_minmax.json \
#--test_case \
#--help_from_command 'pcdl_get_cell_attribute --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

## pcdl_get_cell_attribute_list
#planemo tool_init --force \
#--id pcdl_get_cell_attribute_list \
#--name pcdl_get_cell_attribute_list \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_get_cell_attribute_list output_3d --microenv true --physiboss true --settingxml PhysiCell_settings.xml --verbose false 2> attribute.txt' \
#--example_input output_3d/ \
#--example_output version.txt \
#--test_case \
#--help_from_command 'pcdl_get_cell_attribute_list --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

## pcdl_get_cell_df
#planemo tool_init --force \
#--id pcdl_get_cell_df \
#--name pcdl_get_cell_df \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_get_cell_df output_3d 1 --microenv true --physiboss true --settingxml PhysiCell_settings.xml --verbose false --drop --keep --collapse true' \
#--example_input output_3d/ \
#--example_output output_3d/timeseries_cell.csv \
#--test_case \
#--help_from_command 'pcdl_get_cell_df --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

## pcdl_get_celltype_list
#planemo tool_init --force \
#--id pcdl_get_celltype_list \
#--name pcdl_get_celltype_list \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_get_celltype_list output_3d --settingxml PhysiCell_settings.xml --verbose false 2> celltype.txt' \
#--example_input output_3d/ \
#--example_output celltype.txt \
#--test_case \
#--help_from_command 'pcdl_get_celltype_list --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

## pcdl_get_conc_attribute
#planemo tool_init --force \
#--id pcdl_get_conc_attribute \
#--name pcdl_get_conc_attribute \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_get_conc_attribute output_3d 1 --verbose false --drop --keep --allvalues false' \
#--example_input output_3d/ \
#--example_output output_3d/timeseries_conc_attribute_minmax.json \
#--test_case \
#--help_from_command 'pcdl_get_conc_attribute --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

## pcdl_get_conc_df
#planemo tool_init --force \
#--id pcdl_get_conc_df \
#--name pcdl_get_conc_df \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_get_conc_df output_3d 1 --verbose false --drop --keep --collapse true' \
#--example_input output_3d/ \
#--example_output output_3d/timeseries_conc.csv \
#--test_case \
#--help_from_command 'pcdl_get_conc_df --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

## pcdl_get_substrate_list
#planemo tool_init --force \
#--id pcdl_get_substrate_list \
#--name pcdl_get_substrate_list \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_get_substrate_list output_3d --verbose false 2> substrate.txt' \
#--example_input output_3d/ \
#--example_output substrate.txt \
#--test_case \
#--help_from_command 'pcdl_get_substrate_list --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

## pcdl_get_unit_dict
#planemo tool_init --force \
#--id pcdl_get_unit_dict \
#--name pcdl_get_unit_dict \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_get_unit_dict output_3d --microenv true --settingxml PhysiCell_settings.xml --verbose false' \
#--example_input output_3d/ \
#--example_output output_3d/timeseries_unit.csv \
#--test_case \
#--help_from_command 'pcdl_get_unit_dict --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

# pcdl_get_version
#planemo tool_init --force \
#--id pcdl_get_version \
#--name pcdl_get_version \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_get_version output_3d --verbose false 2> version.txt' \
#--example_input output_3d/ \
#--example_output version.txt \
#--test_case \
#--help_from_command 'pcdl_get_version --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

## pcdl_make_cell_vtk
#rm output_3d/output00000000_cell.vtp
#planemo tool_init --force \
#--id pcdl_make_cell_vtk \
#--name pcdl_make_cell_vtk \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_make_cell_vtk output_3d/output00000000.xml cell_type --custom_data_type --microenv true --physiboss true --settingxml PhysiCell_settings.xml --verbose false' \
#--example_input output_3d/output00000000.xml \
#--example_output output_3d/output00000000_cell.vtp \
#--test_case \
#--help_from_command 'pcdl_make_cell_vtk --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

## pcdl_make_conc_vtk
#rm output_3d/output00000000_conc.vtr
#planemo tool_init --force \
#--id pcdl_make_conc_vtk \
#--name pcdl_make_conc_vtk \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_make_conc_vtk output_3d/output00000000.xml --verbose false' \
#--example_input output_3d/output00000000.xml \
#--example_output output_3d/output00000000_conc.vtr \
#--test_case \
#--help_from_command 'pcdl_make_conc_vtk --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

## pcdl_make_graph_gml
#rm output_3d/output00000012_neighbor.gml
#planemo tool_init --force \
#--id pcdl_make_graph_gml \
#--name pcdl_make_graph_gml \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_make_graph_gml output_3d/output00000012.xml neighbor --custom_data_type --microenv true --physiboss true --settingxml PhysiCell_settings.xml --verbose false --edge_attribute true --node_attribute' \
#--example_input output_3d/output00000012.xml \
#--example_output output_3d/output00000012_neighbor.gml \
#--test_case \
#--help_from_command 'pcdl_get_anndata --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

## pcdl_plot_contour
#rm output_3d/conc_oxygen_z-5.0/output00000012_oxygen.jpeg
#planemo tool_init --force \
#--id pcdl_plot_contour \
#--name pcdl_plot_contour \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_plot_contour output_3d/output00000012.xml oxygen --verbose false --z_slice 0.0 --extrema none --alpha 1.0 --fill true --cmap viridis --title "" --grid true --xlim none --ylim none --xyequal true --figsizepx none --ext jpeg --figbgcolor none' \
#--example_input output_3d/output00000012.xml \
#--example_output output_3d/conc_oxygen_z-5.0/output00000012_oxygen.jpeg \
#--test_case \
#--help_from_command 'pcdl_plot_contour --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

## pcdl_plot_scatter
#rm output_3d/cell_cell_type_z-5.0/output00000012_cell_type.jpeg
#planemo tool_init --force \
#--id pcdl_plot_scatter \
#--name pcdl_plot_scatter \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_plot_scatter output_3d/output00000012.xml cell_type --custom_data_type --microenv true --physiboss true --settingxml PhysiCell_settings.xml --verbose false --z_slice 0.0 --z_axis none --alpha 1.0 --cmap viridis --title "" --grid true --legend_loc "lower left" --xlim none --ylim none --xyequal true --s 1.0 --figsizepx none --ext jpeg --figbgcolor none' \
#--example_input output_3d/output00000012.xml \
#--example_output output_3d/cell_cell_type_z-5.0/output00000012_cell_type.jpeg \
#--test_case \
#--help_from_command 'pcdl_plot_scatter --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

## pcdl_plot_timeseries
#planemo tool_init --force \
#--id pcdl_plot_timeseries \
#--name pcdl_plot_timeseries \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_plot_timeseries output_3d none none mean --custom_data_type --microenv true --physiboss true --settingxml PhysiCell_settings.xml --verbose false --frame cell --z_slice 0.0 --logy false --ylim none --secondary_y false --subplots false --sharex false --sharey false --linestyle - --linewidth none --cmap none --color none --grid true --legend true --yunit none --title none --figsizepx 640 480 --ext jpeg --figbgcolor none' \
#--example_input output_3d/ \
#--example_output output_3d/timeseries_cell_total_count.jpeg \
#--test_case \
#--help_from_command 'pcdl_plot_timeseries --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

## pcdl_make_gif
#planemo tool_init --force \
#--id pcdl_make_gif \
#--name pcdl_make_gif \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_make_gif output_3d/cell_cell_type_z-5.0 jpeg' \
#--example_input output_3d/cell_cell_type_z-5.0 \
#--example_output output_3d/cell_cell_type_z-5.0/cell_cell_type_z-5.0_jpeg.gif \
#--test_case \
#--help_from_command 'pcdl_make_gif --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

## pcdl_make_movie
#planemo tool_init --force \
#--id pcdl_make_movie \
#--name pcdl_make_movie \
#--version 3.0.0 \
#--requirement pcdl@3.3.6 \
#--example_command 'pcdl_make_movie output_3d/cell_cell_type_z-5.0 jpeg' \
#--example_input output_3d/cell_cell_type_z-5.0 \
#--example_output output_3d/cell_cell_type_z-5.0/cell_cell_type_z-5.0_jpeg12.mp4 \
#--test_case \
#--help_from_command 'pcdl_make_movie --help' \
#--cite_url https://github.com/elmbeech/physicelldataloader

