```
usage: pcdl_make_ome_tiff [-h] [--microenv MICROENV] [--physiboss PHYSIBOSS]
                          [--settingxml SETTINGXML] [-v VERBOSE]
                          [--cell_attribute CELL_ATTRIBUTE]
                          [--collapse COLLAPSE]
                          [path]

function to transform chosen mcds output into an 1[um] spaced czyx (channel,
z-axis, y-axis, x-axis) ome tiff file, one substrate or cell_type per channel.
the ome tiff file format can for example be read by the napari
(https://napari.org/stable/) or fiji imagej (https://fiji.sc/) software.

positional arguments:
  path                  path to the PhysiCell output directory or a
                        outputnnnnnnnn.xml file. default is . .

options:
  -h, --help            show this help message and exit
  --microenv MICROENV   should the microenvironment be extracted? setting
                        microenv to False will use less memory and speed up
                        processing, similar to the original pyMCDS_cells.py
                        script. default is True.
  --physiboss PHYSIBOSS
                        if found, should physiboss state data be extracted and
                        loaded into the df_cell dataframe? default is True.
  --settingxml SETTINGXML
                        from which settings.xml should the cell type ID label
                        mapping be extracted? set to None or False if the xml
                        file is missing! default is PhysiCell_settings.xml.
  -v VERBOSE, --verbose VERBOSE
                        setting verbose to False for less text output, while
                        processing. default is True.
  --cell_attribute CELL_ATTRIBUTE
                        mcds.get_cell_df dataframe columns, used for
                        cell_attribute. the column data type has to be numeric
                        (bool, int, float) and can not be string. default is
                        ID, with will result in a segmentation mask.
  --collapse COLLAPSE   should all mcds time steps from the time series be
                        collapsed into one big ome.tiff, or a many ome.tiff,
                        one ome.tiff for each time step?, default is True.

homepage: https://github.com/elmbeech/physicelldataloader
```
