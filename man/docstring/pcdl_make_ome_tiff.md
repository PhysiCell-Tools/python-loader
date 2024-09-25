```
usage: pcdl_make_ome_tiff [-h] [--microenv MICROENV] [--physiboss PHYSIBOSS]
                          [--settingxml SETTINGXML] [-v VERBOSE]
                          [--collapse COLLAPSE]
                          [path] [cell_attribute]

function to transform chosen mcds output into an 1[um] spaced czyx (channel,
z-axis, y-axis, x-axis) ome tiff file, one substrate or cell_type per channel.
the ome tiff file format can for example be read by the napari
(https://napari.org/stable/) or fiji imagej (https://fiji.sc/) software.

positional arguments:
  path                  path to the PhysiCell output directory or a
                        outputnnnnnnnn.xml file. default is . .
  cell_attribute        mcds.get_cell_df dataframe columns, used for
                        cell_attribute. the column data type has to be numeric
                        (bool, int, float) and can not be string. default is
                        ID, with will result in a segmentation mask.

options:
  -h, --help            show this help message and exit
  --microenv MICROENV   should the microenvironment data be loaded? setting
                        microenv to False will use less memory and speed up
                        processing, similar to the original pyMCDS_cells.py
                        script. default is True.
  --physiboss PHYSIBOSS
                        if found, should physiboss state data be extracted and
                        loaded into the df_cell dataframe? default is True.
  --settingxml SETTINGXML
                        the settings.xml that is loaded, from which the cell
                        type ID label mapping, is extracted, if this
                        information is not found in the output xml file. set
                        to None or False if the xml file is missing! default
                        is PhysiCell_settings.xml.
  -v VERBOSE, --verbose VERBOSE
                        setting verbose to False for less text output, while
                        processing. default is True.
  --collapse COLLAPSE   should all mcds time steps from the time series be
                        collapsed into one big ome.tiff, or a many ome.tiff,
                        one ome.tiff for each time step?, default is True.

homepage: https://github.com/elmbeech/physicelldataloader
```
