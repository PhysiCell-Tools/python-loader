```
usage: pcdl_make_cell_vtk [-h] [--custom_data_type [CUSTOM_DATA_TYPE ...]]
                          [--physiboss PHYSIBOSS] [--settingxml SETTINGXML]
                          [-v VERBOSE] [--attribute [ATTRIBUTE ...]]
                          [path]

function that generates 3D glyph vtk file for cells. cells can have specified
attributes like cell_type, pressure, dead, etc. you can post-process this file
in other software like paraview (https://www.paraview.org/).

positional arguments:
  path                  path to the PhysiCell output directory or a
                        outputnnnnnnnn.xml file. default is . .

options:
  -h, --help            show this help message and exit
  --custom_data_type [CUSTOM_DATA_TYPE ...]
                        parameter to specify custom_data variable types other
                        than float (namely: int, bool, str) like this
                        var:dtype myint:int mybool:bool mystr:str . downstream
                        float and int will be handled as numeric, bool as
                        Boolean, and str as categorical data. default is an
                        empty string.
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
  --attribute [ATTRIBUTE ...]
                        listing of mcds.get_cell_df dataframe column names,
                        used for cell attributes. default is a single term:
                        cell_type.

homepage: https://github.com/elmbeech/physicelldataloader
```