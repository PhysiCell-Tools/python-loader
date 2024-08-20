```
usage: pcdl_make_conc_vtk [-h] [--custom_data_type [CUSTOM_DATA_TYPE ...]]
                          [--microenv MICROENV] [--physiboss PHYSIBOSS]
                          [--settingxml SETTINGXML] [-v VERBOSE]
                          [path]

function generates rectilinear grid vtk files, one per mcds time step,
contains distribution of substrates over microenvironment. you can post-
process this files in other software like paraview
(https://www.paraview.org/).

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

homepage: https://github.com/elmbeech/physicelldataloader
```
