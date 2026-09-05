```
usage: pcdl_get_muspan [-h] [--custom_data_type [CUSTOM_DATA_TYPE ...]]
                       [--microenv MICROENV] [--graph GRAPH]
                       [--physiboss PHYSIBOSS] [--settingxml SETTINGXML]
                       [-v VERBOSE] [--values VALUES] [--drop [DROP ...]]
                       [--keep [KEEP ...]]
                       [path] [z_slice]

function to transform mcds time steps into muspan domain objects for
downstream analysis.

positional arguments:
  path                  path to the PhysiCell output directory or a
                        outputnnnnnnnn.xml file. default is . .
  z_slice               z-axis position to slice a 2D xy-plain out of the 3D
                        mesh. if z_slice position numeric but not an exact
                        mesh center coordinate, then z_slice will be adjusted
                        to the nearest mesh center value, the smaller one, if
                        the coordinate lies on a saddle point. default is 0.0.

options:
  -h, --help            show this help message and exit
  --custom_data_type [CUSTOM_DATA_TYPE ...]
                        parameter to specify custom_data variable types other
                        than float (namely: int, bool, str) like this
                        var:dtype myint:int mybool:bool mystr:str . downstream
                        float and int will be handled as numeric, bool as
                        Boolean, and str as categorical data. default is an
                        empty string.
  --microenv MICROENV   should the microenvironment be extracted and loaded
                        into the muspan domain object? setting microenv to
                        False will use less memory and speed up processing.
                        default is True.
  --graph GRAPH         should neighbor graph, attach graph, and attached
                        spring graph be extracted and loaded into the muspan
                        domain object? default is True.
  --physiboss PHYSIBOSS
                        if found, should physiboss state data be extracted and
                        loaded into the muspan domain object? default is True.
  --settingxml SETTINGXML
                        the settings.xml that is loaded, from which the cell
                        type ID label mapping, is extracted, if this
                        information is not found in the output xml file. set
                        to None or False if the xml file is missing! default
                        is False.
  -v, --verbose VERBOSE
                        setting verbose to False for less text output, while
                        processing. default is True.
  --values VALUES       minimal number of values a variable has to have in any
                        of the mcds time steps to be outputted. variables that
                        have only 1 state carry no information. None is a
                        state too. default is 1.
  --drop [DROP ...]     set of column labels to be dropped for the dataframe.
                        don't worry: essential columns like ID, coordinates
                        and time will never be dropped. Attention: when the
                        keep parameter is given, then the drop parameter has
                        to be an empty string! default is an empty string.
  --keep [KEEP ...]     set of column labels to be kept in the dataframe. set
                        values=1 to be sure that all variables are kept. don't
                        worry: essential columns like ID, coordinates and time
                        will always be kept. default is an empty string.

homepage: https://github.com/elmbeech/physicelldataloader
```
