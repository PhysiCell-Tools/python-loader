```
usage: pcdl_make_simularium [-h] [--custom_data_type [CUSTOM_DATA_TYPE ...]]
                            [--physiboss PHYSIBOSS] [--settingxml SETTINGXML]
                            [-v VERBOSE] [--trajectory_title TRAJECTORY_TITLE]
                            [--scale_factor [SCALE_FACTOR]]
                            [path] [focus_cat ...]

function returns a simularium trajectory file which can be viewed with the
simularium viewer (https://simularium.allencell.org/).

positional arguments:
  path                  path to the PhysiCell output directory or a
                        outputnnnnnnnn.xml file. default is . .
  focus_cat             one or two categorical mcds.get_cell_df dataframe
                        column names, used for cell attributes. default is
                        cell_type current_phase.

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
                        the settings.xml that is loaded, from which the cell
                        type ID label mapping, is extracted, if this
                        information is not found in the output xml file. set
                        to None or False if the xml file is missing! default
                        is False.
  -v, --verbose VERBOSE
                        setting verbose to False for less text output, while
                        processing. default is True.
  --trajectory_title, --tt TRAJECTORY_TITLE
                        the trajectory_title will be used for
                        <trajectory_title>.simularium file name and displayed
                        in the simulation.
  --scale_factor, --sf [SCALE_FACTOR]
                        a multiplier for the scene, use if visualization is
                        too large or small. if none is provided, one will be
                        calculated based on the position data. bue 20260822:
                        does not seem to work in current simularium 1.13.0.

homepage: https://github.com/elmbeech/physicelldataloader
```
