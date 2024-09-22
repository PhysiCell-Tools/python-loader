```
usage: pcdl_get_conc_df [-h] [-v VERBOSE] [--drop [DROP ...]]
                        [--keep [KEEP ...]] [--collapse COLLAPSE]
                        [path] [values]

this function extracts dataframes with concentration values for all chemical
species in all voxels and saves them as csv files. additionally, this
dataframe lists voxel and mesh center coordinates.

positional arguments:
  path                  path to the PhysiCell output directory or a
                        outputnnnnnnnn.xml file. default is . .
  values                minimal number of values a variable has to have in any
                        of the mcds time steps to be outputted. variables that
                        have only 1 state carry no information. None is a
                        state too. default is 1.

options:
  -h, --help            show this help message and exit
  -v VERBOSE, --verbose VERBOSE
                        setting verbose to False for less text output, while
                        processing. default is True.
  --drop [DROP ...]     set of column labels to be dropped for the dataframe.
                        don't worry: essential columns like ID, coordinates
                        and time will never be dropped. Attention: when the
                        keep parameter is given, then the drop parameter has
                        to be an empty string! default is an empty string.
  --keep [KEEP ...]     set of column labels to be kept in the dataframe. set
                        values=1 to be sure that all variables are kept. don't
                        worry: essential columns like ID, coordinates and time
                        will always be kept. default is an empty string.
  --collapse COLLAPSE   should all mcds time steps from the time series be
                        collapsed into one big csv, or a many csv, one for
                        each time step? default is True.

homepage: https://github.com/elmbeech/physicelldataloader
```
