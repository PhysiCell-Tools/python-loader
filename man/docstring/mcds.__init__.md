# mcds = pcdl.TimeStep('path/to/outputnnnnnnnn.xml')


## input:
```
            xmlfile: string
                name of the xml file with or without path.
                in the with path case, output_path has to be set to the default!

            output_path: string; default '.'
                relative or absolute path to the directory where
                the PhysiCell output files are stored.

            custom_type: dictionary; default is {}
                variable to specify custom_data variable types
                besides float (int, bool, str) like this: {var: dtype, ...}.
                downstream float and int will be handled as numeric,
                bool as Boolean, and str as categorical data.

            microenv: boole; default True
                should the microenvironment be extracted?
                setting microenv to False will use less memory and speed up
                processing, similar to the original pyMCDS_cells.py script.

            graph: boole; default True
                should neighbor garph and attached graph be extracted?
                setting graph to False will use less memory and speed up processing.

            settingxml: string; default PhysiCell_settings.xml
                from which settings.xml should the cell type ID label mapping
                be extracted?
                set to None or False if the xml file is missing!

            verbose: boole; default True
                setting verbose to False for less text output while processing.

```

## output:
```
            mcds: TimeStep class instance
                all fetched content is stored at mcds.data.

```

## description:
```
            TimeStep.__init__ will call pyMCDS.__init__ that generates a mcds
            class instance, a dictionary of dictionaries data structure that
            contains all output from a single PhysiCell model time step.
            furthermore, the mcds object offers functions to access the stored data.
            the code assumes that all related output files are stored
            in the same directory. data is loaded by reading the xml file for
            a particular time step and the therein referenced files.
        
```