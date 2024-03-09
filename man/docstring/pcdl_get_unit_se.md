```
usage: pcdl_get_unit_se [-h] [--microenv MICROENV] [--settingxml SETTINGXML]
                        [-v VERBOSE]
                        [path]

function returns a csv that lists all tracked variables from metadata, cell,
and microenvironment and maps them to their unit.

positional arguments:
  path                  path to the PhysiCell output directory. default is . .

options:
  -h, --help            show this help message and exit
  --microenv MICROENV   should the microenvironment be extracted? setting
                        microenv to False will use less memory and speed up
                        processing, similar to the original pyMCDS_cells.py
                        script. default is True.
  --settingxml SETTINGXML
                        from which settings.xml should the cell type ID label
                        mapping be extracted? set to None or False if the xml
                        file is missing! default is PhysiCell_settings.xml.
  -v VERBOSE, --verbose VERBOSE
                        setting verbose to False for less text output, while
                        processing. default is True.

homepage: https://github.com/elmbeech/physicelldataloader
```
