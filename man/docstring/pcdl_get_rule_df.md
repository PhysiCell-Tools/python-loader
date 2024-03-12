```
usage: pcdl_get_rule_df [-h] [--settingxml SETTINGXML] [-v VERBOSE] [path]

function returns a csv that lists all tracked variables from metadata, cell,
and microenvironment and maps them to their unit.

positional arguments:
  path                  path to the PhysiCell output directory. default is . .

options:
  -h, --help            show this help message and exit
  --settingxml SETTINGXML
                        from which settings.xml should the rule.csv links be
                        extracted? default is PhysiCell_settings.xml.
  -v VERBOSE, --verbose VERBOSE
                        setting verbose to False for less text output, while
                        processing. default is True.

homepage: https://github.com/elmbeech/physicelldataloader
```
