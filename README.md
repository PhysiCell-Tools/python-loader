# pcDataLoader

## Abstract:

pcDataLoader provides a platform independent, python3 based, [pip](https://en.wikipedia.org/wiki/Pip_(package_manager)) installable interface
to load output, generated with the [PhysiCell](https://github.com/MathCancer/PhysiCell) agent based modeling framework,
into [python3](https://en.wikipedia.org/wiki/Python_(programming_language)).

pcDataLoader is forked from the original [PhysiCell-Tools](https://github.com/PhysiCell-Tools) [python-loader](https://github.com/PhysiCell-Tools/python-loader) implementation.

The pcDataLoader python3 library will maintain two main branches:

+ **Branch version 2** will be strictly compatible with the original PhysiCell-Tools/python-loader code, although pip installable.
+ **Branch version 3** might break with old habits, although tries to be as downwards compatible as possible.
  The aim of the v3 branch is to get a very lean and agile physicell output interface for the ones coming from the python3 world to physicell.

Note: there can only be one version of pcDataLoader installed in each python3 environment.


## Header:
+ Language: python >= 3.6
+ Library dependencies: anndata, matplotlib, numpy, pandas
+ Programmer: Patrick Wall, Elmar Bucher, Randy Heiland, Paul Macklin
+ Date of origin original PhysiCell-Tools python-loader: 2019-09-02
+ Date of origin pcDataLoader fork: 2022-08-30
+ License: [BSD-3-Clause](https://en.wikipedia.org/wiki/BSD_licenses)
+ User manual: this README.md file
+ Source code: [https://github.com/elmbeech/pcDataLoader](https://github.com/elmbeech/pcDataLoader)


## HowTo Guide:

+ Check out: man/HOWTO.md(https://github.com/elmbeech/pcDataLoader/man/HOWTO.md)!


## Tutorial:

+ http://www.mathcancer.org/blog/python-loader/
+ Check out: man/TUTORIAL.md(https://github.com/elmbeech/pcDataLoader/man/TUTORIAL.md)!


## Reference:

This is the technical descriptions of the machinery and how to operate it.
References are maintained in each module`s [docstring](https://en.wikipedia.org/wiki/Docstring).

For each pcDataLoader module, get on the fly reference information with the help command.
```python3
import pcDataLoader as pc

help(pc.pyMCDS)
help(pc.pyMCDS_timeseries)
```


## Discussion:
To be developed.


## About Documentation:
Within the pcDataLoader library, I try to stick to the documentation policy lined out by Daniele Procida in his talk "[what nobody tells you about documentation](https://www.youtube.com/watch?v=azf6yzuJt54)" at PyCon 2017 in Portland, Oregon.


## Contributions:
+ original PhysiCell-Tools python-loader implementation: Patrick Wall, Randy Heiland, Paul Macklin
+ fork pcDataLoader implementation: Elmar Bucher


## Release Notes:
+ version 3.0.0 (2022-++-++): elmbeech/pcDataLoader
    + **pyMCDS** takes new additionally a boolean parameter termed **microenv**, to specify if the microenvironment (substrates) should be read out (for completeness) or not (for speed increase).
    + **pyMCDS** parameter **xml_file** can now handle path/file.xml (unix) or path\file.xml (dos) input, aslong output_path is the default.
    + **mcds.get_cell_df** returns now a pandas dataframe with the cell IDs the index and not as a column. additionaly this dataframe has now voxel and voxel_position collumns.
    + **mcds.get_2D_mesh** was renamed to **mcds.get_mesh_2D** for consistency.
    + **mcds.get_cell_variables** and **mcds.get_substrate_names** return now a strictly alphabetically ordered list.
    + new **mcds.get_concentration_df** function.

+ version 2.0.0 (2022-08-30): elmbeech/pcDataLoader pip installable release, derived from and compatible with PhysiCell-Tools/python-loader release 1.1.0 (2022-07-20).
+ version 1.1.0 (2022-05-09): Physicell-Tools/python-loader release compatible with pre-v1.10.x of PhysiCell
+ version 1.0.1 (2020-01-25): Physicell-Tools/python-loader time-series related bug fix
+ version 1.0.0 (2019-09-28): Physicell-Tools/python-loader first public release!
