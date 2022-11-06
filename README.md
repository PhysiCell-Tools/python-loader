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
+ Check out: [man/HOWTO.md](https://github.com/elmbeech/pcDataLoader/tree/master/man/HOWTO.md)!


## Tutorial:
+ http://www.mathcancer.org/blog/python-loader/
+ Check out: [man/TUTORIAL.md](https://github.com/elmbeech/pcDataLoader/tree/master/man/TUTORIAL.md)!


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
    + **pyMCDS** parameter **xml_file** can now handle path/file.xml (unix) or path\file.xml (dos) input, aslong output_path is the default.
    + **pyMCDS** has a new additionally a boolean **microenv** parameter, to specify if the microenvironment (substrates) should be read (for completeness) or not (for speed increase and less memory usage).
    + **pyMCDS** has a new additionally a boolean **graph** parameter, to specify if the attached and neighbor graph should be read.
    + **pyMCDS** has a new additionally a boolean **verbose** parameter, to specify if there should be text output while processing.
    + pyMCDS **mcds.get_2D_mesh** was renamed to **mcds.get_mesh_2D** for consistency.
    + pyMCDS **mcds.get_linear_voxels** was renamed to **mcds.get_mesh_axis** for consistency.
    + pyMCDS **mcds.get_containing_voxel_ijk** was renamed to **mcds.get_voxel_ijk** for biefness.
    + pyMCDS **mcds.get_voxel_spacing** returns now 3 specific values, one for x, y, and z, insted of 1 general value.
    + pyMCDS **mcds.get_concentrations** if z_slice is not a mesh center value, the function will by default adjust to nearest and no longer break.
    + pyMCDS **mcds.get_cell_variables** and **mcds.get_substrate_names** return now a strictly alphabetically ordered list.
    + pyMCDS **mcds.get_cell_df** returns now a pandas dataframe with the cell IDs the index and not as a column.
      additionaly, this dataframe has now voxel, mesh_center, substrate parameter, substrate concentration, and cell density columns.
    + new pyMCDS **mcds.get_concentration_df** function.
    + new pyMCDS **mcds.get_substrate_df** function.
    + new pyMCDS **mcds.get_unit_df** function.
    + new pyMCDS **mcds.get_physicell_version** function.
    + new pyMCDS **mcds.get_runtime** function.
    + new pyMCDS **mcds.get_timestamp** function.
    + new pyMCDS **mcds.get_voxel_volume** function.
    + new pyMCDS **mcds.get_x_range** function.
    + new pyMCDS **mcds.get_y_range** function.
    + new pyMCDS **mcds.get_z_range** function.
    + new pyMCDS **mcds.get_mesh_m_range** function.
    + new pyMCDS **mcds.get_mesh_n_range** function.
    + new pyMCDS **mcds.get_mesh_p_range** function.
    + new pyMCDS **mcds.get_voxel_i_range** function.
    + new pyMCDS **mcds.get_voxel_j_range** function.
    + new pyMCDS **mcds.get_voxel_k_range** function.
    + new pyMCDS **mcds.is_in_mesh** function.
    + new pyMCDS **mcds.get_attached_graph_dict** function.
    + new pyMCDS **mcds.get_neigbor_graph_dict** function.
    + class **pyMCDS_timeseries** was renamed to **pyMCDSts** and completly rewritten.
    + new pyMCDSts **get_xmlfile_list** function.
    + new pyMCDSts **read_mcds** function.
    + new pyMCDSts **make_gif** function.
    + new pyMCDSts **make_jpeg** function.
    + new pyMCDSts **make_png** function.
    + new pyMCDSts **make_tiff** function.
    + new pyMCDSts **make_movie** function.
    + all **plotting** functions were removed, because pcDataLoader only focus on making the raw data in python easy accessible for in-depth analysis.
    + cell positon coordinats are now constandly labeld as **x,y,z**, mesh center coordinates as **m,n,p**, and voxel coordinates as **i,j,k**.
    + the underling mcds object dictionary structure has changed.
    + [pytest](https://en.wikipedia.org/wiki/Pytest) unit tests exist now for all pyMCDS and pyMCDSts functions.

+ version 2.0.0 (2022-08-30): elmbeech/pcDataLoader pip installable release, derived from and compatible with PhysiCell-Tools/python-loader release 1.1.0 (2022-07-20).
+ version 1.1.0 (2022-05-09): Physicell-Tools/python-loader release compatible with pre-v1.10.x of PhysiCell
+ version 1.0.1 (2020-01-25): Physicell-Tools/python-loader time-series related bug fix
+ version 1.0.0 (2019-09-28): Physicell-Tools/python-loader first public release!
