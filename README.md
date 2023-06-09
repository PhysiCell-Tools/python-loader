![physicellcdataloader logo & title](man/img/pcdataloader_title_v3.0.0.png)

## Abstract:
physicelldataloader (pcdl) provides a platform independent, python3 based, [pip](https://en.wikipedia.org/wiki/Pip_(package_manager)) installable interface
to load output, generated with the [PhysiCell](https://github.com/MathCancer/PhysiCell) agent based modeling framework,
into [python3](https://en.wikipedia.org/wiki/Python_(programming_language)).

pcdl was forked from the original [PhysiCell-Tools](https://github.com/PhysiCell-Tools) [python-loader](https://github.com/PhysiCell-Tools/python-loader) implementation.

The pcdl python3 library will maintain two main branches:

+ **Branch version 2** will be strictly compatible with the original PhysiCell-Tools/python-loader code, although pip installable.
+ **Branch version 3** might break with old habits, although tries to be as much downward compatible as possible.
  The aim of the v3 branch is to get a very lean and agile python3 physicell output interface, for the ones coming from the python world.

Note: there can only be one version of pcdl installed in each python3 environment.
In the HowTo guide is in detail explained, how to install and uninstall pcdl branches.


## Header:
+ Language: python >= 3.6
+ Library dependencies: anndata, matplotlib, numpy, pandas
+ Date of origin original PhysiCell-Tools python-loader: 2019-09-02
+ Date of origin pcdl fork: 2022-08-30
+ License: [BSD-3-Clause](https://en.wikipedia.org/wiki/BSD_licenses)
+ User manual: this README.md file
+ Source code: [https://github.com/elmbeech/physicelldataloader](https://github.com/elmbeech/physicelldataloader)


## HowTo Guide:
+ Check out: [man/HOWTO.md](https://github.com/elmbeech/physicelldataloader/tree/master/man/HOWTO.md)!


## Tutorial:
+ Check out: [man/TUTORIAL.md](https://github.com/elmbeech/physicelldataloader/tree/master/man/TUTORIAL.md)!


## Reference Manual:
+ Check out: [man/REFERENCE.md](https://github.com/elmbeech/physicelldataloader/tree/master/man/REFERENCE.md)!


## Discussion:
To be developed.


## About Documentation:
Within the pcdl library, we tried to stick to the documentation policy lined out by Daniele Procida in his "[what nobody tells you about documentation](https://www.youtube.com/watch?v=azf6yzuJt54)" talk at PyCon 2017 in Portland, Oregon.


## Contributions:
+ original PhysiCell-Tools python-loader implementation: Patrick Wall, Randy Heiland, Paul Macklin
+ fork pcdl implementation: Elmar Bucher
+ fork pcdl co-programmer and tester: Heber Rocha, Marshal Gress


## Release Notes:
+ version 3.1.0: elmbeech/physicelldataloader
    + resolve pypa problem.
    + change the pcDataLoader library name to pcdl.

+ version 3.0.7 (2023-06-08): elmbeech/physicelldataloader
    + pyMCDSts: replaces the svg dependent **mcdsts.make_jpeg**, **mcdsts.make_png**, and **mcdsts.make_tiff** with **mcdsts.make_imgcell** and **mcdsts.make_imgsubs** which generate images straight out of the loaded data. the **mcdsts.make_gif** and **mcdsts.make_movie** functions were adjusted accordingly.
    + pyMCDSts: **mcdsts.read_mcds** loads now automatically all mcds snapshots, if no xmlfile_list is provided (default).

+ version 3.0.6 (2023-04-29): elmbeech/physicelldataloader
    + pyMCDS **\_read_xml** is now able to load time steps with zero cells.
    + pyMCDS **mcds.get_contour** can handle more input parameters.

+ version 3.0.5 (2023-02-26): elmbeech/physicelldataloader pyMCDS **mcds.get_contour**  plots span now the whole domain and not only to the border voxel centers.
+ version 3.0.4 (2023-02-21): elmbeech/physicelldataloader pyMCDS **mcds.get_contour** function, to easily generate for substrates matplotlib contourf and contour plots because they do not exist as pandas plots.
+ version 3.0.3 (2023-02-19): elmbeech/physicelldataloader branch 3 has no longer anndata and as such hdf5 dependency.
+ version 3.0.2 (2023-01-06): elmbeech/physicelldataloader bugfix installing package data.
+ version 3.0.0 (2023-01-06): elmbeech/physicelldataloader
    + **pyMCDS** parameter **xml_file** can now handle path/file.xml (unix) or path\file.xml (dos) input, as long output_path is the default.
    + **pyMCDS** has a new additional boolean **microenv** parameter, to specify if the microenvironment (substrates) should be read (for completeness) or not (for speed increase and less memory usage).
    + **pyMCDS** has a new additional boolean **graph** parameter, to specify if the attached and neighbor graph should be read.
    + **pyMCDS** has a new additional boolean **verbose** parameter, to specify if there should be text output while processing.
    + pyMCDS **mcds.get_2D_mesh** was renamed to **mcds.get_mesh_2D** for consistency.
    + pyMCDS **mcds.get_linear_voxels** was renamed to **mcds.get_mesh_coordinate** for consistency.
    + pyMCDS **mcds.get_containing_voxel_ijk** was renamed to **mcds.get_voxel_ijk** for briefness.
    + pyMCDS **mcds.get_voxel_spacing** returns now 3 specific values, one for x, y, and z, instead of 1 general value.
    + pyMCDS **mcds.get_concentrations** was renamed to **mcds.get_concentration** for consistency
    + pyMCDS **mcds.get_concentrations_at** was renamed to **mcds.get_concentration_at** for consistency
    + pyMCDS **mcds.get_concentration_at** if z_slice is not a mesh center value, the function will by default adjust to nearest and no longer break.
    + pyMCDS **mcds.get_cell_variables** and **mcds.get_substrate_names** return now a strictly alphabetically ordered list.
    + pyMCDS **mcds.get_cell_df** returns now a pandas dataframe with the cell IDs as the index and not as a column.
      additionally, this dataframe contains now voxel, mesh_center, substrate parameter, substrate concentration, and cell density information too.
    + new pyMCDS **mcds.get_concentration_df** function.
    + new pyMCDS **mcds.get_substrate_df** function.
    + new pyMCDS **mcds.get_unit_df** function.
    + new pyMCDS **mcds.get_multicellds_version** function.
    + new pyMCDS **mcds.get_physicell_version** function.
    + new pyMCDS **mcds.get_runtime** function.
    + new pyMCDS **mcds.get_timestamp** function.
    + new pyMCDS **mcds.get_voxel_ijk_range** function.
    + new pyMCDS **mcds.get_voxel_ijk_axis** function.
    + new pyMCDS **mcds.get_voxel_spacing** function.
    + new pyMCDS **mcds.get_voxel_volume** function.
    + new pyMCDS **mcds.get_mesh_mnp_range** function.
    + new pyMCDS **mcds.get_mesh_mnp_axis** function.
    + new pyMCDS **mcds.get_xyz_range** function.
    + new pyMCDS **mcds.is_in_mesh** function.
    + new pyMCDS **mcds.get_attached_graph_dict** function.
    + new pyMCDS **mcds.get_neigbor_graph_dict** function.
    + class **pyMCDS_timeseries** was renamed to **pyMCDSts** and completely rewritten.
    + new pyMCDSts **get_xmlfile_list** function.
    + new pyMCDSts **read_mcds** function.
    + new pyMCDSts **make_jpeg** function.
    + new pyMCDSts **make_png** function.
    + new pyMCDSts **make_tiff** function.
    + new pyMCDSts **make_gif** function.
    + new pyMCDSts **make_movie** function.
    + all **plotting** functions were removed because pcdl only focus on making the raw data in python easy accessible for in-depth analysis.
    + cell position coordinates are now constantly labeled as **x,y,z**, mesh center coordinates as **m,n,p**, and voxel coordinates as **i,j,k**.
    + the underling [mcds object data dictionary structure](https://github.com/elmbeech/physicelldataloader/tree/master/man/img/physicelldataloader_data_dictionary_v3.0.0.png) has changed.
    + [pytest](https://en.wikipedia.org/wiki/Pytest) unit tests exist now for all pyMCDS and pyMCDSts functions.

+ version 2.0.2 (2023-01-06): elmbeech/physicelldataloader reset patch voxel spacing bugfix, so that branch2 is full compatible with branch1 again. use branch3 for a bugfixed version!
+ version 2.0.1 (2022-11-08): elmbeech/physicelldataloader beta release patch voxel spacing bugfix.
+ version 2.0.0 (2022-08-30): elmbeech/physicelldataloader pip installable release, derived from and compatible with PhysiCell-Tools/python-loader release 1.1.0 (2022-07-20).
+ version 1.1.0 (2022-05-09): Physicell-Tools/python-loader release compatible with pre-v1.10.x of PhysiCell
+ version 1.0.1 (2020-01-25): Physicell-Tools/python-loader time-series related bug fix
+ version 1.0.0 (2019-09-28): Physicell-Tools/python-loader first public release!
