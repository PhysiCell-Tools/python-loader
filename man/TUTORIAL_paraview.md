# PhysiCell Data Loader Tutorial: pcdl and Paraview

[Paraview](https://www.paraview.org/) is a free and open source scientific visualization software,
that lets us load and analyze vtk rectilinear grid data and vtk polynomial data.

## Install paraview

For installation, please follow the official documentation.
+ https://www.paraview.org/download/


## Generate vtk files from the command line

```bash
pcdl_make_conc_vtk output
```
```bash
pcdl_make_cell_vtk output
```


## Generate vtk files from within python

```python
import pcdl

mcdsts = pcdl.TimeSeries('output/')
mcdsts.make_conc_vtk()
mcdsts.make_cell_vtk()
```


## Load vtk files with paraview

1. **File** / **Open...** path/to/PhysiCell/output/output..\_conc.vtr [OK]
2. **File** / **Open...** path/to/PhysiCell/output/output..\_cell.vtp [OK]
3. In the **Pipeline Browser**, click the closed **eyes**, so that they open.
4. In the **Pipeline Browser**, click the output00000000\_cell.vtp\* and under **Coloring** select cell\_type in the dropdown menu.
5. In the **Pipeline Browser**, click the output00000000\_conc.vtr\* and under **Coloring** select the oxygen in the dropdown menu.

For learning more about how to run the software,
please work through the official documentation.
+ https://www.paraview.org/resources/


That's it. The rest is analysis within paraview!
