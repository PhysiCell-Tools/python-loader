# PhysiCell Data Loader Tutorial: pcdl and the Matlab programming language  <!-- and Octave-->

<!-- bue 20240903: could someone test and update who actuallty knows matlab or octave? -->

[Matlab](https://www.mathworks.com/products/matlab.html)  <!-- / GNU [Octave](https://octave.org/) -->
is a scientific programming language used by some engineers.
The earliest predecessor of pcdl was actually a Matlab implementation.
Coming full circle, this is how you can load some pcdl data constructs into Matlab.  <!-- and GNU Octave -->


## &#x2728; Handle csv files

### Save pcdl data constructs as csv files from the command line

```bash
cd path/to/PhysiCell
```
```bash
pcdl_get_conc_df output
```
```bash
pcdl_get_cell_df output
```

### Load a csv file into Matlab as a tabel  <!-- or Octave -->

&#x26A0; **bue 2024-09-22:** the readtable function is not yet implemented in Octave.

+ https://www.mathworks.com/help/matlab/ref/readtable.html
+ https://www.mathworks.com/help/matlab/matlab_external/python-pandas-dataframes.html

```matlab
df_conc = readtable("output/timeseries_conc.csv")
```
```matlab
df_cell = readtable("output/timeseries_conc.csv")
```


## &#x2728; Handle json files

### Save pcdl data constructs as json files from the command line

```bash
cd path/to/PhysiCell
```
```bash
pcdl_get_conc_attribute output 2
```
```bash
pcdl_get_cell_attribute output 2
```

### Load json files into Matlab  <!-- or Octave -->

```matlab
struct_conc = jsondecode(fileread("output/timeseries_conc_attribute_minmax.json"))
```
```matlab
struct_conc = jsondecode(fileread("output/timeseries_conc_attribute_minmax.json"))
```


## &#x2728; Handle gml graph data files

### Save pcdl data constructs as gml files from the command line

```bash
cd path/to/PhysiCell
```
```bash
pcdl_make_graph_gml output/output00000024.xml neighbor --node_attribute cell_type dead oxygen pressure
```

### Load gml files into a Matlab data construct  <!-- or Octave -->

We will use the [matlab-igraph](https://www.mathworks.com/matlabcentral/fileexchange/159001-matlab-igraph) toolbox to load gml files into Matlab.  <!-- or Octave -->
+ https://github.com/DavidRConnell/matlab-igraph
+ https://igraph.org/

Install matlab-igaph toolbox.

Load gml file.

```matlab
"output/output00000024_neighbor.gml"
```


## &#x2728; Handle vtk files

### Save pcdl data constructs as vtr rectilinear grid files and vtp polynomial data file from the command line

```bash
cd path/to/PhysiCell
```
```bash
pcdl_make_conc_vtk output/output00000024.xml
```
```bash
pcdl_make_cell_vtk output/output00000024.xml
```

### Load a vtr and vtp file into Matlab data construct  <!-- or Octave -->

+ https://www.mathworks.com/matlabcentral/fileexchange/94993-vtktoolbox
+ https://github.com/KIT-IBT/vtkToolbox

Install vtk toolbox.

Load vtkfile into Matlab.  <!-- or Octave -->

```matlab
output/output00000024_conc.vtr
```
```matlab
output/output00000024_cell.vtp
```
