# PhysiCell Data Loader Tutorial: pcdl and Python and Vtk

[Vtk](https://vtk.org/) the visualization tool kit is a software for image processing, 3D graphics, volume rendering and visualization. \
Naturally, we will work here with the [vtk](https://pypi.org/project/vtk/) python binding.
Additionally, we will have a look at [fury](https://fury.gl/latest/index.html),
a python library for scientific visualization and 3D animations that can read vtk files too.



### How to load and visualize a vtk file with the vtk python library

Installation.

```bash
pip3 install vtk
```

Load and display conc.vtr rectilinear grid files.

<!-- Randy or Furkan, could you help me to write a minimal example? -->
```python3
import vtk

# load file
vtr_reader = vtk.vtkXMLRectilinearGridReader()
vtr_reader.SetFileName('output/output00000012_conc.vtr')
vtr_reader.Update()

# actor
# scene
# show
```

Load and display cell.vtp ploydata files.

<!-- Randy or Furkan, could you help me to write a minimal example? -->
```python3
import vtk

# load file
vtp_reader = vtk.vtkXMLPolyDataReader()
vtp_reader.SetFileName('output/output00000012_cell.vtp')
vtp_reader.Update()

# actor
# scene
# show
```

Official documentation:
+ https://book.vtk.org/en/latest/index.html
+ https://docs.vtk.org/en/latest/index.html
+ https://examples.vtk.org/site/Python/
+ https://examples.vtk.org/site/PythonHowTo/
+ https://examples.vtk.org/site/PythonicAPI/
+ https://examples.vtk.org/site/PythonicAPIComments/


### How to load and visualize a vtk file with the fury python library

Installation.

```bash
pip3 install fury
```

<!-- i have ask Elef if this is true -->
Unfortunately, rectiliniar grid files are not supported, conc.vtr files cannot be loaded.

Load and display cell.vtp ploydata files.

```python
import fury

# load file
v_cell = fury.io.load_polydata('output/output00000012_cell.vtp')

# actor
actor = fury.get_actor_from_polydata(v_cell)

# scene
scene = fury.window.Scene()
scene.add(actor)

# show
showm = fury.window.ShowManager(scene)
showm.initialize()
showm.start()
```

Official documentation:
+ https://fury.gl/latest/index.html

