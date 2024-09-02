#


Official documentation:


### howto load and visualize a  vtk file with the vtk python library
+ https://docs.vtk.org/en/latest/index.html

pip3 install vtk

### howto load and visualize a  vtk file with the vtk python library
+ https://fury.gl/latest/index.html


pip3 install fury


```python
import fury

v_cell = fury.io.load_polydata('output/output00000012_cell.vtk')
actor = fury.get_actor_from_polydata(v_cell)
scene = fury.window.Scene()
showm = fury.window.ShowManager(scene, size=(1024,720), reset_camera=False)
showm.initialize()
scene.add(actor)
showm.start()
```

get_actor_from_polydata(polydata)


Get actor from a vtkPolyData.

get_actor_from_primitive(vertices, triangles, *)


Get actor from a vtkPolyData.



vtk library
fury graph library

https://docs.vtk.org/en/latest/getting_started/index.htmls


```python
import vtk

```
```bash
itk-vtk-viewer
```
+ https://kitware.github.io/itk-vtk-viewer/docs/shortcuts.html






+ https://examples.vtk.org/site/Python/
+ https://examples.vtk.org/site/PythonHowTo/


+ https://docs.vtk.org/en/latest/api/python.html
+ https://examples.vtk.org/site/PythonicAPI/
+ https://examples.vtk.org/site/PythonicAPIComments/
