# PhysiCell Data Loader Tutorial: pcdl and Blender

[Blender](https://www.blender.org/) is a modeling, rigging, animation, simulation, rendering, compositing,  motion tracking, video editing, and game creation software.
Blender is free and open source.
There exists a bioxel nodes plugin, that lets us load ome tiff files.
And there exists a vtk nodes plugin, that lets us load <!-- vtk rectilinear grid data, -->vtk polynomial data files.


## The blender bioxel nodes plugin

Please follow the official bioxelnodes instructions for installation
and to learn how to use the plugin.

+ https://omoolab.github.io/BioxelNodes/latest/


## The blender vtk nodes plugin.

This installation is not for the faint hearted!
Please follow the bvtknodes installation instructions,
including workspace setup.

+ https://bvtknodes.readthedocs.io/en/latest/BVTKNodes.html#installation

### Load vtk polynomial data vtk files

1. In the BVTK Node Tree Workspace, click New to add a Node Tree
2. In the BVTK Node Tree Workspace, right click Add / Reader / vtkXMLPolyDataReader
3. In the BVTK Node Tree Workspace, right click Add / Converters / VTKtoBlenderMesh
4. Connect tkXMLPolyDataReader output with input VTKtoBlenerMesh
5. tkXMLPolyDataReader FileName: path/to/output00000000.vtp
6. VTKtoBlenderMesh: click Update Node

### Load rectilinear grid vtk files

I was not able to brdige that data yet.

<!-- Geometry nodes could be part of the solution,
but I would have to write such a converter.
+ https://www.youtube.com/watch?v=oPenYcM6Usw
+ https://www.youtube.com/watch?v=6LMuT2hN2yw
+ https://www.youtube.com/watch?v=LrEHoaq6QFE
-->

### More about blender and vtk

To learn more about Blender and BVTK Node, please study the official documentation.
+ https://github.com/tkeskita/BVtkNodes/tree/master
+ https://docs.blender.org/manual/en/latest/


That's it!
