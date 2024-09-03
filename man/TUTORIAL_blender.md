# PhysiCell Data Loader Tutorial: pcdl and Blender

[Blender](https://www.blender.org/) is modeling, rigging, animation, simulation, rendering, compositing,  motion tracking, video editing, and game creation software tool.
Blender is free and open source, and there exist a vtk blener plugin,
that let us load <!-- vtk rectilinear grid data, -->vtk polynom data <!--, and even ome.tiff data -->.


## Install bender and the blender vtk node plugin.

This installation is not for the faint hearted!
Please follow the bvtknodes installation instruction,
including workspace setup.

+ https://bvtknodes.readthedocs.io/en/latest/BVTKNodes.html#installation


## Load vtk polynom data vtk files

1. In the BVTK Node Tree Workspace, click New to add a Node Tree
2. In the BVTK Node Tree Workspace, right click Add / Reader / vtkXMLPolyDataReader
3. In the BVTK Node Tree Workspace, right click Add / Converters / VTKtoBlenderMesh
4. Connect tkXMLPolyDataReader output with input VTKtoBlenerMesh
5. tkXMLPolyDataReader FileName: path/to/output00000000.vtp
6. VTKtoBlenderMesh: click Update Node

That's it!


## Load rectilinear grid vtk files

I was not able to brdige that data yet.

<!-- Geometry nodes could be part of the solution,
but I would have to write such a converter.
+ https://www.youtube.com/watch?v=oPenYcM6Usw
+ https://www.youtube.com/watch?v=6LMuT2hN2yw
+ https://www.youtube.com/watch?v=LrEHoaq6QFE
-->

<!-- ## Load ome.tiff vtk files -->


# PhysiCell Data Loader Blender Tutorial

To learn more about Blender and BVTK Node, please study the official documentation.
+ https://github.com/tkeskita/BVtkNodes/tree/master
+ https://docs.blender.org/manual/en/latest/

