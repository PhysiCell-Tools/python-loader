# mcds.make_cell_vtk()


## input:
```

            attributes: list of strings; default is 'cell_type'
                column name within cell dataframe.

            visualize: boolean; default is True
                visualize cells using vtk Renderer.

```

## output:
```
            vtk: 3D Glyph VTK that contains cells.


```

## description:
```
            function that 3D Glyph VTK file for cells. Cells can have
            specificed attributes like 'cell_type', 'pressure', 'dead', etc.
            You can post-process this file in other software like Paraview.
        
```