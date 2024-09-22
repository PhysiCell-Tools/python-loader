# mcds.is_in_mesh()


## input:
```
            x: floating point number
                position x-coordinate.

            y: floating point number
                position y-coordinate.

            z: floating point number
                position z-coordinate.

            halt: boolean; default is False
                should program execution break or just spit out a warning,
                if position is not in mesh?

```

## output:
```
            b_isinmesh: boolean
                declares if the given coordinate is inside the mesh.

```

## description:
```
            function evaluates if the given position coordinate
            is inside the boundaries. if the coordinate is outside the
            mesh, a warning will be printed. if additionally
            halt is set to True, program execution will halt.
        
```