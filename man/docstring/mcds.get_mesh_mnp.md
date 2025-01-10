# mcds.get_mesh_mnp()


## input:
```
            x: floating point number
                position x-coordinate.

            y: floating point number
                position y-coordinate.

            z: floating point number
                position z-coordinate.

            is_in_mesh: boolean; default is True
                should function check, if the given coordinate is in the mesh,
                and only calculate ijk values if is so?

```

## output:
```
            lr_mnp : list of 3 floats
                m, n, p indices for the mesh center,
                for the mesh cell containing the x, y, z position.

```

## description:
```
            function returns the meshgrid indices m, n, p
            for the given position x, y, z.
        
```