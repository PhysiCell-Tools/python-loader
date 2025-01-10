# mcds.get_voxel_ijk()


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
            li_ijk : list of 3 integers
                i, j, k indices for the voxel
                containing the x, y, z position.

```

## description:
```
            function returns the meshgrid indices i, j, k
            for the given position x, y, z.
        
```