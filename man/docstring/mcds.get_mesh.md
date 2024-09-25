# mcds.get_mesh()


## input:
```
            flat : bool; default False
                if flat is True, only the m-axis mesh center
                and n-axis mesh center meshgrids will be returned.
                else the m, n, and p mesh center meshgrids will be returned.

```

## output:
```
            aar_meshgrid : 4-way (3D) or 3-way (2D) numpy arrays tensor of floating point numbers
                meshgrid shaped object, with the mesh center
                coordinate values from the m, n, p-axis or m, n-axis.

```

## description:
```
            function returns a numpy array of meshgrids each of which stores
            the mesh center coordinate values from one particular axis.
            the function can either return meshgrids for the full
            m, n, p 3D cube, or only the 2D planes along the p-axis.
        
```