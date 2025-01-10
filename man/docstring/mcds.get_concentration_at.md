# mcds.get_concentration_at()


## input:
```
            x: floating point number
                position x-coordinate.

            y: floating point number
                position y-coordinate.

            z: floating point number; default is 0
                position z-coordinate.

```

## output:
```
            ar_concs: numpy array of floating point numbers
                array of substrate concentrations in the order
                given by get_substrate_list().

```

## description:
```
            function return concentrations of each chemical species
            inside a particular voxel that contains the point specified
            in the arguments.
        
```