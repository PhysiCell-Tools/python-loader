# mcds.get_conc_df()


## input:
```
            z_slice: floating point number; default is None
                z-axis position to slice a 2D xy-plain out of the
                3D substrate concentration mesh. if None the
                whole 3D mesh will be returned.

            halt: boolean; default is False
                should program execution break or just spit out a warning,
                if z_slice position is not an exact mesh center coordinate?
                if False, z_slice will be adjusted to the nearest
                mesh center value, the smaller one, if the coordinate
                lies on a saddle point.

            values: integer; default is 1
                minimal number of values a variable has to have to be outputted.
                variables that have only 1 state carry no information.
                None is a state too.

            drop: set of strings; default is an empty set
                set of column labels to be dropped for the dataframe.
                don't worry: essential columns like ID, coordinates
                and time will never be dropped.
                Attention: when the keep parameter is given, then
                the drop parameter has to be an empty set!

            keep: set of strings; default is an empty set
                set of column labels to be kept in the dataframe.
                set values=1 to be sure that all variables are kept.
                don't worry: essential columns like ID, coordinates
                and time will always be kept.

```

## output:
```
            df_conc : pandas dataframe
                dataframe stores all substrate concentrations in each voxel.

```

## description:
```
            function returns a dataframe with concentration values
            for all chemical species in all voxels. additionally, this
            dataframe lists voxel and mesh center coordinates.
        
```