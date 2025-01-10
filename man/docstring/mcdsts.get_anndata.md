# mcdsts.get_anndata()


## input:
```
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
                don't worry: essential columns like ID, coordinates
                and time will always be kept.

            scale: string; default 'maxabs'
                specify how the data should be scaled.
                possible values are None, maxabs, minmax, std.
                for more input, check out: help(pcdl.scaler)

            collapse: boole; default True
                should all mcds time steps from the time series be collapsed
                into one single anndata object, or a list of anndata objects
                for each time step?

            keep_mcds: boole; default True
                should the loaded original mcds be kept in memory
                after transformation?

```

## output:
```
            annmcds or self.l_annmcds: anndata object or list of anndata objects.
                what is returned depends on the collapse setting.

```

## description:
```
            function to transform mcds time steps into one or many
            anndata objects for downstream analysis.
        
```