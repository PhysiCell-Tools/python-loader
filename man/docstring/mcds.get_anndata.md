# mcds.get_anndata()


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
                set values=1 to be sure that all variables are kept.
                don't worry: essential columns like ID, coordinates
                and time will always be kept.

            scale: string; default 'maxabs'
                specify how the data should be scaled.
                possible values are None, maxabs, minmax, std.
                for more input, check out: help(pcdl.scaler)

```

## output:
```
            annmcds: anndata object
                for this one time step.

```

## description:
```
            function to transform a mcds time step into an anndata object
            for downstream analysis.
        
```