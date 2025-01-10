# mcdsts.get_annmcds_list()


## input:
```
            self: TimeSeries class instance.

```

## output:
```
            self.l_annmcds: list of chronologically ordered anndata mcds objects.
                watch out, this is a pointer to the
                self.l_annmcds list of anndata mcds objects, not a copy of self.l_annmcds!

```

## description:
```
            function returns a binding to the self.l_annmcds list of anndata mcds objects.
        
```