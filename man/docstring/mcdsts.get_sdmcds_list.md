# mcdsts.get_sdmcds_list()


## input:
```
    self: TimeSeries class instance.

```

## output:
```
    self.l_sdmcds: list of chronologically ordered spatialdata mcds objects.
        watch out, this is a pointer to the
        self.l_sdmcds list of spdata mcds objects, not a copy of self.l_sdmcds!

```

## description:
```
    function returns a binding to the self.l_sdmcds list of spdata mcds objects.

```