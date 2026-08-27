# mcdsts.make_simularium()


## input:
```
    focus_cat: list of 1 or 2 string; default is ['cell_type','current_phase']
        specify 1 or 2 categorical column labels, to be found
        in mcdsts.get_cell_df().

    trajectory_title: string; default 'timeseries'
        the trajectory_title will be used as
        <trajectory_title>.simularium file name and displayed
        in the simulation.

    scale_factor: float number; default is None
        A multiplier for the scene, use if visualization is
        too large or small. If None is provided, one will
        be calculated based on the position data.
        bue 20260822: does not seem to work in current simularium 1.13.0.

    camera_defaults: simulariumio.CameraData object; default is None
        camera's initial settings which it also returns to
        when reset.

    model_meta_data: simulariumio.ModelMetaData; default is None
        Metadata for the model that produced this
        trajectory.

```

## output:
```
    <trajectory_title>.simularium file

```

## description:
```
    function returns a simularium trajectory file that can be run with the
    online Simularium Viewer.
    + https://simularium.allencell.org/

```