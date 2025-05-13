```
usage: pcdl_render_neuroglancer [-h] [--timestep TIMESTEP]
                                [--intensity_cmap INTENSITY_CMAP]
                                [tiffpathfile]

function to load a time step from an ome tiff files, generated with
make_ome_tiff, into neuroglancer.

positional arguments:
  tiffpathfile          path to ome tiff file.

options:
  -h, --help            show this help message and exit
  --timestep TIMESTEP   time step, within a possibly collapsed ome tiff file,
                        to render. the default will work with single time step
                        ome tiff files.
  --intensity_cmap INTENSITY_CMAP
                        matlab color map label, used to display expression
                        intensity values. if None, no intensity layers will be
                        generated. https://matplotlib.org/stable/users/explain
                        /colors/colormaps.html

homepage: https://github.com/elmbeech/physicelldataloader
```
