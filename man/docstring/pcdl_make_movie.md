```
usage: pcdl_make_movie [-h] [--framerate FRAMERATE] [path] [interface]

this function generates an mp4 movie file from all image files found in the
path directory in the specified interface file format.

positional arguments:
  path                  relative or absolute path to where the images are from
                        which the mp4 movie will be generated. default is . .
  interface             specify the image format from which the mp4 movie will
                        be generated. these images have to exist under the
                        given path. they can be generated with the
                        plot_scatter or plot_contour function. default is
                        jpeg.

options:
  -h, --help            show this help message and exit
  --framerate FRAMERATE
                        specifies how many images per second will be used.
                        humans are capable of processing 12 images per second
                        and seeing them individually. higher rates are seen as
                        motion. default is 12.

homepage: https://github.com/elmbeech/physicelldataloader
```
