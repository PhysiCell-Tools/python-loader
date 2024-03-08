# mcdsts.make_movie('path/to/plots')


## input:
```
        path: string
            relative or absolute path to where the images are
            from which the movie will be generated.

        interface: string; default jpeg
            this images, from which the mp4 movie will be generated
            have to exist under the given path.
            they can be generated with the plot_scatter or plot_contour
            function.

        framerate: integer; default 12
            specifies how many images per second will be used.
            humans are capable of processing 12 images per second and
            seeing them individually. higher rates are seen as motion.

```

## output:
```
        mp4 move file in the path directory.
`           additionally, the function will return the mp4's path and filename.

```

## description:
```
        this function generates a movie from all interface image files
        found in the path directory.
        https://en.wikipedia.org/wiki/MP4_file_format
        https://en.wikipedia.org/wiki/Making_Movies
    
```