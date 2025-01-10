# mcdsts.make_gif('path/to/plots')


## input:
```
        path: string
            relative or absolute path to where the images are
            from which the gif will be generated.

        interface: string; default jpeg
            this images, from which the gif will be generated
            have to exist under the given path.
            they can be generated with the plot_scatter or plot_contour
            function.

```

## output:
```
        gif file in the path directory.
            additionally, the function will return the gif's path and filename.

```

## description:
```
        this function generates a gif image from all interface image files
        found in the path directory.
        https://en.wikipedia.org/wiki/GIF
    
```