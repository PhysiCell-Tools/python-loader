# mcds.plot_contour()


## input:
```
            substrate: string
                substrate name.

            z_slice: floating point number; default is 0
                z-axis position to slice a 2D xy-plain out of the
                3D substrate concentration mesh. if z_slice position
                is not an exact mesh center coordinate, then z_slice
                will be adjusted to the nearest mesh center value,
                the smaller one, if the coordinate lies on a saddle point.

            vmin: floating point number; default is None
                color scale min value.
                None will take the min value found in the data.

            vmax: floating point number; default is None
                color scale max value.
                None will take the max value found in the data.

            alpha: floating point number; default is 1
                alpha channel transparency value
                between 1 (not transparent at all) and 0 (totally transparent).

            fill: boolean; default is True
                True generates a matplotlib contourf plot.
                False generates a matplotlib contour plot.

            cmap: string; default is viridis
                matplotlib color map color label.
                https://matplotlib.org/stable/tutorials/colors/colormaps.html

            title: string; default None
                possible plot title string.

            grid: boolean; default True
                should be plotted on  a grid or on a blank page?
                True will plot on a grid.

            xlim: tuple of two floating point numbers; default is None
                to specify min and max x axis value.
                None will extract agreeable values from the data.

            ylim: tuple of two floating point numbers; default is None
                to specify min and max y axis value.
                None will extract agreeable values from the data.

            xyequal: boolean; default True
                to specify equal axis spacing for x and y axis.

            figsize: tuple of floating point numbers; default is None
                the specif the figure x and y measurement in inch.
                None result in the default matplotlib setting, which is [6.4, 4.8].

            ax: matplotlib axis object; default setting is None
                the ax object, which will be used as a canvas for plotting.
                None will generate a figure and ax object from scratch.

```

## output:
```
            fig: matplotlib figure, containing the ax axis object,
                with contour plot and color bar.

```

## description:
```
            function returns a matplotlib contour (or contourf) plot,
            inclusive color bar, for the substrate specified.
        
```