```
usage: pcdl_plot_scatter [-h] [--custom_data_type [CUSTOM_DATA_TYPE ...]]
                         [--microenv MICROENV] [--physiboss PHYSIBOSS]
                         [--settingxml SETTINGXML] [-v VERBOSE]
                         [--z_slice Z_SLICE] [--z_axis Z_AXIS [Z_AXIS ...]]
                         [--alpha ALPHA] [--cmap CMAP] [--title TITLE]
                         [--grid GRID] [--legend_loc LEGEND_LOC]
                         [--xlim XLIM [XLIM ...]] [--ylim YLIM [YLIM ...]]
                         [--xyequal XYEQUAL] [--s S]
                         [--figsizepx FIGSIZEPX [FIGSIZEPX ...]] [--ext EXT]
                         [--figbgcolor FIGBGCOLOR]
                         [path] [focus]

function generates pandas scatter plots, under the returned path.

positional arguments:
  path                  path to the PhysiCell output directory or a
                        outputnnnnnnnn.xml file. default is . .
  focus                 column name within conc dataframe. default is
                        cell_type.

options:
  -h, --help            show this help message and exit
  --custom_data_type [CUSTOM_DATA_TYPE ...]
                        parameter to specify custom_data variable types other
                        than float (namely: int, bool, str) like this
                        var:dtype myint:int mybool:bool mystr:str . downstream
                        float and int will be handled as numeric, bool as
                        Boolean, and str as categorical data. default is an
                        empty string.
  --microenv MICROENV   should the microenvironment data be loaded? setting
                        microenv to False will use less memory and speed up
                        processing, similar to the original pyMCDS_cells.py
                        script. default is True.
  --physiboss PHYSIBOSS
                        if found, should physiboss state data be extracted and
                        loaded into the df_cell dataframe? default is True.
  --settingxml SETTINGXML
                        the settings.xml that is loaded, from which the cell
                        type ID label mapping, is extracted, if this
                        information is not found in the output xml file. set
                        to None or False if the xml file is missing! default
                        is PhysiCell_settings.xml.
  -v VERBOSE, --verbose VERBOSE
                        setting verbose to False for less text output, while
                        processing. default is True.
  --z_slice Z_SLICE     z-axis position to slice a 2D xy-plain out of the 3D
                        mesh. if z_slice position numeric but not an exact
                        mesh center coordinate, then z_slice will be adjusted
                        to the nearest mesh center value, the smaller one, if
                        the coordinate lies on a saddle point. default is 0.0.
  --z_axis Z_AXIS [Z_AXIS ...]
                        for a categorical focus: list of labels; for a numeric
                        focus: list of two floats; None, depending on the
                        focus column variable dtype, extracts labels or min
                        and max values from data. default is None
  --alpha ALPHA         alpha channel transparency value between 1 (not
                        transparent at all) and 0 (totally transparent).
                        default is 1.0.
  --cmap CMAP           matplotlib colormap string from https://matplotlib.org
                        /stable/tutorials/colors/colormaps.html . default is
                        viridis.
  --title TITLE         title prefix. default is an empty string.
  --grid GRID           plot axis grid lines. default is True.
  --legend_loc LEGEND_LOC
                        the location of the categorical legend, if applicable.
                        possible strings are: best, 'upper right', 'upper
                        center', 'upper left', 'center left', 'lower left',
                        'lower center', 'lower right', 'center right', center.
                        default is 'lower left'
  --xlim XLIM [XLIM ...]
                        two floats. x axis min and max value. None takes min
                        and max from mesh x axis range. default is None.
  --ylim YLIM [YLIM ...]
                        two floats. y axis min and max value. None takes min
                        and max from mesh y axis range. default is None.
  --xyequal XYEQUAL     to specify equal axis spacing for x and y axis.
                        default is True.
  --s S                 scatter plot dot size in pixel. typographic points are
                        1/72 inch. the marker size s is specified in
                        points**2. plt.rcParams['lines.markersize']**2 is in
                        my case 36. None tries to take the value from the
                        initial.svg file. fall back setting is 36. default is
                        None.
  --figsizepx FIGSIZEPX [FIGSIZEPX ...]
                        size of the figure in pixels (integer), x y. the given
                        x and y will be rounded to the nearest even number, to
                        be able to generate movies from the images. None tries
                        to take the values from the initial.svg file. fall
                        back setting is 640 480. default is None.
  --ext EXT             output image format. possible formats are jpeg, png,
                        and tiff. default is jpeg.
  --figbgcolor FIGBGCOLOR
                        figure background color. None is transparent (png) or
                        white (jpeg, tiff). default is None.

homepage: https://github.com/elmbeech/physicelldataloader
```
