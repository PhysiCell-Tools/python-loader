#########
# title: pyMCDSts.py
#
# language: python3
# date: 2022-08-22
# license: BSD-3-Clause
# authors: Patrick Wall, Randy Heiland, Paul Macklin, Elmar Bucher
#
# description:
#     pyMCDSts.py defines an object class, able to load and access
#     within python a time series of mcds objects loaded from a single
#     PhysiCell model output directory. pyMCDSts.py was first forked from
#     PhysiCell-Tools python-loader, where it was implemented as
#     pyMCDS_timeseries.py, then totally rewritten and further developed.
#
#     the make_image and make_movie functions are cloned from the PhysiCell
#     Makefile. note on difference image magick convert and mogrify:
#     + https://graphicsmagick-tools.narkive.com/9Sowc4HF/gm-tools-mogrify-vs-convert
#########


# load libraries
import matplotlib.pyplot as plt
from glob import glob
import numpy as np
import os
import pathlib
from pcdl import pdplt
from pcdl.pyMCDS import pyMCDS, es_coor_cell, es_coor_conc
import platform
import xml.etree.ElementTree as ET


# classes
class pyMCDSts:
    """
    input:
        output_path: string, default '.'
            relative or absolute path to the directory where
            the PhysiCell output files are stored.

        custom_type: dictionary; default is {}
            variable to specify custom_data variable types
            other than float (int, bool, str) like this: {var: dtype, ...}.
            down stream float and int will be handled as numeric,
            bool as Boolean, and str as categorical data.

        load: boole; default True
            should the whole time series data, all time steps, straight at
            object initialization be read and stored to mcdsts.l_mcds?

        microenv: boole; default True
            should the microenvironment be extracted?
            setting microenv to False will use less memory and speed up
            processing, similar to the original pyMCDS_cells.py script.

        graph: boole; default True
            should the graphs be extracted?
            setting graph to False will use less memory and speed up processing.

        settingxml: boole; default True
            should the substrate and cell type ID label mapping defined in
            PhysiCell_settings.xml be extracted?
            only set to False if the xml file is missing!

        verbose: boole; default True
            setting verbose to False for less text output, while processing.

    output:
        mcdsts: pyMCDSts class instance
            this instance offers functions to process all stored time steps
            from a simulation.

    description:
        pyMCDSts.__init__ generates a class instance the instance offers
        functions to process all time steps in the output_path directory.
    """
    def __init__(self, output_path='.', custom_type={}, load=True, microenv=True, graph=True, settingxml=True, verbose=True):
        output_path = output_path.replace('\\','/')
        if (output_path[-1] != '/'):
            output_path = output_path + '/'
        if not os.path.isdir(output_path):
            print(f'Error @ pyMCDSts.__init__ : this is not a path! could not load {output_path}.')
        self.output_path = output_path
        self.custom_type = custom_type
        self.microenv = microenv
        self.graph = graph
        self.settingxml = settingxml
        self.verbose = verbose
        if load:
            self.read_mcds()


    ## LOAD DATA
    def get_xmlfile_list(self):
        """
        input:
            self: pyMCDSts class instance.

        output:
            xmlfile_list: list of strings
                alphanumerical sorted list of /path/to/output*.xml strings.

        description:
            function returns an alphanumerical (and as such chronological)
            ordered list of physicell xml path and output file names. the
            list can be manipulated and used as input for the
            mcdsts.read_mcds function.
        """
        # bue 2022-10-22: is output*.xml always the correct pattern?
        ls_pathfile = [o_pathfile.as_posix() for o_pathfile in sorted(pathlib.Path(self.output_path).glob('output*.xml'))]
        return ls_pathfile


    def read_mcds(self, xmlfile_list=None):
        """
        input:
            self: pyMCDSts class instance.

            xmlfile_list: list of strings; default None
                list of physicell output /path/to/output*.xml strings.

        output:
            l_mcds: list of mcds objects

        description:
            the function returns a list of mcds objects loaded by
            pyMCDS calls.
        """
        # handle input
        if (xmlfile_list is None):
            xmlfile_list = self.get_xmlfile_list()

        # load mcds objects into list
        l_mcds = []
        for s_pathfile in xmlfile_list:
            mcds = pyMCDS(
                xmlfile = s_pathfile,
                custom_type = self.custom_type,
                microenv = self.microenv,
                graph = self.graph,
                settingxml = self.settingxml,
                verbose = self.verbose
            )
            l_mcds.append(mcds)
            if self.verbose:
                print() # carriage return

        # output
        self.l_mcds = l_mcds
        return l_mcds


    ## TRIAGE DATA
    def get_cell_df_columns_min_states(self, states=2):
        """
        input:
            states: integer; default is 2
                minimal number of states a variable has to have,
                in any of the mcds time steps, to be outputted.
                variables that have only 1 state carry no information.
                None is a state too.

        output:
            ls_variable: list of strings
                list of all non-coordinate column names, that at least in
                one of the time steps or in between time steps reach
                the given minimal state count.

        description:
            function to detect informative variables in a time series.
            this function detects even variables, which have in each
            time step less than the minimal state count, but different values
            from time step to time step.
        """
        # processing
        des_variable = {}
        for mcds in self.l_mcds:
            df_cell = mcds.get_cell_df()
            for s_column in df_cell.columns:
                if not (s_column in es_coor_cell):
                    es_state = set(df_cell.loc[:,s_column])
                    try:
                        des_variable[s_column] = des_variable[s_column].union(es_state)
                    except KeyError:
                        des_variable.update({s_column: es_state})
        es_variable = set()
        for s_column, es_state in des_variable.items():
            if len(es_state) >= states:
                es_variable.add(s_column)
        # output
        ls_variable = sorted(es_variable)
        return ls_variable


    def get_conc_df_columns_min_states(self, states=2):
        """
        input:
            states: integer; default is 2
                minimal number of states a variable has to have,
                in any of the mcds time steps, to be outputted.
                variables that have only 1 state carry no information.
                None is a state too.

        output:
            ls_variable: list of strings
                list of all non-coordinate column names, that at least in
                one of the time steps or in between time steps reach
                the given minimal state count.

        description:
            function to detect informative substrate concentration variables
            in a time series. this function detects even variables, which have
            in each time step less than the minimal state count, but
            different values from time step to time step.
        """
        # processing
        des_variable = {}
        for mcds in self.l_mcds:
            df_conc = mcds.get_concentration_df()
            for s_column in df_conc.columns:
                if not (s_column in es_coor_conc):
                    es_state = set(df_conc.loc[:,s_column])
                    try:
                        des_variable[s_column] = des_variable[s_column].union(es_state)
                    except KeyError:
                        des_variable.update({s_column: es_state})
        es_variable = set()
        for s_column, es_state in des_variable.items():
            if len(es_state) >= states:
                es_variable.add(s_column)
        # output
        ls_variable = sorted(es_variable)
        return ls_variable


    ## GENERATE AND TRANSFORM IMAGES
    def _handle_magick(self):
        """
        input:
            self: pyMCDSts class instance.

        output:
            s_magick: string
                image magick command line command call.

        description:
            internal function manipulates the command line command call,
            so that the call as well works on linux systems, which all
            too often run image magick < 7.0.
        """
        s_magick = 'magick '
        if (platform.system() in {'Linux'}) and (os.system('magick --version') != 0) and (os.system('convert --version') == 0):
            s_magick = ''
        return s_magick


    def make_imgcell(self, focus='cell_type', z_slice=0, z_axis=None, cmap='viridis', grid=True, xlim=None, ylim=None, xyequal=True, s=None, figsizepx=None, ext='jpeg', figbgcolor=None):
        """
        input:
            self: pyMCDSts class instance

            focus: string; default is 'cell_type'
                column name within cell dataframe.

            z_slice: floating point number; default is 0
                z-axis position to slice a 2D xy-plain out of the
                3D substrate concentration mesh. if z_slize position
                is not an exact mesh center coordinate, then z_slice
                will be adjusted to the nearest mesh center value,
                the smaller one, if the coordinate lies on a saddle point.

            z_axis: for a categorical focus: set of labels;
               for a numeric focus: tuple of two floats; default is None
               depending on the focus column variable dtype, default extracts
               labels or min and max values from data.

            cmap: string; default viridis.
                matplotlib colormap.
                https://matplotlib.org/stable/tutorials/colors/colormaps.html

            grid: boolean default True.
                plot axis grid lines.

            xlim: tuple of two floats; default is None
                x axis min and max value.
                default takes min and max from mesh x axis range.

            ylim: tuple of two floats; default is None
                y axis min and max value.
                default takes min and max from mesh y axis range.

            xyequal: boolean; default True
                to specify equal axis spacing for x and y axis.

            s: integer; default is None
                scatter plot dot size in pixel.
                typographic points are 1/72 inch.
                the marker size s is specified in points**2.
                plt.rcParams['lines.markersize']**2 is in my case 36.
                None tries to take the value from the initial.svg file.
                fall back setting is 36.

            figsizepx: list of two integers, default is None
                size of the figure in pixels, (x, y).
                the given x and y will be rounded to the nearest even number,
                to be able to generate movies from the images.
                None tries to take the values from the initial.svg file.
                fall back setting is [640, 480].

            ext: string
                output image format. possible formats are jpeg, png, and tiff.

            figbgcolor: string; default is None which is transparent (png)
                or white (jpeg, tiff).
                figure background color.

        output:
            image files under the returned path.

        description:
            this function generates image time series
            based on the physicell data output.
            jpeg is by definition a lossy compressed image format.
            png is by definition a lossless compressed image format.
            tiff can by definition be a lossy or lossless compressed format.
            https://en.wikipedia.org/wiki/JPEG
            https://en.wikipedia.org/wiki/Portable_Network_Graphics
            https://en.wikipedia.org/wiki/TIFF
        """
        # handle initial.svg for s and figsizepx
        if (s is None) or (figsizepx is None):
            s_pathfile = self.output_path + 'initial.svg'
            try:
                tree = ET.parse(s_pathfile)
                root = tree.getroot()
                if s is None:
                    circle_element = root.find('.//{*}circle')
                    r_radius = float(circle_element.get('r')) # px
                    s = int(round((r_radius)**2))
                    if self.verbose:
                        print(f's set to {s}.')
                if figsizepx is None:
                    i_width = int(np.ceil(float(root.get('width')))) # px
                    i_height = int(np.ceil(float(root.get('height'))))  # px
                    figsizepx = [i_width, i_height]
            except FileNotFoundError:
                print(f'Warning @ pyMCDSts.make_imgcell : could not load {s_pathfile}.')
                if s is None:
                    s = plt.rcParams['lines.markersize']**2
                    if self.verbose:
                        print(f's set to {s}.')
                if figsizepx is None:
                    figsizepx = [640, 480]

        # handle z_slice
        _, _, ar_p_axis = self.l_mcds[0].get_mesh_mnp_axis()
        if not (z_slice in ar_p_axis):
            z_slice = ar_p_axis[(ar_p_axis - z_slice).argmin()]
            print(f'z_slice set to {z_slice}.')

        # handle z_axis categorical cases
        df_cell = self.l_mcds[0].get_cell_df()
        if (str(df_cell.loc[:,focus].dtype) in {'bool', 'object'}):
            lr_extrema = [None, None]
            if (z_axis is None):
                # extract set of labels from data
                es_label = set()
                for mcds in self.l_mcds:
                    df_cell = mcds.get_cell_df()
                    es_label = es_label.union(set(df_cell.loc[:,focus]))
            else:
                es_label = z_axis

        # handle z_axis numerical cases
        else:  # df_cell.loc[:,focus].dtype is numeric
            es_label = None
            if (z_axis is None):
                # extract min and max values from data
                lr_extrema = [None, None]
                for mcds in self.l_mcds:
                    df_cell = mcds.get_cell_df()
                    r_min = df_cell.loc[:,focus].min()
                    r_max = df_cell.loc[:,focus].max()
                    if (lr_extrema[0] is None) or (lr_extrema[0] > r_min):
                        lr_extrema[0] = np.floor(r_min)
                    if (lr_extrema[1] is None) or (lr_extrema[1] < r_max):
                        lr_extrema[1] = np.ceil(r_max)
            else:
                lr_extrema = z_axis

        # handle z_axis summary
        print(f'labels found: {es_label}.')
        print(f'min max extrema set to: {lr_extrema}.')

        # handle xlim and ylim
        if xlim is None:
            xlim = self.l_mcds[0].get_mesh_mnp_range()[0]
        if ylim is None:
            ylim = self.l_mcds[0].get_mesh_mnp_range()[1]

        # handle figure size
        figsizepx[0] = figsizepx[0] - (figsizepx[0] % 2)  # enforce even pixel number
        figsizepx[1] = figsizepx[1] - (figsizepx[1] % 2)
        r_px = 1 / plt.rcParams['figure.dpi']  # translate px to inch
        figsize = [None, None]
        figsize[0] = figsizepx[0] * r_px
        figsize[1] = figsizepx[1] * r_px
        if self.verbose:
            print(f'px figure size set to {figsizepx}.')

        # handle figure background color
        if figbgcolor is None:
            figbgcolor = 'auto'

        # handle output path
        s_path = f'{self.output_path}cell_{focus}_z{z_slice}/'

        # plotting
        for mcds in self.l_mcds:
            fig, ax = plt.subplots(figsize=figsize)
            df_cell = mcds.get_cell_df()
            df_cell = df_cell.loc[(df_cell.position_z == z_slice),:]
            if not (es_label is None):
                ds_color = pdplt.df_label_to_color(
                    df_abc = df_cell,
                    s_label = focus,
                    es_label = es_label,
                    s_cmap = cmap,
                    b_shuffle = False,
                )
                pdplt.ax_colorlegend(
                    ax = ax,
                    ds_color = ds_color,
                    r_x_figure2legend_space = 0.01,
                    s_fontsize = 'small',
                )
                s_focus_color = focus + '_color'
                c = list(df_cell.loc[:, s_focus_color].values)
                s_cmap = None
            else:
                c = focus
                s_cmap = cmap
            df_cell.plot(
                kind = 'scatter',
                x = 'position_x',
                y = 'position_y',
                c = c,
                vmin = lr_extrema[0],
                vmax = lr_extrema[1],
                cmap = s_cmap,
                xlim = xlim,
                ylim = ylim,
                s = s,
                grid = grid,
                title = f'{mcds.get_time()}[min] {df_cell.shape[0]}[agent]',
                ax = ax,
            )
            if xyequal:
                ax.axis('equal')
            os.makedirs(s_path, exist_ok=True)
            s_pathfile = f'{s_path}{focus}_{str(mcds.get_time()).zfill(11)}.{ext}'
            fig.savefig(s_pathfile, facecolor=figbgcolor)
            plt.close(fig)

        # output
        return s_path


    def make_imgsubs(self, focus, z_slice=0, extrema=None, alpha=1, fill=True, cmap='viridis', grid=True, xlim=None, ylim=None, xyequal=True, figsizepx=None, ext='jpeg', figbgcolor=None):
        """
        input:
            self: pyMCDSts class instance

            focus: string; default is 'cell_type'
                column name within cell dataframe.

            z_slice: floating point number; default is 0
                z-axis position to slice a 2D xy-plain out of the
                3D substrate concentration mesh. if z_slize position
                is not an exact mesh center coordinate, then z_slice
                will be adjusted to the nearest mesh center value,
                the smaller one, if the coordinate lies on a saddle point.

            extrema: tuple of two floats; default is None
                default takes min and max from data.

            alpha: floating point number; default is 1
                alpha channel transparency value
                between 1 (not transparent at all) and 0 (totally transparent).

            fill: boolean
                True generates a matplotlib contourf plot.
                False generates a matplotlib contour plot.

            cmap: string; default viridis.
                matplotlib colormap.
                https://matplotlib.org/stable/tutorials/colors/colormaps.html

            grid: boolean default True.
                plot axis grid lines.

            xlim: tuple of two floats; default is None
                x axis min and max value.
                default takes min and max from mesh x axis range.

            ylim: tuple of two floats; default is None
                y axis min and max value.
                default takes min and max from mesh y axis range.

            xyequal: boolean; default True
                to specify equal axis spacing for x and y axis.

            figsizepx: list of two integers, default is None
                size of the figure in pixels, (x, y).
                the given x and y will be rounded to the nearest even number,
                to be able to generate movies from the images.
                None tries to take the values from the initial.svg file.
                fall back setting is [640, 480].

            ext: string
                output image format. possible formats are jpeg, png, and tiff.

            figbgcolor: string; default is None which is transparent (png)
                or white (jpeg, tiff).
                figure background color.

        output:
            image files under the returned path.

        description:
            this function generates image time series
            based on the physicell data output.
            jpeg is by definition a lossy compressed image format.
            png is by definition a lossless compressed image format.
            tiff can by definition be a lossy or lossless compressed format.
            https://en.wikipedia.org/wiki/JPEG
            https://en.wikipedia.org/wiki/Portable_Network_Graphics
            https://en.wikipedia.org/wiki/TIFF
        """
        # handle initial.svg for s and figsizepx
        if (figsizepx is None):
            s_pathfile = self.output_path + 'initial.svg'
            try:
                tree = ET.parse(s_pathfile)
                root = tree.getroot()
                i_width = int(np.ceil(float(root.get('width')))) # px
                i_height = int(np.ceil(float(root.get('height'))))  # px
                figsizepx = [i_width, i_height]
            except FileNotFoundError:
                print(f'Warning @ pyMCDSts.make_imgsubs : could not load {s_pathfile}.')
                figsizepx = [640, 480]

        # handle z_slice
        _, _, ar_p_axis = self.l_mcds[0].get_mesh_mnp_axis()
        if not (z_slice in ar_p_axis):
            z_slice = ar_p_axis[(ar_p_axis - z_slice).argmin()]
            print(f'z_slice set to {z_slice}.')

        # handle extrema
        if extrema == None:
            extrema = [None, None]
            for mcds in self.l_mcds:
                df_cell = mcds.get_cell_df()
                r_min = df_cell.loc[:,focus].min()
                r_max = df_cell.loc[:,focus].max()
                if (extrema[0] is None) or (extrema[0] > r_min):
                    extrema[0] = np.floor(r_min)
                if (extrema[1] is None) or (extrema[1] < r_max):
                    extrema[1] = np.ceil(r_max)
            print(f'min max extrema set to {extrema}.')

        # handle xlim and ylim
        if (xlim is None):
            xlim = self.l_mcds[0].get_mesh_mnp_range()[0]
        if (ylim is None):
            ylim = self.l_mcds[0].get_mesh_mnp_range()[1]

        # handle figure size
        figsizepx[0] = figsizepx[0] - (figsizepx[0] % 2)  # enforce even pixel number
        figsizepx[1] = figsizepx[1] - (figsizepx[1] % 2)
        r_px = 1 / plt.rcParams['figure.dpi']  # translate px to inch
        figsize = [None, None]
        figsize[0] = figsizepx[0] * r_px
        figsize[1] = figsizepx[1] * r_px
        if self.verbose:
            print(f'px figure size set to {figsizepx}.')

        # handle figure background color
        if figbgcolor is None:
            figbgcolor = 'auto'

        # handle output path
        s_path = f'{self.output_path}substrate_{focus}_z{z_slice}/'

        # plotting
        for mcds in self.l_mcds:
            fig = mcds.get_contour(
                substrate = focus,
                z_slice = z_slice,
                vmin = extrema[0],
                vmax = extrema[1],
                alpha = alpha,
                fill = fill,
                cmap = cmap,
                title = f'{mcds.get_time()}[min] {df_cell.shape[0]}[agent]',
                grid = grid,
                xlim = xlim,
                ylim = ylim,
                xyequal = xyequal,
                figsize = figsize,
                ax = None,
            )
            os.makedirs(s_path, exist_ok=True)
            s_pathfile = f'{s_path}{focus}_{str(mcds.get_time()).zfill(11)}.{ext}'
            fig.savefig(s_pathfile, facecolor=figbgcolor)
            plt.close(fig)

        # output
        return s_path


    def make_gif(self, path, interface='jpeg'):
        """
        input:
            self: pyMCDSts class instance.

            path: string
                relative or absolute path to where the images are
                from which the gif will be generated.

            interface: string; default jpeg
                this images, from which the gif will be generated
                have to exist under the given path.
                they can be generated with the make_imgcell or make_imgsubs
                function.

        output:
            gif file in the path directory.
`               additionally, the function will return the gif's path and filename.

        description:
            this function generates a gif image from all interface image files
            found in the path directory.
            https://en.wikipedia.org/wiki/GIF
        """
        s_magick = self._handle_magick()
        # handle path and file name
        path = path.replace('\\','/')
        if path.endswith('/'): path = path[:-1]
        s_file = path.split('/')[-1]
        if s_file.startswith('.'): s_file = s_file[1:]
        if (len(s_file) == 0): s_file = 'movie'
        s_file += f'_{interface}.gif'
        s_opathfile = f'{path}/{s_file}'
        s_ipathfiles = f'{path}/*.{interface}'
        # genaerate gif
        os.system(f'{s_magick}convert {s_ipathfiles} {s_opathfile}')

        # output
        return s_opathfile


    def make_movie(self, path, interface='jpeg', framerate=12):
        """
        input:
            self: pyMCDSts class instance.

            path: string
                relative or absolute path to where the images are
                from which the movie will be generated.

            interface: string; default jpeg
                this images, from which the mp4 movie will be generated
                have to exist under the given path.
                they can be generated with the make_imgcell or make_imgsubs
                function.

            framerate: integer; default 24
                specifies how many images per second will be used.

        output:
            mp4 move file in the path directory.
`               additionally, the function will return the mp4's path and filename.

        description:
            this function generates a movie from all interface image files
            found in the path directory.
            https://en.wikipedia.org/wiki/MP4_file_format
            https://en.wikipedia.org/wiki/Making_Movies
        """
        # handle path and file name
        path = path.replace('\\','/')
        if path.endswith('/'): path = path[:-1]
        s_file = path.split('/')[-1]
        if s_file.startswith('.'): s_file = s_file[1:]
        if (len(s_file) == 0): s_file = 'movie'
        s_file += f'_{interface}{framerate}.mp4'
        s_opathfile = f'{path}/{s_file}'
        s_ipathfiles = f'{path}/*.{interface}'

        # generate movie
        s_cmd = f'ffmpeg -r {framerate} -f image2 -pattern_type glob -i "{s_ipathfiles}" -vcodec libx264 -pix_fmt yuv420p -strict -2 -tune animation -crf 15 -acodec none {s_opathfile}'
        os.system(s_cmd)

        # output
        return s_opathfile
