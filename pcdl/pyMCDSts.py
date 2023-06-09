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
import platform
from .pyMCDS import pyMCDS
#import shutil

# classes
class pyMCDSts:
    """
    input:
        output_path: string, default '.'
            relative or absolute path to the directory where
            the PhysiCell output files are stored.

        microenv: booler; default True
            should the microenvironment be extracted?
            setting microenv to False will use less memory and speed up
            processing, similar to the original pyMCDS_cells.py script.

        graph: boole; default True
            should the graphs be extracted?
            setting graph to False will use less memory and speed up processing.

        verbose: boole; default True
            setting verbose to False for less text output, while processing.

    output:
        mcdsts: pyMCDSts class instance
            this instance offers functions to process all stored time steps
            from a simulation. no data is fetched by initialization.

    description:
        pyMCDSts.__init__ generates a class instance and stores
        the input parameters. no data is fetched at initialization.
        the instance offers functions to process all time steps
        in the output_path directory.
    """
    def __init__(self, output_path='.', microenv=True, graph=True, verbose=True):
        self.output_path = output_path
        self.microenv = microenv
        self.graph = graph
        self.verbose = False
        self.l_mcds = self.read_mcds()
        self.verbose = verbose


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
        return(ls_pathfile)


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
                microenv = self.microenv,
                graph = self.graph,
                verbose = self.verbose
            )
            l_mcds.append(mcds)
            if self.verbose:
                print() # carriage return

        # output
        self.l_mcds = l_mcds
        return(l_mcds)


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
        return(s_magick)


    def make_imgcell(self, focus='cell_type', z_slice=0, extrema=None, cmap='viridis', grid=True, xlim=None, ylim=None, s=plt.rcParams['lines.markersize']**2, figsizepx=[640, 480], ext='jpeg', figbgcolor=None):
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

            s: integer; default is plt.rcParams['lines.markersize']**2
                scatter plot dot size in pixel.
                Typographic points are 1/72 inch.
                The marker size s is specified in points**2.
                plt.rcParams['lines.markersize']**2 is in my case 36.

            figsizepx: list of two integers, default is [640, 480]
                size of the figure in pixels, (x, y).
                the given x and y will be rounded to the nearest even number,
                to be able to generate movies from the images.

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
        s_path = f'{self.output_path}cell_{focus}_z{z_slice}/'
        #if os.path.exists(s_path):
        #    shutil.rmtree(s_path)

        # plotting
        for mcds in self.l_mcds:
            df_cell = mcds.get_cell_df()
            df_cell = df_cell.loc[(df_cell.position_z == z_slice), :]
            fig, ax = plt.subplots(figsize=figsize)
            df_cell.plot(
                kind = 'scatter',
                x = 'position_x',
                y = 'position_y',
                c = focus,
                vmin = extrema[0],
                vmax = extrema[1],
                cmap = cmap,
                xlim = xlim,
                ylim = ylim,
                s = s,
                grid = grid,
                title = f'{mcds.get_time()}[min] {df_cell.shape[0]}[agent]',
                ax = ax,
            )
            os.makedirs(s_path, exist_ok=True)
            s_pathfile = f'{s_path}{focus}_{str(mcds.get_time()).zfill(11)}.{ext}'
            fig.savefig(s_pathfile, facecolor=figbgcolor)
            plt.close(fig)

        # output
        return(s_path)


    def make_imgsubs(self, focus, z_slice=0, extrema=None, alpha=1, fill=True, cmap='viridis', grid=True, xlim=None, ylim=None, figsizepx=[640, 480], ext='jpeg', figbgcolor=None):
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

            figsizepx: tuple of two integers, default is (640, 480)
                size of the figure in pixels, (x, y).
                the given x and y will be rounded to the nearest even number,
                to be able to generate movies from the images.

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
        #if os.path.exists(s_path):
        #    shutil.rmtree(s_path)

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
                figsize = figsize,
                ax = None,
            )
            os.makedirs(s_path, exist_ok=True)
            s_pathfile = f'{s_path}{focus}_{str(mcds.get_time()).zfill(11)}.{ext}'
            fig.savefig(s_pathfile, facecolor=figbgcolor)
            plt.close(fig)

        # output
        return(s_path)


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
        return(s_opathfile)


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
        return(s_opathfile)
