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
import os
import pathlib
import platform
from .pyMCDS import pyMCDS
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from glob import glob

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


    def read_mcds(self, xmlfile_list=None) -> list[pyMCDS]:
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
        return(l_mcds)


    ## TRANSFORM SVG
    def _handle_magick(self):
        """
        input:
            self: pyMCDSts class instance.

        output:
            s_magick: string
                image magick command line command call

        description:
            internal function manipulates the command line command call,
            so that the call as well works on linux systems, which all
            too often run image magick < 7.0
        """
        s_magick = 'magick '
        if (platform.system() in {'Linux'}) and (os.system('magick --version') != 0) and (os.system('convert --version') == 0):
            s_magick = ''
        return(s_magick)
    
    # def make_jpeg(self, resize_factor=1):
    #     """
    #     input:
    #         self: pyMCDSts class instance.
    #         glob: string
    #             wildcard filename pattern.
    #         resize_factor: floating point number; default 1
    #             to specify image magnification or scale down.
    #             the resize parameter will in any case be adjusted,
    #             so that the resulting image's height and width are
    #             integer divisible by 2. this is because of a
    #             ffmpeg constrain for generating a movie out of images.
    #     output:
    #         jpeg files in output_path directory.
    #     description:
    #         this function generates jpeg image equivalents from all svg files
    #         found in the output_path directory.
    #         jpeg is by definition a lossy compressed image format.
    #         https://en.wikipedia.org/wiki/JPEG
    #     """
    #     # bue: use mogrify, convert might cause troubles here!
    #     s_magick = self._handle_magick()
    #     s_resize = self._handle_resize(resize_factor=resize_factor)
    #     for s_glob in ls_glob:
    #         if (len(set(pathlib.Path(self.output_path).glob(s_glob))) > 0):
    #             if (s_glob in es_resize):
    #                 os.system(f'{s_magick}mogrify {s_resize} -format jpeg {self.output_path}/{s_glob} &')
    #             else:
    #                 os.system(f'{s_magick}mogrify -format jpeg {self.output_path}/{s_glob} &')


    # def make_png(self, resize_factor=1, addargs='-transparent white'):
    #     """
    #     input:
    #         self: pyMCDSts class instance.
    #         resize_factor: floating point number; default 1
    #             to specify image magnification or scale down.
    #             the resize parameter will in any case be adjusted,
    #             so that the resulting image's height and width are
    #             integer divisible by 2. this is because of a
    #             ffmpeg constrain for generating a movie out of images.
    #         addargs: string; default '-transparent white'
    #             sting to additional image magick parameters.
    #             by default, alpha channel transparency is set to white.
    #     output:
    #         png files in output_path directory.
    #     description:
    #         this function generates png image equivalents from all svg files
    #         found in the output_path directory.
    #         png is by definition a lossless compressed image format.
    #         https://en.wikipedia.org/wiki/Portable_Network_Graphics
    #     """
    #     # bue: use mogrify, convert might cause troubles here!
    #     s_magick = self._handle_magick()
    #     s_resize = self._handle_resize(resize_factor=resize_factor)
    #     for s_glob in ls_glob:
    #         if (len(set(pathlib.Path(self.output_path).glob(s_glob))) > 0):
    #             if (s_glob in es_resize):
    #                 os.system(f'{s_magick}mogrify {s_resize} {addargs} -format png {self.output_path}/{s_glob} &')
    #             else:
    #                 os.system(f'{s_magick}mogrify {addargs} -format png {self.output_path}/{s_glob} &')


    # def make_tiff(self, resize_factor=1):
    #     """
    #     input:
    #         self: pyMCDSts class instance.
    #         resize_factor: floating point number; default 1
    #             to specify image magnification or scale down.
    #             the resize parameter will in any case be adjusted,
    #             so that the resulting image's height and width are
    #             integer divisible by 2. this is because of a
    #             ffmpeg constrain for generating a movie out of images.
    #     output:
    #         tiff files in output_path directory.
    #     decription:
    #         this function generates tiff image equivalents from all svg files
    #         found in the output_path directory.
    #         https://en.wikipedia.org/wiki/TIFF
    #     """
    #     # bue: use mogrify, convert might cause troubles here!
    #     s_magick = self._handle_magick()
    #     s_resize = self._handle_resize(resize_factor=resize_factor)
    #     for s_glob in ls_glob:
    #         if (len(set(pathlib.Path(self.output_path).glob(s_glob))) > 0):
    #             if (s_glob in es_resize):
    #                 os.system(f'{s_magick}mogrify {s_resize} -format tiff {self.output_path}/{s_glob} & ')
    #             else:
    #                 os.system(f'{s_magick}mogrify -format tiff {self.output_path}/{s_glob} & ')
    
    def make_imgcell(self, cmap, s, figsize, ext, focus='cellt_type', range=None, figbgcolor=None):
        """
        input:
            self: pyMCDSts class instance

            cmap: string
                matplotlib colormap

            s: integer
                scatter plot dot size in pixel
            
            figsize: tuple (int, int)
                size of the figure in pixels, (x, y)
                (given x and y will be rounded to nearest even number)

            ext: string
                output image interface (jpeg, png, tiff)
            
            focus: string; default is 'cellt_type'
                column name within cell dataframe;
                returns error if string not a valid name in df
            

            range: tuple (float, float); default is None
                default takes min and max from data
            
            figbgcolor: string; default is None (transparent)
                figure background color

        output:

        description:
        """
        # generate image of cell of given extension type
        ls_xmlfilelist = self.get_xmlfile_list()
        l_mcds = self.read_mcds(ls_xmlfilelist)

        for mcds in l_mcds:
            df_cell = mcds.get_cell_df()
            if focus not in df_cell.columns:
                raise Exception(f'make_imgcell : given focus string "{focus}" not a valid key in dataframe')
            
            # filter data
            if range is not None:
                df_filter = range[0] <= df_cell <= range[1]
                df_cell = df_cell[df_filter]
            
            # plot
            fig, ax = plt.subplots()
            df_cell.loc[:, focus].plot(
                kind='scatter', 
                x='position_x',
                y='position_y',
                c=focus,
                cmap=cmap,
                s=s,
                ax=ax
            )

            # set figure background color
            fig.set_facecolor(figbgcolor)

            # set figure size
            i_x = figsize[0] - (figsize[0] % 2) # enforce even number
            i_y = figsize[1] - (figsize[1] % 2)
            fig.set_figwidth(i_x)
            fig.set_figheight(i_y)

            # generate file name and save image file
            s_filename = f'{self.output_path}/{focus}_{mcds.get_time()}.{ext}'
            fig.savefig(fname=s_filename)

    def make_imgsubstrate(self, focus, cmap, s, figsize, ext, range=None, figbgcolor=None):
        """
        input:
            self: pyMCDSts class instance
            
            focus: string
                column name within cell dataframe;
                returns error if string not a valid name in df

            cmap: string
                matplotlib colormap

            s: integer
                scatter plot dot size in pixel
            
            figsize: tuple (int, int)
                size of the figure in pixels, (x, y)
                (given x and y will be rounded to nearest even number)

            ext: string
                output image interface (jpeg, png, tiff)
            
            range: tuple (float, float); default is None
                default takes min and max from data
            
            figbgcolor: string; default is None (transparent)
                figure background color
                

        output:

        description:
        """
        # generate image of substrate of given extension type
        ls_xmlfilelist = self.get_xmlfile_list()
        l_mcds = self.read_mcds(ls_xmlfilelist)

        for mcds in l_mcds:
            if focus not in mcds.get_substrate_names():
                raise Exception(f'make_imgsubstrate : given focus string "{focus}" not a valid substrate name')
            
            fig, ax = plt.subplots()
            fig = mcds.get_contour(substrate=focus, cmap=cmap, ax=ax)

            # set figure background color and size
            fig.set_facecolor(figbgcolor)
            
             # set figure size
            i_x = figsize[0] - (figsize[0] % 2) # enforce even number
            i_y = figsize[1] - (figsize[1] % 2)
            fig.set_figwidth(i_x)
            fig.set_figheight(i_y)

            # generate file name and save image file
            s_filename = f'{self.output_path}/{focus}_{mcds.get_time()}.{ext}'
            fig.savefig(fname=s_filename)

        

    def make_gif(self, focus, interface='jpeg'):
        """
        input:
            self: pyMCDSts class instance.

            focus: string
                specifies either cell variable or substrate

            interface: string; default jpeg
                ffmpeg cannot directly translate svg image into a move.
                the interface image format will be used to bridge the gap.
                this images, from which the movie will be generated, have to exist.
                they can be generated with the make_jpeg, make_png, or make_tiff
                function.

        output:
            gif file in output_path directory.
`               additionally, the function will return the path and filename.

        description:
            this function generates a gif image from all snapshot svg files
            found in the output_path directory.
        """
        s_magick = self._handle_magick()
        # bue: use convert, mogrify will cause troubles here!
        s_opathfile = f'{self.output_path}/{focus}_{interface}.gif'
        s_ipathfiles = f'{self.output_path}/{focus}_*.{interface}'

        # check files exist
        if not glob(s_ipathfiles): # list returned by glob (files matching wildcard) is empty
            raise Exception('make_gif : images not found in output directory for given focus string')
        
        # make system call
        os.system(f'{s_magick}convert {s_ipathfiles} {s_opathfile}')

        # output
        return(s_opathfile)

    def make_movie(self, focus, interface='jpeg', frame_rate=24):
        """
        input:
            self: pyMCDSts class instance.

            focus: string
                specifies either cell variable or substrate

            interface: string; default jpeg
                ffmpeg cannot directly translate svg image into a move.
                the interface image format will be used to bridge the gap.
                this images, from which the movie will be generated, have to exist.
                they can be generated with the make_jpeg, make_png, or make_tiff
                function.

            frame_rate: integer; default 24
                specifies how many images per second will be used.

        output:
            mp4 move file in output_path directory.
                interface image files in output_path directory.
`               additionally, the function will return the movie path and filename.

        description:
            this function generates a movie from all interface image files
            found in the output_path directory.
        """
        s_opathfile = f'{self.output_path}/{focus}_{interface}.mp4'
        s_ipathfiles = f'{self.output_path}/{focus}_*.{interface}'

        # check files exist
        if not glob(s_ipathfiles): # list returned by glob (files matching wildcard) is empty
            raise Exception('make_gif : images not found in output directory for given focus string')
        
        # system call
        os.system(f'ffmpeg -r {frame_rate} -f image2 -i {s_ipathfiles} -vcodec libx264 -pix_fmt yuv420p -strict -2 -tune animation -crf 15 -acodec none {s_opathfile}')

        # output
        return(s_opathfile)