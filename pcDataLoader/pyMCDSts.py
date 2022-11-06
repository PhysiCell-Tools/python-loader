#########
# title: pyMCDSts.py
#
# language: python3
# date: 2022-08-22
# license: BSD-3-Clause
# authors: Patrick Wall, Randy Heiland, Paul Macklin, Elmar Bucher
#
# description:
#     pyMCDSts.py defineds an object class, able to load and access
#     within python a time series of mcds objects loaded form a single
#     PhysiCell model output folder. pyMCDSts.py was a froked from
#     PhysiCell-Tools python-loader as pyMCDS_timeseries.py, then
#     totally rewritten and further developed. the make_image and
#     make_movie functions are cloned from PhysiCell Makefile.
#########

# load libraries
import os
import pathlib
import platform
from .pyMCDS import pyMCDS
import xml.etree.ElementTree as ET

# classes
class pyMCDSts:
    '''
    input:
    output_path : string
        String containing the path (relative or absolute) to the directory
        containing the PhysiCell output files

    output:
    timeseries : array-like (pyMCDS) [n_timesteps,]
        Numpy array of pyMCDS objects sorted by time.

    description:
        This class contains a np.array of pyMCDS objects as well as functions for
        extracting information from that list.
    '''
    def __init__(self, output_path='.', microenv=True, graph=True, verbose=True):
        self.output_path = output_path
        self.graph = graph
        self.microenv = microenv
        self.verbose = verbose


    ## LOAD DATA
    def get_xmlfile_list(self):
        '''
        input:
            self: pyMCDSts class instance.

        '''
        # get a generator of output xml files sorted alphanumerically
        # bue 2022-10-22: is the output*.xml always the correct pattern?
        ls_pathfile = [o_pathfile.as_posix() for o_pathfile in sorted(pathlib.Path(self.output_path).glob('output*.xml'))]

        return(ls_pathfile)

    def read_mcds(self, xmlfile_list=None):
        """
        input:
            self: pyMCDSts class instance.

        Internal function. Does the actual work of initializing MultiCellDS by parsing the xml
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
                print()

        # output
        return(l_mcds)


    ## TRANSFORM SVG
    def _handle_magick(self):
        '''
        input:
            self: pyMCDSts class instance.

        '''
        s_magick = 'magick '
        if (platform.system() in {'Linux'}) and (os.system('magick --version') != 0) and (os.system('convert --version') == 0):
            s_magick = ''
        return(s_magick)

    def _handle_resize(self, resize_factor=1, movie=False):
        '''
        input:
            self: pyMCDSts class instance.
        '''
        s_resize = ''
        if movie or (resize_factor != 1):
            # extract information form svg
            tree = ET.parse(f'{self.output_path}/initial.svg')
            root = tree.getroot()
            r_width = float(root.get('width'))
            r_hight = float(root.get('height'))
            if movie:
                r_width = int(r_width / 2) * 2
                r_hight = int(r_hight / 2) * 2
            s_resize = f"-resize '{r_width * resize_factor}!x{r_hight * resize_factor}!'"
        return(s_resize)

    def make_gif(self, giffile='timeseries.gif', resize_factor=1):
        '''
        input:
            self: pyMCDSts class instance.

        gif
        '''
        s_magick = self._handle_magick()
        s_resize = self._handle_resize(resize_factor=resize_factor, movie=False)

        # generate gif
        s_opathfile = f'{self.output_path}/{giffile}'
        os.system(f'{s_magick}convert {s_resize} {self.output_path}/snapshot*.svg {s_opathfile}')

        # output
        return(s_opathfile)

    def make_jpeg(self, resize_factor=1, movie=False):
        '''
        input:
            self: pyMCDSts class instance.

        jpeg
        '''
        s_magick = self._handle_magick()
        s_resize = self._handle_resize(resize_factor=resize_factor, movie=movie)
        os.system(f'{s_magick}mogrify {s_resize} -format jpeg {self.output_path}/*.svg')

    def make_png(self, resize_factor=1, addargs='-transparent white', movie=False):
        '''
        input:
            self: pyMCDSts class instance.

        png
        '''
        s_magick = self._handle_magick()
        s_resize = self._handle_resize(resize_factor=resize_factor, movie=movie)
        os.system(f'{s_magick}mogrify {s_resize} {addargs} -format png {self.output_path}/*.svg')

    def make_tiff(self, resize_factor=1, movie=False):
        '''
        input:
            self: pyMCDSts class instance.

        tiff
        '''
        s_magick = self._handle_magick()
        s_resize = self._handle_resize(resize_factor=resize_factor, movie=movie)
        os.system(f'{s_magick}mogrify {s_resize} -format tiff {self.output_path}/*.svg')

    def make_movie(self, moviefile='movie.mp4', frame_rate=24, resize_factor=1, interface='jpeg'):
        """
        input:
            self: pyMCDSts class instance.


        generates a movie from all svg files found in the PhysiCell output directory.

        Parameters
        ----------
        output_path: str, optional
            String containing the path (relative or absolute) to the directory
            where PhysiCell output image files are stored (default= "output/*.svg")

        moviefile: str, optional

        Returns
        -------
        mp4 moviefile move, made from the svg images.
        """
        # gererate interface images
        if interface in {'JPEG', 'jpeg', 'jpg', 'jpe'}:
            interface = 'jpeg'
            self.make_jpeg(resize_factor=resize_factor, movie=True)

        elif interface in {'PNG', 'png'}:
            interface = 'png'
            self.make_png(resize_factor=resize_factor, addargs='', movie=True)

        elif interface in {'TIFF', 'tiff', 'tif'}:
            interface = 'tiff'
            self.make_tiff(resize_factor=resize_factor, movie=True)

        else:
            sys.exit(f'Error @ pyMCDSts.make_movie : unknown interface format {interface}.\nknoen are jpeg, png, and tiff.')

        # generate movie
        s_opathfile = f'{self.output_path}/{moviefile}'
        os.system(f'ffmpeg -r {frame_rate} -f image2 -i {self.output_path}/snapshot%08d.{interface} -vcodec libx264 -pix_fmt yuv420p -strict -2 -tune animation -crf 15 -acodec none {s_opathfile}')

        # output
        return(s_opathfile)

