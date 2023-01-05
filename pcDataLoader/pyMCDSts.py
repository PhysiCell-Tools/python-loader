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
#     PhysiCell model output directory. pyMCDSts.py was first froked from
#     PhysiCell-Tools python-loader, where it was implemented as
#     pyMCDS_timeseries.py, then totally rewritten and further developed.
#     the make_image and make_movie functions are cloned from the PhysiCell
#     Makefile.
#
#     note on difference image magick convert and mogrify
#     + https://graphicsmagick-tools.narkive.com/9Sowc4HF/gm-tools-mogrify-vs-convert
#########


# load libraries
import os
import pathlib
import platform
from .pyMCDS import pyMCDS
import xml.etree.ElementTree as ET

# constants
es_resize = {'*0.svg','*1.svg','*2.svg','*3.svg','*4.svg','*5.svg','*6.svg','*7.svg','*8.svg','*9.svg'} # only those this will be resized
ls_glob = sorted(es_resize) + ['initial.svg','legend.svg']

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
            setting grap to False will use less memory and speed up processing.

        verbose: boole; default True
            setting verbose to False for less text output, while processing.

    output:
        mcdsts: pyMCDSts class instance
            this instance offers functions to process all stored time steps
            for a simulation. no data is fetched by initialization.

    description:
        pyMCDSts.__init__ generates a class instance and stores there
        the input parameters. no data is fetched at initialization.
        the instance offers functions to process all timesteps
        in the output_path directory.
    """
    def __init__(self, output_path='.', graph=True, microenv=True, verbose=True):
        self.output_path = output_path
        self.graph = graph
        self.microenv = microenv
        self.verbose = verbose


    ## LOAD DATA
    def get_xmlfile_list(self):
        """
        input:
            self: pyMCDSts class instance.

        output:
            xmlfile_list: list of strings
            alphanumerical sorted list of physicell xml output files.

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
            list of physicell output xml pathfiles.

        output:
            l_mcds: listof mcds objetcs

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
                print() # carriage retun

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


    def _handle_resize(self, resize_factor=1):
        """
        input:
            self: pyMCDSts class instance.

            resize_factor: floating point number; default 1
            to specify image maginfied or scale down.
            the resize parameter will in any case be adjusted,
            so that the resulting image's hight and width are
            integer divisable by 2. this is becasue of a
            ffmpeg constrain for generating a movie out of images.

        output:
            s_resize: string
            image magick command resize parameter setting.

        description:
            internal function returns image magick command
            resize parameter setting, which in any case, even when
            resize_factor is 1, will generate ffmpeg compatible images.
        """
        # extract information form svg and resize
        tree = ET.parse(f'{self.output_path}/initial.svg')
        root = tree.getroot()
        r_width = float(root.get('width')) * resize_factor
        r_height = float(root.get('height')) * resize_factor
        # movie treat
        r_width = int(round(r_width / 2)) * 2
        r_height = int(round(r_height / 2)) * 2
        # output
        s_resize = f"-resize '{r_width}!x{r_height}!'"
        return(s_resize)


    def make_gif(self, resize_factor=1, giffile='timeseries.gif'):
        """
        input:
            self: pyMCDSts class instance.

            giffile: string; default 'timeseries.gif'
            gif image filename.

            resize_factor: floating point number; defaulat 1
            to specify image maginfied or scale down.

        output:
            gif file in output_path directory.
`            additionally, function will retun the path and filename.

        description:
            this function generate a gif image from all snapshot svg files
            found in the output_path directory.
        """
        s_magick = self._handle_magick()
        s_resize = self._handle_resize(resize_factor=resize_factor)

        # generate gif
        # bue: use convert, mogrify will cause troubles here!
        s_opathfile = f'{self.output_path}/{giffile}'
        os.system(f'{s_magick}convert {s_resize} {self.output_path}/snapshot*.svg {s_opathfile}')

        # output
        return(s_opathfile)


    def make_jpeg(self, resize_factor=1):
        """
        input:
            self: pyMCDSts class instance.

            glob: string
            wildcard filename pattern.

            resize_factor: floating point number; default 1
            to specify image maginfied or scale down.
            the resize parameter will in any case be adjusted,
            so that the resulting image's hight and width are
            integer divisable by 2. this is becasue of a
            ffmpeg constrain for generating a movie out of images.

        output:
            jpeg files in output_path directory.

        description:
            this function generate jpeg image equivalents from all svg files
            found in the output_path directory.
            jpeg is by definition a lossy compressed image format.
            https://en.wikipedia.org/wiki/JPEG
        """
        # bue: use mogrify, convert might cause troubles here!
        s_magick = self._handle_magick()
        s_resize = self._handle_resize(resize_factor=resize_factor)
        for s_glob in ls_glob:
            if (len(set(pathlib.Path(self.output_path).glob(s_glob))) > 0):
                if (s_glob in es_resize):
                    os.system(f'{s_magick}mogrify {s_resize} -format jpeg {self.output_path}/{s_glob} &')
                else:
                    os.system(f'{s_magick}mogrify -format jpeg {self.output_path}/{s_glob} &')


    def make_png(self, resize_factor=1, addargs='-transparent white'):
        """
        input:
            self: pyMCDSts class instance.

            resize_factor: floating point number; defaulat 1
            to specify image maginfied or scale down.
            the resize parameter will in any case be adjusted,
            so that the resulting image's hight and width are
            integer divisable by 2. this is becasue of a
            ffmpeg constrain for generating a movie out of images.

            addargs: string; default '-transparent white'
            sting to additional image mackig parameters.
            by default alpha channel transparence is set to white.

        output:
            png files in output_path directory.

        description:
            this function generate png image equivalents from all svg files
            found in the output_path directory.
            png is by definition a lossless compressed image format.
            https://en.wikipedia.org/wiki/Portable_Network_Graphics
        """
        # bue: use mogrify, convert might cause troubles here!
        s_magick = self._handle_magick()
        s_resize = self._handle_resize(resize_factor=resize_factor)
        for s_glob in ls_glob:
            if (len(set(pathlib.Path(self.output_path).glob(s_glob))) > 0):
                if (s_glob in es_resize):
                    os.system(f'{s_magick}mogrify {s_resize} {addargs} -format png {self.output_path}/{s_glob} &')
                else:
                    os.system(f'{s_magick}mogrify {addargs} -format png {self.output_path}/{s_glob} &')


    def make_tiff(self, resize_factor=1):
        """
        input:
            self: pyMCDSts class instance.

            resize_factor: floating point number; defaulat 1
            to specify image maginfied or scale down.
            the resize parameter will in any case be adjusted,
            so that the resulting image's hight and width are
            integer divisable by 2. this is becasue of a
            ffmpeg constrain for generating a movie out of images.

        output:
            tiff files in output_path directory.

        decription:
            this function generate tiff image equivalents from all svg files
            found in the output_path directory.
            https://en.wikipedia.org/wiki/TIFF
        """
        # bue: use mogrify, convert might cause troubles here!
        s_magick = self._handle_magick()
        s_resize = self._handle_resize(resize_factor=resize_factor)
        for s_glob in ls_glob:
            if (len(set(pathlib.Path(self.output_path).glob(s_glob))) > 0):
                if (s_glob in es_resize):
                    os.system(f'{s_magick}mogrify {s_resize} -format tiff {self.output_path}/{s_glob} & ')
                else:
                    os.system(f'{s_magick}mogrify -format tiff {self.output_path}/{s_glob} & ')


    def make_movie(self, interface='jpeg', moviefile='movie.mp4', frame_rate=24):
        """
        input:
            self: pyMCDSts class instance.

            interface: string; default jpeg
            ffmpeg can not directely translate svg image into a move.
            the interface image format will be used to bridge the gap.
            this images, from which the movie will be gererated, have to exist.
            they can be generated with the make_jpeg, make_png, or make_tiff
            function.

            moviefile: sting; default 'movie.mp4'
            mp4 movie file name.

            frame_rate: integer; default 24
            specifies how many images per secound will be used.

        output:
            mp4 move file in output_path directory.
            interface image files in output_path directory.
`           additionally, function will retun the movie path and filename.

        description:
            this function generates a movie from all interface image files
            found in the output_path directory.
        """
        # generate movie
        s_opathfile = f'{self.output_path}/{moviefile}'
        os.system(f'ffmpeg -r {frame_rate} -f image2 -i {self.output_path}/snapshot%08d.{interface} -vcodec libx264 -pix_fmt yuv420p -strict -2 -tune animation -crf 15 -acodec none {s_opathfile}')

        # output
        return(s_opathfile)

