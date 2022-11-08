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
#     PhysiCell model output folder. pyMCDSts.py was first froked from
#     PhysiCell-Tools python-loader, where it was implemented as 
#     pyMCDS_timeseries.py, then totally rewritten and further developed. 
#     the make_image and  make_movie functions are cloned from the PhysiCell 
#     Makefile.
#########

# load libraries
import os
import pathlib
import platform
from .pyMCDS import pyMCDS
import xml.etree.ElementTree as ET

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
            for an simulation. no data is fetched by initialization.
                              
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

    def _handle_resize(self, resize_factor=1, movie=False):
        """
        input:
            self: pyMCDSts class instance.

            resize_factor: floating point number; defaulat 1
            to specify image maginfied or scale down. 

            movie: boolean; default False
            if set to True, the resize parameter will be adjusted,
            so that the resulting image's hight and width are 
            integer divisable by 2. this is an ffmpeg constrain for 
            generating a movie out of images.

        output:
            s_resize: string
            image magick command resize parameter setting.

        description:
            internal function returns a working image magick command
            resize parameter setting.
        """
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

    # BUE HERE I AM
    def make_gif(self, giffile='timeseries.gif', resize_factor=1):
        """
        input:
            self: pyMCDSts class instance.

            giffile: string; default 'timeseries.gif'


            resize_factor: floating point number; defaulat 1
            to specify image maginfied or scale down. 

        output:
            gif file in output_path folder.
            s_opathfile

        description:

        gif
        """
        s_magick = self._handle_magick()
        s_resize = self._handle_resize(resize_factor=resize_factor, movie=False)

        # generate gif
        s_opathfile = f'{self.output_path}/{giffile}'
        os.system(f'{s_magick}convert {s_resize} {self.output_path}/snapshot*.svg {s_opathfile}')

        # output
        return(s_opathfile)

    def make_jpeg(self, resize_factor=1, movie=False):
        """
        input:
            self: pyMCDSts class instance.

            resize_factor: floating point number; defaulat 1
            to specify image maginfied or scale down. 

            movie: boolean; default False
            if set to True, the resize parameter will be adjusted,
            so that the resulting image's hight and width are 
            integer divisable by 2. this is an ffmpeg constrain for 
            generating a movie out of images.

        output:
            jpeg files in output_path folder.

        description:

        jpeg
        """
        s_magick = self._handle_magick()
        s_resize = self._handle_resize(resize_factor=resize_factor, movie=movie)
        os.system(f'{s_magick}mogrify {s_resize} -format jpeg {self.output_path}/*.svg')

    def make_png(self, resize_factor=1, addargs='-transparent white', movie=False):
        """
        input:
            self: pyMCDSts class instance.

            resize_factor: floating point number; defaulat 1
            to specify image maginfied or scale down. 

            addargs: string; default '-transparent white'
            

            movie: boolean; default False
            if set to True, the resize parameter will be adjusted,
            so that the resulting image's hight and width are 
            integer divisable by 2. this is an ffmpeg constrain for 
            generating a movie out of images.

        output:
            png files in output_path folder.

        description:

        png
        """
        s_magick = self._handle_magick()
        s_resize = self._handle_resize(resize_factor=resize_factor, movie=movie)
        os.system(f'{s_magick}mogrify {s_resize} {addargs} -format png {self.output_path}/*.svg')

    def make_tiff(self, resize_factor=1, movie=False):
        """
        input:
            self: pyMCDSts class instance.

            resize_factor: floating point number; defaulat 1
            to specify image maginfied or scale down. 

            movie: boolean; default False
            if set to True, the resize parameter will be adjusted,
            so that the resulting image's hight and width are 
            integer divisable by 2. this is an ffmpeg constrain for 
            generating a movie out of images.

        output:
            tiff files in output_path folder.

        decription:

        tiff
        """
        s_magick = self._handle_magick()
        s_resize = self._handle_resize(resize_factor=resize_factor, movie=movie)
        os.system(f'{s_magick}mogrify {s_resize} -format tiff {self.output_path}/*.svg')

    def make_movie(self, moviefile='movie.mp4', frame_rate=24, resize_factor=1, interface='jpeg'):
        """
        input:
            self: pyMCDSts class instance.

            moviefile: sting; default 'movie.mp4'

            frame_rate: integer; default 24

            resize_factor: floating point number; defaulat 1
            to specify image maginfied or scale down. 

            interface: string; default 'jpeg'

        output:
            mp4 move file in output_path.
            s_opathfile

        description:


        generates a movie from all svg files found in the PhysiCell output directory.
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

