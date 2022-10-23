########
#
#
#
########

# load libraries
import os
from pathlib import Path
from .pyMCDS import pyMCDS


# classes
class pyMCDS_ts:
    '''
    This class contains a np.array of pyMCDS objects as well as functions for
    extracting information from that list.

    Parameters
    ----------
    output_path : string
        String containing the path (relative or absolute) to the directory
        containing the PhysiCell output files

    Attributes
    ----------
    timeseries : array-like (pyMCDS) [n_timesteps,]
        Numpy array of pyMCDS objects sorted by time.
    '''
    def __init__(self, output_path='.', microenv=True, graph=True, verbose=True):
        self.output_path = output_path
        self.graph = graph
        self.microenv = microenv
        self.verbose = verbose

    ## LOAD DATA
    def get_xmlfile_list(self):
        '''
        '''
        # get a generator of output xml files sorted alphanumerically
        # bue 2022-10-22: is the output*.xml always the correct pattern?
        ls_pathfile = [o_pathfile.as_posix() for o_pathfile in sorted(Path(self.output_path).glob('output*.xml'))]

        return(ls_pathfile)

    def read_mcds(self, xmlfile_list=None):
        """
        Internal function. Does the actual work of initializing MultiCellDS by parsing the xml
        """
        # handle input
        if (xmlfile_list is None):
            xmlfile_list = self.get_xmlfile_list()

        # load mcds objects into list
        l_mcds = []
        for s_pathfile in xmlfile_list:
            mcds = pyMCDS(
                xml_file = s_pathfile,
                microenv = self.microenv,
                graph = self.graph,
                verbose = self.verbose
            )
            l_mcds.append(mcds)
            if self.verbose:
                print()

        # output
        return(l_mcds)

        #def make_movie(self, movie_file='movie.mp4'):
        """
        generates a movie from all svg files found in the PhysiCell output directory.

        Parameters
        ----------
        output_path: str, optional
            String containing the path (relative or absolute) to the directory
            where PhysiCell output image files are stored (default= "output/*.svg")

        movie_file: str, optional

        Returns
        -------
        mp4 movie_file move, made from the svg images.
        """
        # generate jpeg
        #os.system(f'identify -format "%h" $({self.output_path})/initial.svg > __H.txt')
        #os.system(f'identify -format "%w" $({self.output_path})/initial.svg > __W.txt')
        #os.system(f'expr 2 * \( $$(grep . __H.txt) / 2 \) > __H1.txt')
        #os.system(f'expr 2 * \( $$(grep . __W.txt) / 2 \) > __W1.txt')
        #os.system(f'echo "$$(grep . __W1.txt)!x$$(grep . __H1.txt)!" > __resize.txt')
        #os.system(f'mogrify -format jpg -resize $$(grep . __resize.txt) $({self.output_path})/snapshot*.svg')
        #os.system(f'rm -f __H*.txt __W*.txt __resize.txt')

        # generate mp4 from jpeg
        #s_opathfile = f'{self.output_path}/{movie_file}'
        #os.system(f'ffmpeg -r 24 -f image2 -i $({self.output_path})/snapshot%08d.jpg -vcodec libx264 -pix_fmt yuv420p -strict -2 -tune animation -crf 15 -acodec none $({s_opathfile})')
        #return(s_opathfile)

