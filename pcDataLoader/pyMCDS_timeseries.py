########
#
#
#
########

# load libraries
from pathlib import Path
from .pyMCDS import pyMCDS


# classes
class pyMCDS_timeseries:
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
    def __init__(self, output_path='.', microenv=True, verbose=True):
        self.microenv = microenv
        self.verbose = verbose
        self.data = self._read_mcds(output_path=output_path)

    ## LOAD DATA
    def _read_mcds(self, output_path='.', microenv=True, verbose=True):
        """
        Internal function. Does the actual work of initializing MultiCellDS by parsing the xml
        """
        # get a generator of output xml files sorted alphanumerically
        # bue 2022-10-22: is the output*.xml always the correct pattern?
        ls_pathfile = [o_pathfile.as_posix() for o_pathfile in sorted(Path(output_path).glob('output*.xml'))]
        
        # load mcds objects into list
        l_mcds = []
        for s_pathfile in ls_pathfile;
            mcds = pyMCDS(xml_file=s_pathfile, microenv=True, verbose=True)
            l_mcds.append(mcds)

        # output
        return(l_mcds)


class pyMCDS_movie:
    """
    This class generates a movie from all svg files found in the PhysiCell output directory.

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
    def __init__(self, output_path='.', movie_file='movie.mp4'):
        self.data = self._making_movie(output_path, movie_file)

    def _making_movie(self, output_path='.',  movie_file='movie.mp4'):
        """ gererate move """
        s_ipathfile = f'{output_path}/*.svg'
        s_opathfile = f'{output_path}/{movie_file}'
        os.system(f"ffmpeg -f image2 -framerate 2/1 -pattern_type glob -i '{s_ipathfile}' -c:v libx264 -pix_fmt yuv420p {s_opathfile}")

