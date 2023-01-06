####
# title: test_timeseries.py
#
# language: python3
# author: bue
# date: 2022-10-15
# license: BSD 3-Clause
#
# description:
#   pytest unit test library for heat.py
#   + https://docs.pytest.org/
#
#   note:
#   assert actual == expected, message
#   == value equality
#   is reference equality
#   pytest.approx for real values
#####

# load library
import os
import pathlib
import pcDataLoader as pc
import time

# const
i_sleep = 7
s_path_2d = str(pathlib.Path(pc.__file__).parent.resolve()/'data_timeseries_2d')

# load physicell data time series
class TestPyMcdsTs(object):
    ''' test for pc.pyMCDSts data loader. '''
    mcds = pc.pyMCDSts(s_path_2d, verbose=False)

    ### ffmpeg command ###
    def test_ffmpeg(self):
        b_ok = (os.system('ffmpeg -version') == 0)
        if not b_ok:
            print('Warning @ pcDataLoader : ffmpeg is not installed!')
        assert b_ok

    ## making movies with jpeg as interface ##
    def test_mcds_make_movie(self, mcds=mcds):
        # initialize
        s_pathfile = f'{s_path_2d}/movie_jpeg.mp4'
        if os.path.exists(s_pathfile):
            os.remove(s_pathfile)
        # generate jpeg interface images
        mcds.make_jpeg()
        time.sleep(i_sleep)
        # generate movie
        mcds.make_movie(moviefile='movie_jpeg.mp4', )
        assert os.path.getsize(s_pathfile) > 0
        # clean up
        os.remove(s_pathfile)
        for s_file in os.listdir(s_path_2d):
            if s_file.endswith('.jpeg'):
                os.remove(f'{s_path_2d}/{s_file}')
