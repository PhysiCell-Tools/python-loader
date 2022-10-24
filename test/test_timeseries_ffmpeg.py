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

# const
s_path_2d = str(pathlib.Path(pc.__file__).parent.resolve()/'data_timeseries_2d')

# load physicell data timeseries
class TestPyMcdsTs(object):
    ''' test for pc.pyMCDStimeseries data loader. '''
    mcds = pc.pyMCDSts(s_path_2d, verbose=False)

    ## movie ##
    def test_mcds_make_movie_jpeg(self, mcds=mcds):
        s_pathfile = f'{s_path_2d}/movie_jpeg.mp4'
        if os.path.exists(s_pathfile):
            os.remove(s_pathfile)
        mcds.make_movie(moviefile='movie_jpeg.mp4')
        assert os.path.getsize(s_pathfile) > 0
        os.remove(s_pathfile)
        for s_file in os.listdir(s_path_2d):
            if s_file.endswith('.jpeg'):
                os.remove(f'{s_path_2d}/{s_file}')

    def test_mcds_make_movie_png(self, mcds=mcds):
        s_pathfile = f'{s_path_2d}/movie_png.mp4'
        if os.path.exists(s_pathfile):
            os.remove(s_pathfile)
        mcds.make_movie(moviefile='movie_png.mp4', interface='png')
        assert os.path.getsize(s_pathfile) > 0
        os.remove(s_pathfile)
        for s_file in os.listdir(s_path_2d):
            if s_file.endswith('.png'):
                os.remove(f'{s_path_2d}/{s_file}')

    def test_mcds_make_movie_tiff(self, mcds=mcds):
        s_pathfile = f'{s_path_2d}/movie_tiff.mp4'
        if os.path.exists(s_pathfile):
            os.remove(s_pathfile)
        mcds.make_movie(moviefile='movie_tiff.mp4', interface='tiff')
        assert os.path.getsize(s_pathfile) > 0
        os.remove(s_pathfile)
        for s_file in os.listdir(s_path_2d):
            if s_file.endswith('.tiff'):
                os.remove(f'{s_path_2d}/{s_file}')

