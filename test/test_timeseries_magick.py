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
import platform
import time

# const
i_sleep = 6
s_path_2d = str(pathlib.Path(pc.__file__).parent.resolve()/'data_timeseries_2d')

# load physicell data timeseries
class TestPyMcdsTs(object):
    ''' test for pc.pyMCDStimeseries data loader. '''
    mcds = pc.pyMCDSts(s_path_2d, verbose=False)

    ## magick command ##
    def test_mcds_handle_magick(self, mcds=mcds):
        s_magick = mcds._handle_magick()
        if not((os.system('magick --version') == 0) or ((platform.system() in {'Linux'}) and (os.system('convert --version') == 0))):
            s_magick = None
            print('Error @ pyMCDSts._handle_magick : image magick installation version >= 7.0 missing!')
        assert s_magick in {'', 'magick '}

    ## gif ##
    def test_mcds_make_gif(self, mcds=mcds):
        s_pathfile = f'{s_path_2d}/timeseries.gif'
        if os.path.exists(s_pathfile):
            os.remove(s_pathfile)
        mcds.make_gif()
        time.sleep(i_sleep)
        assert os.path.exists(s_pathfile)
        os.remove(s_pathfile)

    def test_mcds_make_gif_03(self, mcds=mcds):
        s_pathfile = f'{s_path_2d}/timeseries.gif'
        if os.path.exists(s_pathfile):
            os.remove(s_pathfile)
        mcds.make_gif(resize_factor=1/3)
        time.sleep(i_sleep)
        assert os.path.exists(s_pathfile)
        os.remove(s_pathfile)

    ## jpeg ##
    def test_mcds_make_jpeg(self, mcds=mcds):
        s_pathfile = f'{s_path_2d}/initial.jpeg'
        if os.path.exists(s_pathfile):
            os.remove(s_pathfile)
        mcds.make_jpeg()
        time.sleep(i_sleep)
        assert os.path.isfile(s_pathfile)
        for s_file in os.listdir(s_path_2d):
            if s_file.endswith('.jpeg'):
                os.remove(f'{s_path_2d}/{s_file}')

    def test_mcds_make_jpeg_03(self, mcds=mcds):
        s_pathfile = f'{s_path_2d}/initial.jpeg'
        if os.path.exists(s_pathfile):
            os.remove(s_pathfile)
        mcds.make_jpeg(resize_factor=1/3)
        time.sleep(i_sleep)
        assert os.path.exists(s_pathfile)
        for s_file in os.listdir(s_path_2d):
            if s_file.endswith('.jpeg'):
                os.remove(f'{s_path_2d}/{s_file}')

    ## png ##
    def test_mcds_make_png(self, mcds=mcds):
        s_pathfile = f'{s_path_2d}/initial.png'
        if os.path.exists(s_pathfile):
            os.remove(s_pathfile)
        mcds.make_png()
        time.sleep(i_sleep)
        assert os.path.exists(s_pathfile)
        for s_file in os.listdir(s_path_2d):
            if s_file.endswith('.png'):
                os.remove(f'{s_path_2d}/{s_file}')

    def test_mcds_make_png_03(self, mcds=mcds):
        s_pathfile = f'{s_path_2d}/initial.png'
        if os.path.exists(s_pathfile):
            os.remove(s_pathfile)
        mcds.make_png(resize_factor=1/3)
        time.sleep(i_sleep)
        assert os.path.exists(s_pathfile)
        for s_file in os.listdir(s_path_2d):
            if s_file.endswith('.png'):
                os.remove(f'{s_path_2d}/{s_file}')

    ## tiff ##
    def test_mcds_make_tiff(self, mcds=mcds):
        s_pathfile = f'{s_path_2d}/initial.tiff'
        if os.path.exists(s_pathfile):
            os.remove(s_pathfile)
        mcds.make_tiff()
        time.sleep(i_sleep)
        assert os.path.exists(s_pathfile)
        for s_file in os.listdir(s_path_2d):
            if s_file.endswith('.tiff'):
                os.remove(f'{s_path_2d}/{s_file}')

    def test_mcds_make_tiff_03(self, mcds=mcds):
        s_pathfile = f'{s_path_2d}/initial.tiff'
        if os.path.exists(s_pathfile):
            os.remove(s_pathfile)
        mcds.make_tiff(resize_factor=1/3)
        time.sleep(i_sleep)
        assert os.path.exists(s_pathfile)
        for s_file in os.listdir(s_path_2d):
            if s_file.endswith('.tiff'):
                os.remove(f'{s_path_2d}/{s_file}')
