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
import pcdl
import platform
import shutil

# const
s_path_2d = str(pathlib.Path(pcdl.__file__).parent.resolve()/'data_timeseries_2d')

# load physicell data time series
class TestPyMcdsTs(object):
    ''' test for pcdl.pyMCDSts data loader. '''
    mcdsts = pcdl.pyMCDSts(s_path_2d, verbose=False)

    ## get_xmlfile and read_mcds command ##
    def test_mcdsts_get_xmlfile_list(self, mcdsts=mcdsts):
        ls_xmlfile = mcdsts.get_xmlfile_list()
        assert len(ls_xmlfile) == 25

    def test_mcdsts_get_xmlfile_list_read_mcds(self, mcdsts=mcdsts):
        ls_xmlfile = mcdsts.get_xmlfile_list()
        ls_xmlfile = ls_xmlfile[-3:]
        ls_mcds = mcdsts.read_mcds(ls_xmlfile)
        assert len(ls_xmlfile) == 3 and \
               len(ls_mcds) == 3 and \
               ls_mcds[2].get_time() == 1440

    def test_mcdsts_read_mcds(self, mcdsts=mcdsts):
        ls_mcds = mcdsts.read_mcds()
        assert len(ls_mcds) == 25 and \
               ls_mcds[-1].get_time() == 1440


    ## magick command ##
    def test_mcdsts_handle_magick(self, mcdsts=mcdsts):
        s_magick = mcdsts._handle_magick()
        if not((os.system('magick --version') == 0) or ((platform.system() in {'Linux'}) and (os.system('convert --version') == 0))):
            s_magick = None
            print('Error @ pyMCDSts._handle_magick : image magick installation version >= 7.0 missing!')
        assert s_magick in {'', 'magick '}


    ## make_imgcell command ##
    def test_mcdsts_make_imgcell(self, mcdsts=mcdsts):
        s_path = mcdsts.make_imgcell(
            focus='cell_type',
            z_slice = -3.333,   # test if
            extrema = None,  # test if
            #cmap = 'viridis',  # matplotlib
            #grid = True,  # matplotlib
            xlim = None,  # test if
            ylim = None,  # test if
            #s = 36,  # matplotlib
            figsizepx = [641, 481],  # test non even pixel number
            ext = 'jpeg',
            figbgcolor = None,  # test if
        )
        assert os.path.exists(s_path + 'cell_type_000000000.0.jpeg') and \
               os.path.exists(s_path + 'cell_type_000001440.0.jpeg')
        shutil.rmtree(s_path)

    ## make_imgsubs command ##
    def test_mcdsts_make_imgsubs(self, mcdsts=mcdsts):
        s_path = mcdsts.make_imgsubs(
            focus = 'oxygen',
            z_slice = -3.333,  # test if
            extrema = None,  # test if
            #alpha = 1,  # matplotlib
            #fill = True,  # mcds.get_contour
            #cmap = 'viridis',  # matplotlib
            #grid = True,  # matplotlib
            xlim = None,  # test if
            ylim = None,  # test if
            figsizepx = [641, 481],  # test non even pixel number
            ext = 'jpeg',
            figbgcolor = None,  # test if
        )
        assert os.path.exists(s_path + 'oxygen_000000000.0.jpeg') and \
               os.path.exists(s_path + 'oxygen_000001440.0.jpeg')
        shutil.rmtree(s_path)

    ## make_gif command ##
    def test_mcdsts_make_gif(self, mcdsts=mcdsts):
        s_path = mcdsts.make_imgcell()
        s_opathfile = mcdsts.make_gif(
            path = s_path,
            #interface = 'jpeg',
        )
        assert os.path.exists(s_opathfile) and \
            (s_opathfile == s_path+'data_timeseries_2dcell_cell_type_z0_jpeg.gif')
        #shutil.rmtree(s_path)

    ## make_movie command ##
    def test_mcdsts_make_movie(self, mcdsts=mcdsts):
        s_path = mcdsts.make_imgcell()
        s_opathfile = mcdsts.make_movie(
            path = s_path,
            #interface = 'jpeg',
            #framerate = 12,
        )
        assert os.path.exists(s_opathfile) and \
            (s_opathfile == s_path+'data_timeseries_2dcell_cell_type_z0_jpeg12.mp4')
        shutil.rmtree(s_path)

