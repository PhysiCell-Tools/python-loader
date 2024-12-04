####
# title: test_timeseries_2d.py
#
# language: python3
# author: Elmar Bucher
# date: 2022-10-15
# license: BSD 3-Clause
#
# description:
#   pytest unit test library for the pcdl library pyMCDSts class.
#   + https://docs.pytest.org/
#
#   note:
#   assert actual == expected, message
#   == value equality
#   is reference equality
#   pytest.approx for real values
#####


# load library
import matplotlib.pyplot as plt
import numpy  as np
import os
import pathlib
import pcdl
import shutil


# const
s_path_2d = str(pathlib.Path(pcdl.__file__).parent.resolve()/'data_timeseries_2d')


## download test data ##
if not os.path.exists(s_path_2d):
    pcdl.install_data()


## making movies related functions ##

class TestPyMcdsTsMovies(object):
    ''' tests for loading a pcdl.pyMCDS data set. '''
    mcdsts = pcdl.pyMCDSts(s_path_2d, verbose=True)

    ## make_gif and magick ommand ##
    def test_mcdsts_make_gif_jpeg(self, mcdsts=mcdsts):
        ls_opathfile = mcdsts.plot_scatter()
        s_opath = '/'.join(ls_opathfile[0].replace('\\','/').split('/')[:-1])
        s_opathfile = mcdsts.make_gif(
            path = s_opath,
            #interface = 'jpeg',
        )
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (os.path.exists(s_opathfile)) and \
              (s_opathfile.endswith('pcdl/data_timeseries_2d/cell_cell_type_z0.0/cell_cell_type_z0.0_jpeg.gif'))
        #os.remove(s_opathfile)
        shutil.rmtree(s_opath)

    def test_mcdsts_make_gif_tiff(self, mcdsts=mcdsts):
        ls_opathfile = mcdsts.plot_scatter(ext='tiff')
        s_opath = '/'.join(ls_opathfile[0].replace('\\','/').split('/')[:-1])
        s_opathfile = mcdsts.make_gif(
            path = s_opath,
            interface = 'tiff',
        )
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (os.path.exists(s_opathfile)) and \
              (s_opathfile.endswith('pcdl/data_timeseries_2d/cell_cell_type_z0.0/cell_cell_type_z0.0_tiff.gif'))
        #os.remove(s_opathfile)
        shutil.rmtree(s_opath)

    ## make_movie and magick command ##
    def test_mcdsts_make_movie_jpeg12(self, mcdsts=mcdsts):
        ls_opathfile = mcdsts.plot_scatter()
        s_opath = '/'.join(ls_opathfile[0].replace('\\','/').split('/')[:-1])
        s_opathfile = mcdsts.make_movie(
            path = s_opath,
            #interface = 'jpeg',
            #framerate = 12,
        )
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (os.path.exists(s_opathfile)) and \
              (s_opathfile.endswith('pcdl/data_timeseries_2d/cell_cell_type_z0.0/cell_cell_type_z0.0_jpeg12.mp4'))
        #os.remove(s_opathfile)
        shutil.rmtree(s_opath)

    def test_mcdsts_make_movie_tiff12(self, mcdsts=mcdsts):
        ls_opathfile = mcdsts.plot_scatter(ext='tiff')
        s_opath = '/'.join(ls_opathfile[0].replace('\\','/').split('/')[:-1])
        s_opathfile = mcdsts.make_movie(
            path = s_opath,
            interface = 'tiff',
            #framerate = 12,
        )
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (os.path.exists(s_opathfile)) and \
              (s_opathfile.endswith('pcdl/data_timeseries_2d/cell_cell_type_z0.0/cell_cell_type_z0.0_tiff12.mp4'))
        #os.remove(s_opathfile)
        shutil.rmtree(s_opath)

    def test_mcdsts_make_movie_jpeg6(self, mcdsts=mcdsts):
        ls_opathfile = mcdsts.plot_scatter()
        s_opath = '/'.join(ls_opathfile[0].replace('\\','/').split('/')[:-1])
        s_opathfile = mcdsts.make_movie(
            path = s_opath,
            #interface = 'jpeg',
            framerate = 6,
        )
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (os.path.exists(s_opathfile)) and \
              (s_opathfile.endswith('pcdl/data_timeseries_2d/cell_cell_type_z0.0/cell_cell_type_z0.0_jpeg6.mp4'))
        #os.remove(s_opathfile)
        shutil.rmtree(s_opath)


## data loading related functions ##

class TestPyMcdsTsInit(object):
    ''' tests for loading a pcdl.pyMCDSts data set. '''

    def test_mcdsts_set_verbose_true(self):
        mcdsts = pcdl.pyMCDSts(s_path_2d, load=False, verbose=False)
        mcdsts.set_verbose_true()
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (mcdsts.verbose)

    def test_mcdsts_set_verbose_false(self):
        mcdsts = pcdl.pyMCDSts(s_path_2d, load=False, verbose=True)
        mcdsts.set_verbose_false()
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (not mcdsts.verbose)

    ## get_xmlfile and read_mcds command and get_mcds_list ##
    def test_mcdsts_get_xmlfile_list(self):
        mcdsts = pcdl.pyMCDSts(s_path_2d, load=False, verbose=True)
        ls_xmlfile = mcdsts.get_xmlfile_list()
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (ls_xmlfile[0] == 'output00000000.xml') and \
              (ls_xmlfile[-1] == 'output00000024.xml') and \
              (len(ls_xmlfile) == 25)

    def test_mcdsts_get_mcds_list(self):
        mcdsts = pcdl.pyMCDSts(s_path_2d, load=True, verbose=True)
        l_mcds = mcdsts.get_mcds_list()
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(mcdsts.l_mcds[0])) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(mcdsts.l_mcds[-1])) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (mcdsts.l_mcds[0].get_time() == 0) and \
              (mcdsts.l_mcds[-1].get_time() == 1440) and \
              (len(mcdsts.l_mcds) == 25) and \
              (mcdsts.l_mcds == l_mcds)

    def test_mcdsts_read_mcds(self):
        mcdsts = pcdl.pyMCDSts(s_path_2d, load=False, verbose=True)
        l_mcds_loadfalse  = mcdsts.get_mcds_list()
        mcdsts.read_mcds()
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(mcdsts.l_mcds[0])) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(mcdsts.l_mcds[-1])) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (mcdsts.l_mcds[0].get_time() == 0) and \
              (mcdsts.l_mcds[-1].get_time() == 1440) and \
              (len(mcdsts.l_mcds) == 25) and \
              (l_mcds_loadfalse is None)

    def test_mcdsts_read_mcds_xmlfilelist(self):
        mcdsts = pcdl.pyMCDSts(s_path_2d, load=False, verbose=True)
        ls_xmlfile = mcdsts.get_xmlfile_list()
        ls_xmlfile = ls_xmlfile[-3:]
        l_mcds = mcdsts.read_mcds(ls_xmlfile)
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(mcdsts.l_mcds[0])) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(mcdsts.l_mcds[-1])) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (mcdsts.l_mcds[0].get_time() == 1320) and \
              (mcdsts.l_mcds[-1].get_time() == 1440) and \
              (len(ls_xmlfile) == 3) and \
              (len(mcdsts.l_mcds) == 3) and \
              (mcdsts.l_mcds == l_mcds)


## micro environment related functions ##

class TestPyMcdsTsMicroenv(object):
    ''' tests for pcdl.pyMCDS micro environment related functions. '''
    mcdsts = pcdl.pyMCDSts(s_path_2d, verbose=True)

    def test_mcdsts_get_conc_df(self, mcdsts=mcdsts):
        ldf_conc = mcdsts.get_conc_df(values=2, drop=set(), keep=set(), collapse=False)
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(ldf_conc)) == "<class 'list'>") and \
              (str(type(ldf_conc[0])) == "<class 'pandas.core.frame.DataFrame'>") and \
              (ldf_conc[0].shape == (121, 9)) and \
              (ldf_conc[-1].shape == (121, 10)) and \
              (len(ldf_conc) == 25)

    def test_mcdsts_get_conc_df_collapse(self, mcdsts=mcdsts):
        df_conc = mcdsts.get_conc_df(values=2, drop=set(), keep=set(), collapse=True)
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(df_conc)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_conc.shape == (3025, 10))

    def test_mcdsts_get_conc_attribute(self, mcdsts=mcdsts):
        dl_conc = mcdsts.get_conc_attribute(values=1, drop=set(), keep=set(), allvalues=False)
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(dl_conc)) == "<class 'dict'>") and \
              (str(type(dl_conc['oxygen'])) == "<class 'list'>") and \
              (str(type(dl_conc['oxygen'][0])) == "<class 'float'>") and \
              (len(dl_conc.keys()) == 1) and \
              (len(dl_conc['oxygen']) == 2)

    def test_mcdsts_get_conc_attribute_values(self, mcdsts=mcdsts):
        dl_conc = mcdsts.get_conc_attribute(values=2, drop=set(), keep=set(), allvalues=False)
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(dl_conc)) == "<class 'dict'>") and \
              (str(type(dl_conc['oxygen'])) == "<class 'list'>") and \
              (str(type(dl_conc['oxygen'][0])) == "<class 'float'>") and \
              (len(dl_conc.keys()) == 1) and \
              (len(dl_conc['oxygen']) == 2)

    def test_mcdsts_get_conc_attribute_allvalues(self, mcdsts=mcdsts):
        dl_conc = mcdsts.get_conc_attribute(values=1, drop=set(), keep=set(), allvalues=True)
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(dl_conc)) == "<class 'dict'>") and \
              (str(type(dl_conc['oxygen'])) == "<class 'list'>") and \
              (str(type(dl_conc['oxygen'][0])) == "<class 'float'>") and \
              (len(dl_conc.keys()) == 1) and \
              (len(dl_conc['oxygen']) > 2)

    ## plot_contour command ##
    def test_mcdsts_plot_contour_if(self, mcdsts=mcdsts):
        ls_pathfile = mcdsts.plot_contour(
            focus = 'oxygen',
            z_slice = -3.333,  # test if
            extrema = None,  # test if and for loop
            #alpha = 1,  # pyMCD
            #fill = True,  # pyMCD
            #cmap = 'viridis',  # pyMCD
            #title = '',  # test default
            #grid = True,  # pyMCD
            xlim = None,  # test if
            ylim = None,  # test if
            #xyequal = True,  # pyMCD
            figsizepx = None,  # test if
            ext = 'jpeg',  # test file case
            figbgcolor = None,  # test if
        )
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (ls_pathfile[0].replace('\\','/').endswith('/pcdl/data_timeseries_2d/conc_oxygen_z0.0/output00000000_oxygen.jpeg')) and \
              (os.path.exists(ls_pathfile[0])) and \
              (os.path.getsize(ls_pathfile[0]) > 2**10) and \
              (ls_pathfile[-1].replace('\\','/').endswith('/pcdl/data_timeseries_2d/conc_oxygen_z0.0/output00000024_oxygen.jpeg')) and \
              (os.path.exists(ls_pathfile[-1])) and \
              (os.path.getsize(ls_pathfile[-1]) > 2**10) and \
              (len(ls_pathfile) == 25)
        for s_pathfile in ls_pathfile:
            os.remove(s_pathfile)

    def test_mcdsts_plot_contour_else(self, mcdsts=mcdsts):
        l_fig = mcdsts.plot_contour(
            focus = 'oxygen',
            z_slice = 0.0,  # jump over if
            extrema = [0, 38],  # jump over if
            #alpha = 1,  # pyMCDS
            #fill = True,  # pyMCDS
            #cmap = 'viridis',  # pyMCDS
            title = 'abc',  # test non default
            #grid = True,  # pyMCDS
            xlim = [-31, 301],  # jump over if
            ylim = [-21, 201],  # jump over if
            #xyequal = True,  # pyMCDS
            figsizepx = [641, 481],  # test non even pixel
            ext = None,  # test fig case
            figbgcolor = 'yellow',  # not a file
        )
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(l_fig[0])) == "<class 'matplotlib.figure.Figure'>") and \
              (str(type(l_fig[-1])) == "<class 'matplotlib.figure.Figure'>") and \
              (len(l_fig) == 25)
        plt.close()

    def test_mcdsts_make_conc_vtk(self, mcdsts=mcdsts):
        ls_pathfile = mcdsts.make_conc_vtk()
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (ls_pathfile[0].endswith('/pcdl/data_timeseries_2d/output00000000_conc.vtr')) and \
              (ls_pathfile[-1].endswith('/pcdl/data_timeseries_2d/output00000024_conc.vtr')) and \
              (os.path.exists(ls_pathfile[0])) and \
              (os.path.exists(ls_pathfile[-1])) and \
              (os.path.getsize(ls_pathfile[0]) > 2**10) and\
              (os.path.getsize(ls_pathfile[-1]) > 2**10) and\
              (len(ls_pathfile) == 25)
        for s_pathfile in ls_pathfile:
            os.remove(s_pathfile)


## cell related functions ##

class TestPyMcdsCell(object):
    ''' tests for pcdl.pyMCDS cell related functions. '''
    mcdsts = pcdl.pyMCDSts(s_path_2d, verbose=False)

    def test_mcdsts_get_cell_df(self, mcdsts=mcdsts):
        ldf_cell = mcdsts.get_cell_df(values=2, drop=set(), keep=set(), collapse=False)
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(ldf_cell)) == "<class 'list'>") and \
              (str(type(ldf_cell[0])) == "<class 'pandas.core.frame.DataFrame'>") and \
              (ldf_cell[0].shape == (889, 19)) and \
              (ldf_cell[-1].shape == (1099, 40)) and \
              (len(ldf_cell) == 25)

    def test_mcdsts_get_cell_df_collapse(self, mcdsts=mcdsts):
        df_cell = mcdsts.get_cell_df(values=2, drop=set(), keep=set(), collapse=True)
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_cell.shape == (24758, 41))

    def test_mcdsts_get_cell_attribute(self, mcdsts=mcdsts):
        dl_cell = mcdsts.get_cell_attribute(values=1, drop=set(), keep=set(), allvalues=False)
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(dl_cell)) == "<class 'dict'>") and \
              (str(type(dl_cell['dead'])) == "<class 'list'>") and \
              (str(type(dl_cell['dead'][0])) == "<class 'bool'>") and \
              (str(type(dl_cell['cell_count_voxel'])) == "<class 'list'>") and \
              (str(type(dl_cell['cell_count_voxel'][0])) == "<class 'int'>") and \
              (str(type(dl_cell['cell_density_micron3'])) == "<class 'list'>") and \
              (str(type(dl_cell['cell_density_micron3'][0])) == "<class 'float'>") and \
              (str(type(dl_cell['cell_type'])) == "<class 'list'>") and \
              (str(type(dl_cell['cell_type'][0])) == "<class 'str'>") and \
              (len(dl_cell.keys()) == 83) and \
              (len(dl_cell['dead']) == 2) and \
              (len(dl_cell['cell_count_voxel']) == 2) and \
              (len(dl_cell['cell_density_micron3']) == 2) and \
              (len(dl_cell['cell_type']) == 1)

    def test_mcdsts_get_cell_attribute_values(self, mcdsts=mcdsts):
        dl_cell = mcdsts.get_cell_attribute(values=2, drop=set(), keep=set(), allvalues=False)
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(dl_cell)) == "<class 'dict'>") and \
              (len(dl_cell.keys()) == 28)

    def test_mcdsts_get_cell_attribute_allvalues(self, mcdsts=mcdsts):
        dl_cell = mcdsts.get_cell_attribute(values=1, drop=set(), keep=set(), allvalues=True)
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(dl_cell)) == "<class 'dict'>") and \
              (str(type(dl_cell['dead'])) == "<class 'list'>") and \
              (str(type(dl_cell['dead'][0])) == "<class 'bool'>") and \
              (str(type(dl_cell['cell_count_voxel'])) == "<class 'list'>") and \
              (str(type(dl_cell['cell_count_voxel'][0])) == "<class 'int'>") and \
              (str(type(dl_cell['cell_density_micron3'])) == "<class 'list'>") and \
              (str(type(dl_cell['cell_density_micron3'][0])) == "<class 'float'>") and \
              (str(type(dl_cell['cell_type'])) == "<class 'list'>") and \
              (str(type(dl_cell['cell_type'][0])) == "<class 'str'>") and \
              (len(dl_cell.keys()) == 83) and \
              (len(dl_cell['dead']) == 2) and \
              (len(dl_cell['cell_count_voxel']) > 2) and \
              (len(dl_cell['cell_density_micron3']) > 2) and \
              (len(dl_cell['cell_type']) == 1)

    ## plot_scatter command ##
    def test_mcdsts_plot_scatter_num(self, mcdsts=mcdsts):
        ls_pathfile = mcdsts.plot_scatter(
            focus='pressure',  # case numeric
            z_slice = -3.333,   # test if
            z_axis = None,  # test iff numeric
            #alpha = 1,  # matplotlib
            #cmap = 'viridis',  # matplotlib
            #title = '', test default
            #grid = True,  # matplotlib
            #legend_loc='lower left',  # matplotlib
            xlim = None,  # test if
            ylim = None,  # test if
            #xyequal = True,  # pyMCDS
            s = None,  # test if
            figsizepx = None,  # case extract from initial.svg
            ext = 'jpeg',  # generate file case
            figbgcolor = None,  # test if
        )
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (ls_pathfile[0].replace('\\','/').endswith('/pcdl/data_timeseries_2d/cell_pressure_z0.0/output00000000_pressure.jpeg')) and \
              (os.path.exists(ls_pathfile[0])) and \
              (os.path.getsize(ls_pathfile[0]) > 2**10) and \
              (ls_pathfile[-1].replace('\\','/').endswith('/pcdl/data_timeseries_2d/cell_pressure_z0.0/output00000024_pressure.jpeg')) and \
              (os.path.exists(ls_pathfile[-1])) and \
              (os.path.getsize(ls_pathfile[-1]) > 2**10) and \
              (len(ls_pathfile) == 25)
        for s_pathfile in ls_pathfile:
            os.remove(s_pathfile)

    def test_mcdsts_plot_scatter_cat(self, mcdsts=mcdsts):
        l_fig = mcdsts.plot_scatter(
            focus='cell_type',  # case categorical
            z_slice = 0.0,   # jump over if
            z_axis = None,  # test iff  categorical
            #alpha = 1,  # pyMCDS
            #cmap = 'viridis',  # pyMCDS
            title = 'abc', # test non default
            #grid = True,  # pyMCDS
            #legend_loc='lower left',  # pyMCDS
            xlim = None,  # test if
            ylim = None,  # test if
            #xyequal = True,  # pyMCDS
            s = None,  # test if
            figsizepx = [641, 481],  # test case non even pixel number
            ext = None,  # test fig case
            figbgcolor = None,  # not a file
        )
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(l_fig[0])) == "<class 'matplotlib.figure.Figure'>") and \
              (str(type(l_fig[-1])) == "<class 'matplotlib.figure.Figure'>") and \
              (len(l_fig) == 25)
        plt.close()

    def test_mcdsts_make_cell_vtk(self, mcdsts=mcdsts):
        ls_pathfile = mcdsts.make_cell_vtk(visualize=False)
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (ls_pathfile[0].endswith('/pcdl/data_timeseries_2d/output00000000_cell.vtp')) and \
              (ls_pathfile[-1].endswith('/pcdl/data_timeseries_2d/output00000024_cell.vtp')) and \
              (os.path.exists(ls_pathfile[0])) and \
              (os.path.exists(ls_pathfile[-1])) and \
              (os.path.getsize(ls_pathfile[0]) > 2**10) and\
              (os.path.getsize(ls_pathfile[-1]) > 2**10) and\
              (len(ls_pathfile) == 25)
        for s_pathfile in ls_pathfile:
            os.remove(s_pathfile)


## graph related functions ##
class TestPyMcdsGraph(object):
    ''' tests for pcdl.pyMCDS graph related functions. '''
    mcdsts = pcdl.pyMCDSts(s_path_2d, verbose=False)

    ## graph related functions ##
    def test_mcdsts_get_graph_gml_attached_defaultattr(self, mcdsts=mcdsts):
        ls_pathfile = mcdsts.make_graph_gml(graph_type='attached', edge_attribute=True, node_attribute=[])
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (ls_pathfile[0].endswith('/pcdl/data_timeseries_2d/output00000000_attached.gml')) and \
              (ls_pathfile[-1].endswith('/pcdl/data_timeseries_2d/output00000024_attached.gml')) and \
              (os.path.exists(ls_pathfile[0])) and \
              (os.path.exists(ls_pathfile[-1])) and \
              (len(ls_pathfile) == 25)
        for s_pathfile in ls_pathfile:
            os.remove(s_pathfile)

    def test_mcdsts_get_graph_gml_neighbor_noneattr(self, mcdsts=mcdsts):
        ls_pathfile = mcdsts.make_graph_gml(graph_type='neighbor', edge_attribute=False, node_attribute=[])
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (ls_pathfile[0].endswith('/pcdl/data_timeseries_2d/output00000000_neighbor.gml')) and \
              (ls_pathfile[-1].endswith('/pcdl/data_timeseries_2d/output00000024_neighbor.gml')) and \
              (os.path.exists(ls_pathfile[0])) and \
              (os.path.exists(ls_pathfile[-1])) and \
              (len(ls_pathfile) == 25)
        for s_pathfile in ls_pathfile:
            os.remove(s_pathfile)

    def test_mcdsts_get_graph_gml_neighbor_allattr(self, mcdsts=mcdsts):
        ls_pathfile = mcdsts.make_graph_gml(graph_type='neighbor', edge_attribute=True, node_attribute=['dead','cell_count_voxel','cell_density_micron3','cell_type'])
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (ls_pathfile[0].endswith('/pcdl/data_timeseries_2d/output00000000_neighbor.gml')) and \
              (ls_pathfile[-1].endswith('/pcdl/data_timeseries_2d/output00000024_neighbor.gml')) and \
              (os.path.exists(ls_pathfile[0])) and \
              (os.path.exists(ls_pathfile[-1])) and \
              (len(ls_pathfile) == 25)
        for s_pathfile in ls_pathfile:
            os.remove(s_pathfile)

## ome tiff related functions ##
class TestPyMcdsOmeTiff(object):
    ''' tests for pcdl.pyMCDS ome tiff related functions. '''
    mcdsts = pcdl.pyMCDSts(s_path_2d, verbose=False)

    ## ome tiff related functions ##
    def test_mcdsts_make_ome_tiff_defaultattr_00(self, mcdsts=mcdsts):
        la_ometiff = mcdsts.make_ome_tiff(cell_attribute='ID', conc_cutoff={}, focus=None, file=False, collapse=False)
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (type(la_ometiff) is list) and \
              (type(la_ometiff[0]) is np.ndarray) and \
              (type(la_ometiff[-1]) is np.ndarray) and \
              (la_ometiff[0].dtype == np.float32) and \
              (la_ometiff[-1].dtype == np.float32) and \
              (la_ometiff[0].shape == (2, 1, 200, 300)) and \
              (la_ometiff[-1].shape ==  (2, 1, 200, 300)) and \
              (len(la_ometiff) == 25)

    def test_mcdsts_make_ome_tiff_defaultattr_01(self, mcdsts=mcdsts):
        a_ometiff = mcdsts.make_ome_tiff(cell_attribute='ID', conc_cutoff={}, focus=None, file=False, collapse=True)
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (type(a_ometiff) is np.ndarray) and \
              (a_ometiff.dtype == np.float32) and \
              (a_ometiff.shape == (25, 2, 1, 200, 300))

    def test_mcdsts_make_ome_tiff_defaultattr_10(self, mcdsts=mcdsts):
        ls_pathfile = mcdsts.make_ome_tiff(cell_attribute='ID', conc_cutoff={}, focus=None, file=True, collapse=False)
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (ls_pathfile[0].endswith('pcdl/data_timeseries_2d/output00000000_oxygen_cancer_cell_ID.ome.tiff')) and \
              (ls_pathfile[-1].endswith('pcdl/data_timeseries_2d/output00000024_oxygen_cancer_cell_ID.ome.tiff')) and \
              (os.path.exists(ls_pathfile[0])) and \
              (os.path.exists(ls_pathfile[-1])) and \
              (os.path.getsize(ls_pathfile[0]) > 2**10) and\
              (os.path.getsize(ls_pathfile[-1]) > 2**10) and\
              (len(ls_pathfile) == 25)
        for s_pathfile in ls_pathfile:
            os.remove(s_pathfile)

    def test_mcdsts_make_ome_tiff_defaultattr_11(self, mcdsts=mcdsts):
        s_pathfile = mcdsts.make_ome_tiff(cell_attribute='ID', conc_cutoff={}, focus=None, file=True, collapse=True)
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (s_pathfile.endswith('pcdl/data_timeseries_2d/timeseries_oxygen_cancer_cell_ID.ome.tiff')) and \
              (os.path.exists(s_pathfile)) and \
              (os.path.getsize(s_pathfile) > 2**10 )
        os.remove(s_pathfile)


## timeseries related functions ##

class TestPyMcdsTimeseries(object):
    ''' tests for pcdl.pyMCDS graph related functions. '''
    mcdsts = pcdl.pyMCDSts(s_path_2d, verbose=False)

    ## plot_timeseries command ##
    def test_mcdsts_plot_timeseries_none_none_none_cell_ax_jpeg(self, mcdsts=mcdsts):
        fig, ax = plt.subplots()
        s_pathfile = mcdsts.plot_timeseries(
            focus_cat = None,  # test if {None/total, 'cell_type'}
            focus_num = None,  # test if {None/count, 'oxygen'}
            #aggregate_num = np.mean,  # pandas
            frame = 'cell',  # test if else {'df_cell', 'df_conc'}
            z_slice = None,  # test timeseries
            #logy = False,  # pandas
            #ylim = None,  # pandas
            #secondary_y = None,  # pandas
            #subplots = False,  # pandas
            #sharex = False,  # pandas
            #sharey = False,  # pandas
            #linestyle = '-',  # pandas
            #linewidth = None,  # pandas
            #cmap = None,  # pandas
            #color = None,  # pandas
            #grid = True,  # pandas
            #legend = True, # pandas
            yunit = None,  # test if {None, 'mmHg'}
            #title = None, pandas
            ax = ax,  # test if else {None, ax}
            figsizepx = [641, 481],  # test non even pixel number
            ext = 'jpeg',  # test if else {'jpeg', None}
            figbgcolor = None  # test if
        )
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (s_pathfile.endswith('/pcdl/data_timeseries_2d/timeseries_cell_total_count.jpeg')) and \
              (os.path.exists(s_pathfile))
        os.remove(s_pathfile)

    def test_mcdsts_plot_timeseries_cat_none_yunit_cell(self, mcdsts=mcdsts):
        fig = mcdsts.plot_timeseries(
            focus_cat = 'cell_type',  # test if {None/total, 'cell_type'}
            focus_num = None,  # test if {None/count, 'oxygen'}
            #aggregate_num = np.mean,  # pandas
            frame = 'df_cell',  # test if else {'df_cell', 'df_conc'}
            z_slice = -0.3,  # test if if
            #logy = False,  # pandas
            #ylim = None,  # pandas
            #secondary_y = None,  # pandas
            #subplots = False,  # pandas
            #sharex = False,  # pandas
            #sharey = False,  # pandas
            #linestyle = '-',  # pandas
            #linewidth = None,  # pandas
            #cmap = None,  # pandas
            #color = None,  # pandas
            #grid = True,  # pandas
            #legend = True,
            yunit = 'mmHg',  # test if {None, 'mmHg'}
            #title = None, pandas
            ax = None,  # test if else {None, ax}
            figsizepx = [641, 481],  # test non even pixel number
            ext = None,  # test if else {'jpeg', None}
            figbgcolor = None  # test if
        )
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(fig)) == "<class 'matplotlib.figure.Figure'>")
        plt.close()

    def test_mcdsts_plot_timeseries_none_num_yunit_cell(self, mcdsts=mcdsts):
        fig = mcdsts.plot_timeseries(
            focus_cat = None,  # test if {None/total, 'cell_type'}
            focus_num = 'oxygen',  # test if {None/count, 'oxygen'}
            #aggregate_num = np.mean,  # pandas
            frame = 'df_cell',  # test if else {'df_cell', 'df_conc'}
            z_slice = -0.3,  # test if if
            #logy = False,  # pandas
            #ylim = None,  # pandas
            #secondary_y = None,  # pandas
            #subplots = False,  # pandas
            #sharex = False,  # pandas
            #sharey = False,  # pandas
            #linestyle = '-',  # pandas
            #linewidth = None,  # pandas
            #cmap = None,  # pandas
            #color = None,  # pandas
            #grid = True,  # pandas
            #legend = True,
            yunit = 'mmHg',  # test if {None, 'mmHg'}
            #title = None, pandas
            ax = None,  # test if else {None, ax}
            figsizepx = [641, 481],  # test non even pixel number
            ext = None,  # test if else {'jpeg', None}
            figbgcolor = None  # test if
        )
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(fig)) == "<class 'matplotlib.figure.Figure'>")
        plt.close()

    def test_mcdsts_plot_timeseries_cat_num_none_cell(self, mcdsts=mcdsts):
        fig = mcdsts.plot_timeseries(
            focus_cat = 'cell_type',  # test if {None/total, 'cell_type'}
            focus_num = 'oxygen',  # test if {None/count, 'oxygen'}
            #aggregate_num = np.mean,  # pandas
            frame = 'cell',  # test if else {'df_cell', 'df_conc'}
            z_slice = -0.3,  # test if if
            #logy = False,  # pandas
            #ylim = None,  # pandas
            #secondary_y = None,  # pandas
            #subplots = False,  # pandas
            #sharex = False,  # pandas
            #sharey = False,  # pandas
            #linestyle = '-',  # pandas
            #linewidth = None,  # pandas
            #cmap = None,  # pandas
            #color = None,  # pandas
            #grid = True,  # pandas
            #legend = True,
            yunit = None,  # test if {None, 'mmHg'}
            #title = None, pandas
            ax = None,  # test if else {None, ax}
            figsizepx = [641, 481],  # test non even pixel number
            ext = None,  # test if else {'jpeg', None}
            figbgcolor = None  # test if
        )
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(fig)) == "<class 'matplotlib.figure.Figure'>")
        plt.close()

    def test_mcdsts_plot_timeseries_none_none_none_conc_ax_jpeg(self, mcdsts=mcdsts):
        fig, ax = plt.subplots()
        s_pathfile = mcdsts.plot_timeseries(
            focus_cat = None,  # test if {None/total, 'voxel_i'}
            focus_num = None,  # test if {None/count, 'oxygen'}
            #aggregate_num = np.mean,  # pandas
            frame = 'conc',  # test if else {'df_cell', 'df_conc'}
            z_slice = None,  # test timeseries
            #logy = False,  # pandas
            #ylim = None,  # pandas
            #secondary_y = None,  # pandas
            #subplots = False,  # pandas
            #sharex = False,  # pandas
            #sharey = False,  # pandas
            #linestyle = '-',  # pandas
            #linewidth = None,  # pandas
            #cmap = None,  # pandas
            #color = None,  # pandas
            #grid = True,  # pandas
            #legend = True,
            yunit = None,  # test if {None, 'mmHg'}
            #title = None, pandas
            ax = ax,  # test if else {None, ax}
            figsizepx = [641, 481],  # test non even pixel number
            ext = 'jpeg',  # test if else {'jpeg', None}
            figbgcolor = None  # test if
        )
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (s_pathfile.endswith('/pcdl/data_timeseries_2d/timeseries_conc_total_count.jpeg')) and \
              (os.path.exists(s_pathfile))
        os.remove(s_pathfile)
        plt.close()

    def test_mcdsts_plot_timeseries_cat_none_yunit_conc(self, mcdsts=mcdsts):
        fig = mcdsts.plot_timeseries(
            focus_cat = 'voxel_i',  # test if {None/total, 'voxel_i'}
            focus_num = None,  # test if {None/count, 'oxygen'}
            #aggregate_num = np.mean,  # pandas
            frame = 'df_conc',  # test if else {'df_cell', 'df_conc'}
            z_slice = -0.3,  # test if if
            #logy = False,  # pandas
            #ylim = None,  # pandas
            #secondary_y = None,  # pandas
            #subplots = False,  # pandas
            #sharex = False,  # pandas
            #sharey = False,  # pandas
            #linestyle = '-',  # pandas
            #linewidth = None,  # pandas
            #cmap = None,  # pandas
            #color = None,  # pandas
            #grid = True,  # pandas
            #legend = True,
            yunit = 'mmHg',  # test if {None, 'mmHg'}
            #title = None, pandas
            ax = None,  # test if else {None, ax}
            figsizepx = [641, 481],  # test non even pixel number
            ext = None,  # test if else {'jpeg', None}
            figbgcolor = None  # test if
        )
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(fig)) == "<class 'matplotlib.figure.Figure'>")
        plt.close()

    def test_mcdsts_plot_timeseries_none_num_yunit_conc(self, mcdsts=mcdsts):
        fig = mcdsts.plot_timeseries(
            focus_cat = None,  # test if {None/total, 'voxel_i'}
            focus_num = 'oxygen',  # test if {None/count, 'oxygen'}
            #aggregate_num = np.mean,  # pandas
            frame = 'df_conc',  # test if else {'df_cell', 'df_conc'}
            z_slice = -0.3,  # test if if
            #logy = False,  # pandas
            #ylim = None,  # pandas
            #secondary_y = None,  # pandas
            #subplots = False,  # pandas
            #sharex = False,  # pandas
            #sharey = False,  # pandas
            #linestyle = '-',  # pandas
            #linewidth = None,  # pandas
            #cmap = None,  # pandas
            #color = None,  # pandas
            #grid = True,  # pandas
            #legend = True,
            yunit = 'mmHg',  # test if {None, 'mmHg'}
            #title = None, pandas
            ax = None,  # test if else {None, ax}
            figsizepx = [641, 481],  # test non even pixel number
            ext = None,  # test if else {'jpeg', None}
            figbgcolor = None  # test if
        )
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(fig)) == "<class 'matplotlib.figure.Figure'>")
        plt.close()

    def test_mcdsts_plot_timeseries_cat_num_none_conc(self, mcdsts=mcdsts):
        fig = mcdsts.plot_timeseries(
            focus_cat = 'voxel_i',  # test if {None/total, 'voxel_i'}
            focus_num = 'oxygen',  # test if {None/count, 'oxygen'}
            #aggregate_num = np.mean,  # pandas
            frame = 'conc',  # test if else {'df_cell', 'df_conc'}
            z_slice = -0.3,  # test if if
            #logy = False,  # pandas
            #ylim = None,  # pandas
            #secondary_y = None,  # pandas
            #subplots = False,  # pandas
            #sharex = False,  # pandas
            #sharey = False,  # pandas
            #linestyle = '-',  # pandas
            #linewidth = None,  # pandas
            #cmap = None,  # pandas
            #color = None,  # pandas
            #grid = True,  # pandas
            #legend = True,
            yunit = None,  # test if {None, 'mmHg'}
            #title = None, pandas
            ax = None,  # test if else {None, ax}
            figsizepx = [641, 481],  # test non even pixel number
            ext = None,  # test if else {'jpeg', None}
            figbgcolor = None  # test if
        )
        assert(str(type(mcdsts)) == "<class 'pcdl.pyMCDSts.pyMCDSts'>") and \
              (str(type(fig)) == "<class 'matplotlib.figure.Figure'>")
        plt.close()

