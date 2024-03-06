####
# title: test_timeseries.py
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
import os
import pathlib
import pcdl
import platform
import shutil


# const
s_path_2d = str(pathlib.Path(pcdl.__file__).parent.resolve()/'data_timeseries_2d')


## download test data ##
if not os.path.exists(s_path_2d):
    pcdl.install_data()


## making movies related functions ##

class TestPyMcdsTsMovies(object):
    ''' tests for loading a pcdl.pyMCDS data set. '''
    mcdsts = pcdl.pyMCDSts(s_path_2d, verbose=False)

    ## make_gif and magick ommand ##
    def test_mcdsts_make_gif_jpeg(self, mcdsts=mcdsts):
        s_opath = mcdsts.plot_scatter()
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
        s_opath = mcdsts.plot_scatter(ext='tiff')
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
        s_opath = mcdsts.plot_scatter()
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
        s_opath = mcdsts.plot_scatter(ext='tiff')
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
        s_opath = mcdsts.plot_scatter()
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

"""
## data loading related functions ##

class TestPyMcdsTsInit(object):
    ''' tests for loading a pcdl.pyMCDSts data set. '''
    mcdsts = pcdl.pyMCDSts(s_path_2d, verbose=False)

    def test_mcdsts_set_verbose_true(self, mcdsts=mcdsts):
        mcdsts.set_verbose_true()
        assert mcdsts.verbose

    def test_mcdsts_set_verbose_false(self, mcdsts=mcdsts):
        mcdsts.set_verbose_false()
        assert not mcdsts.verbose

    ## get_xmlfile and read_mcds command and get_mcds_list ##
    def test_mcdsts_get_xmlfile_list(self, mcdsts=mcdsts):
        ls_xmlfile = mcdsts.get_xmlfile_list()
        assert len(ls_xmlfile) == 25

    def test_mcdsts_get_mcds_list(self, mcdsts=mcdsts):
        l_mcds = mcdsts.get_mcds_list()
        assert l_mcds == mcdsts.l_mcds

    def test_mcdsts_read_mcds(self, mcdsts=mcdsts):
        l_mcds = mcdsts.read_mcds()
        assert len(l_mcds) == 25 and \
               len(mcdsts.l_mcds) == 25 and \
               mcdsts.l_mcds[-1].get_time() == 1440

    def test_mcdsts_get_xmlfile_list_read_mcds(self):
        mcdsts = pcdl.pyMCDSts(s_path_2d, load=False, verbose=False)
        ls_xmlfile = mcdsts.get_xmlfile_list()
        ls_xmlfile = ls_xmlfile[-3:]
        l_mcds = mcdsts.read_mcds(ls_xmlfile)
        assert len(ls_xmlfile) == 3 and \
               len(l_mcds) == 3 and \
               len(mcdsts.l_mcds) == 3 and \
               mcdsts.l_mcds[2].get_time() == 1440


## micro environment related functions ##

class TestPyMcdsTsMicroenv(object):
    ''' tests for pcdl.pyMCDS micro environment related functions. '''
    mcdsts = pcdl.pyMCDSts(s_path_2d, verbose=False)

    def test_mcdsts_get_conc_df(self, mcdsts=mcdsts):
        ldf_conc = mcdsts.get_conc_df(values=2, drop=set(), keep=set(), collapse=False)
        assert len(ldf_conc) == 25 and \
               ldf_conc[0].shape == (121, 9) and \
               ldf_conc[-1].shape == (121, 10)

    def test_mcdsts_get_conc_df_collapse(self, mcdsts=mcdsts):
        df_conc = mcdsts.get_conc_df(values=2, drop=set(), keep=set(), collapse=True)
        assert df_conc.shape == (3025, 10)

    def test_mcdsts_get_conc_df_features(self, mcdsts=mcdsts):
        dl_conc = mcdsts.get_conc_df_features(values=2, drop=set(), keep=set(), allvalues=False)
        assert len(dl_conc.keys()) == 1 and \
               len(dl_conc['oxygen']) == 2

    def test_mcdsts_get_conc_df_features_allvalues(self, mcdsts=mcdsts):
        dl_conc = mcdsts.get_conc_df_features(values=2, drop=set(), keep=set(), allvalues=True)
        assert len(dl_conc.keys()) == 1 and \
               len(dl_conc['oxygen']) > 2

    ## plot_contour command ##
    def test_mcdsts_plot_contour(self, mcdsts=mcdsts):
        s_path = mcdsts.plot_contour(
            focus = 'oxygen',
            z_slice = -3.333,  # test if
            extrema = None,  # test if
            #alpha = 1,  # matplotlib
            #fill = True,  # mcds.plot_contour
            #cmap = 'viridis',  # matplotlib
            #grid = True,  # matplotlib
            xlim = None,  # test if
            ylim = None,  # test if
            xyequal = True,  # test if
            figsizepx = [641, 481],  # test non even pixel number
            ext = 'jpeg',
            figbgcolor = None,  # test if
        )
        assert os.path.exists(s_path + 'output00000000_oxygen.jpeg') and \
               os.path.exists(s_path + 'output00000024_oxygen.jpeg')
        shutil.rmtree(s_path)


## cell related functions ##

class TestPyMcdsCell(object):
    ''' tests for pcdl.pyMCDS cell related functions. '''
    mcdsts = pcdl.pyMCDSts(s_path_2d, verbose=False)

    def test_mcdsts_get_cell_df(self, mcdsts=mcdsts):
        ldf_cell = mcdsts.get_cell_df(values=2, drop=set(), keep=set(), collapse=False)
        assert len(ldf_cell) == 25 and \
               ldf_cell[0].shape == (889, 19) and \
               ldf_cell[-1].shape == (1099, 40)

    def test_mcdsts_get_cell_df_collapse(self, mcdsts=mcdsts):
        df_cell = mcdsts.get_cell_df(values=2, drop=set(), keep=set(), collapse=True)
        assert df_cell.shape == (24758, 41)

    def test_mcdsts_get_cell_df_features(self, mcdsts=mcdsts):
        dl_cell = mcdsts.get_cell_df_features(values=2, drop=set(), keep=set(), allvalues=False)
        assert len(dl_cell.keys()) == 28 and \
               len(dl_cell['oxygen']) == 2

    def test_mcdsts_get_cell_df_features_allvalues(self, mcdsts=mcdsts):
        dl_cell = mcdsts.get_cell_df_features(values=2, drop=set(), keep=set(), allvalues=True)
        assert len(dl_cell.keys()) == 28 and \
               len(dl_cell['oxygen']) > 2

    ## plot_scatter command ##
    def test_mcdsts_plot_scatter_cat(self, mcdsts=mcdsts):
        s_path = mcdsts.plot_scatter(
            focus='cell_type',  # case categorical
            z_slice = -3.333,   # test if
            z_axis = None,  # test if categorical
            #alpha = 1,  # matplotlib
            #cmap = 'viridis',  # matplotlib
            #grid = True,  # matplotlib
            #legend_loc='lower left',  # matplotlib
            xlim = None,  # test if
            ylim = None,  # test if
            xyequal = True,  # test if
            s = None,  # test if
            figsizepx = [641, 481],  # case non even pixel number
            ext = 'jpeg', # test if
            figbgcolor = None,  # test if
        )
        assert os.path.exists(s_path + 'output00000000_cell_type.jpeg') and \
               os.path.exists(s_path + 'output00000000_cell_type.jpeg')
        shutil.rmtree(s_path)

    def test_mcdsts_plot_scatter_cat_cmap(self, mcdsts=mcdsts):
        s_path = mcdsts.plot_scatter(
            focus='cell_type',  # case categorical
            z_slice = -3.333,   # test if
            z_axis = None,  # test if categorical
            #alpha = 1,  # matplotlib
            cmap = {'cancer_cell': 'maroon'},  # matplotlib
            #grid = True,  # matplotlib
            #legend_loc='lower left',  # matplotlib
            xlim = None,  # test if
            ylim = None,  # test if
            xyequal = True,  # test if
            s = None,  # test if
            figsizepx = [641, 481],  # case non even pixel number
            ext = 'jpeg', # test if
            figbgcolor = None,  # test if
        )
        assert os.path.exists(s_path + 'output00000000_cell_type.jpeg') and \
               os.path.exists(s_path + 'output00000024_cell_type.jpeg')
        shutil.rmtree(s_path)

    def test_mcdsts_plot_scatter_num(self, mcdsts=mcdsts):
        s_path = mcdsts.plot_scatter(
            focus='pressure',  # case numeric
            z_slice = -3.333,   # test if
            z_axis = None,  # test if numeric
            #alpha = 1,  # matplotlib
            #cmap = 'viridis',  # matplotlib
            #grid = True,  # matplotlib
            #legend_loc='lower left',  # matplotlib
            xlim = None,  # test if
            ylim = None,  # test if
            xyequal = True,  # test if
            s = None,  # test if
            figsizepx = None,  # case extract from initial.svg
            ext = 'jpeg', # test if
            figbgcolor = None,  # test if
        )
        assert os.path.exists(s_path + 'output00000000_pressure.jpeg') and \
               os.path.exists(s_path + 'output00000024_pressure.jpeg')
        shutil.rmtree(s_path)


## graph related functions ##

class TestPyMcdsGraph(object):
    ''' tests for pcdl.pyMCDS graph related functions. '''
    mcdsts = pcdl.pyMCDSts(s_path_2d, verbose=False)

    ## graph related functions ##
    def test_mcdsts_get_graph_gml_attached_defaultattr(self, mcdsts=mcdsts):
        ls_pathfile = mcdsts.make_graph_gml(graph_type='attached', edge_attr=True, node_attr=['cell_type'])
        s_pathfile = mcdsts.output_path + 'output00000012_attached.gml'
        assert (len(ls_pathfile) == 25) and \
               (os.path.exists(s_pathfile))
        for s_pathfile in ls_pathfile:
            os.remove(s_pathfile)

    def test_mcdsts_get_graph_gml_neighbor_noneattr(self, mcdsts=mcdsts):
        ls_pathfile = mcdsts.make_graph_gml(graph_type='neighbor', edge_attr=False, node_attr=[])
        s_pathfile = mcdsts.output_path + 'output00000012_neighbor.gml'
        assert (len(ls_pathfile) == 25) and \
               (os.path.exists(s_pathfile))
        for s_pathfile in ls_pathfile:
            os.remove(s_pathfile)

    def test_mcdsts_get_graph_gml_neighbor_allattr(self, mcdsts=mcdsts):
        ls_pathfile = mcdsts.make_graph_gml(graph_type='neighbor', edge_attr=True, node_attr=['cell_type','dead','cell_count_voxel','cell_density_micron3'])
        s_pathfile = mcdsts.output_path + 'output00000012_neighbor.gml'
        assert (len(ls_pathfile) == 25) and \
               (os.path.exists(s_pathfile))
        for s_pathfile in ls_pathfile:
            os.remove(s_pathfile)


## Timeseries related functions ##

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
            frame = 'cell_df',  # test if else {'df_cell', 'df_conc'}
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
            ax = ax,  # test if else {None, ax}
            figsizepx = [641, 481],  # test non even pixel number
            ext = 'jpeg',  # test if else {'jpeg', None}
            figbgcolor = None  # test if
        )
        assert os.path.exists(s_pathfile)
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
        assert(str(type(fig)) == "<class 'matplotlib.figure.Figure'>")

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
        assert(str(type(fig)) == "<class 'matplotlib.figure.Figure'>")

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
        assert(str(type(fig)) == "<class 'matplotlib.figure.Figure'>")

    def test_mcdsts_plot_timeseries_none_none_none_conc_ax_jpeg(self, mcdsts=mcdsts):
        fig, ax = plt.subplots()
        s_pathfile = mcdsts.plot_timeseries(
            focus_cat = None,  # test if {None/total, 'voxel_i'}
            focus_num = None,  # test if {None/count, 'oxygen'}
            #aggregate_num = np.mean,  # pandas
            frame = 'conc_df',  # test if else {'df_cell', 'df_conc'}
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
            ax = ax,  # test if else {None, ax}
            figsizepx = [641, 481],  # test non even pixel number
            ext = 'jpeg',  # test if else {'jpeg', None}
            figbgcolor = None  # test if
        )
        assert os.path.exists(s_pathfile)
        os.remove(s_pathfile)

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
        assert(str(type(fig)) == "<class 'matplotlib.figure.Figure'>")

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
        assert(str(type(fig)) == "<class 'matplotlib.figure.Figure'>")

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
        assert(str(type(fig)) == "<class 'matplotlib.figure.Figure'>")
"""
