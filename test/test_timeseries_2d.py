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
import glob
import matplotlib.pyplot as plt
import os
import pathlib
import pcdl
import platform
import shutil


# const
s_path_2d = str(pathlib.Path(pcdl.__file__).parent.resolve()/'data_timeseries_2d')


# test data
if not os.path.exists(s_path_2d):
    pcdl.install_data()


# load physicell data time series
class TestPyMcdsTs(object):
    ''' test for pcdl.pyMCDSts data loader. '''
    mcdsts = pcdl.pyMCDSts(s_path_2d, verbose=False)

    ## get_xmlfile and read_mcds command and get_mcds_list ##
    def test_mcdsts_get_xmlfile_list(self, mcdsts=mcdsts):
        ls_xmlfile = mcdsts.get_xmlfile_list()
        assert len(ls_xmlfile) == 25

    def test_mcdsts_get_xmlfile_list_read_mcds(self, mcdsts=mcdsts):
        ls_xmlfile = mcdsts.get_xmlfile_list()
        ls_xmlfile = ls_xmlfile[-3:]
        l_mcds = mcdsts.read_mcds(ls_xmlfile)
        assert len(ls_xmlfile) == 3 and \
               len(l_mcds) == 3 and \
               len(mcdsts.l_mcds) == 3 and \
               mcdsts.l_mcds[2].get_time() == 1440

    def test_mcdsts_read_mcds(self, mcdsts=mcdsts):
        l_mcds = mcdsts.read_mcds()
        assert len(l_mcds) == 25 and \
               len(mcdsts.l_mcds) == 25 and \
               mcdsts.l_mcds[-1].get_time() == 1440

    def test_mcdsts_get_mcds_list(self, mcdsts=mcdsts):
        l_mcds = mcdsts.get_mcds_list()
        assert l_mcds == mcdsts.l_mcds

    ## data triage command ##
    def test_mcdsts_get_cell_df_features(self, mcdsts=mcdsts):
        dl_cell = mcdsts.get_cell_df_features(values=2, drop=set(), keep=set(), allvalues=False)
        assert len(dl_cell.keys()) == 28 and \
               len(dl_cell['oxygen']) == 2

    def test_mcdsts_get_cell_df_features_allvalues(self, mcdsts=mcdsts):
        dl_cell = mcdsts.get_cell_df_features(values=2, drop=set(), keep=set(), allvalues=True)
        assert len(dl_cell.keys()) == 28 and \
               len(dl_cell['oxygen']) > 2

    def test_mcdsts_get_conc_df_features(self, mcdsts=mcdsts):
        dl_conc = mcdsts.get_conc_df_features(values=2, drop=set(), keep=set(), allvalues=False)
        assert len(dl_conc.keys()) == 1 and \
               len(dl_conc['oxygen']) == 2

    def test_mcdsts_get_conc_df_features_allvalues(self, mcdsts=mcdsts):
        dl_conc = mcdsts.get_conc_df_features(values=2, drop=set(), keep=set(), allvalues=True)
        assert len(dl_conc.keys()) == 1 and \
               len(dl_conc['oxygen']) > 2

    ## magick command ##
    def test_mcdsts_handle_magick(self, mcdsts=mcdsts):
        s_magick = mcdsts._handle_magick()
        if not((os.system('magick --version') == 0) or ((platform.system() in {'Linux'}) and (os.system('convert --version') == 0))):
            s_magick = None
            print('Error @ pyMCDSts._handle_magick : image magick installation version >= 7.0 missing!')
        assert s_magick in {'', 'magick '}

    ## plot_scatter command ##
    def test_mcdsts_plot_scatter_cat(self, mcdsts=mcdsts):
        s_path = mcdsts.plot_scatter(
            focus='cell_type',  # case categorical
            z_slice = -3.333,   # test if
            z_axis = None,  # test if categorical
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
        assert os.path.exists(s_path + 'cell_type_000000000.0.jpeg') and \
               os.path.exists(s_path + 'cell_type_000001440.0.jpeg')
        shutil.rmtree(s_path)

    def test_mcdsts_plot_scatter_cat_cmap(self, mcdsts=mcdsts):
        s_path = mcdsts.plot_scatter(
            focus='cell_type',  # case categorical
            z_slice = -3.333,   # test if
            z_axis = None,  # test if categorical
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
        assert os.path.exists(s_path + 'cell_type_000000000.0.jpeg') and \
               os.path.exists(s_path + 'cell_type_000001440.0.jpeg')
        shutil.rmtree(s_path)

    def test_mcdsts_plot_scatter_num(self, mcdsts=mcdsts):
        s_path = mcdsts.plot_scatter(
            focus='pressure',  # case numeric
            z_slice = -3.333,   # test if
            z_axis = None,  # test if numeric
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
        assert os.path.exists(s_path + 'pressure_000000000.0.jpeg') and \
               os.path.exists(s_path + 'pressure_000001440.0.jpeg')
        shutil.rmtree(s_path)

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
        assert os.path.exists(s_path + 'oxygen_000000000.0.jpeg') and \
               os.path.exists(s_path + 'oxygen_000001440.0.jpeg')
        shutil.rmtree(s_path)

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

    ## make_gif command ##
    def test_mcdsts_make_gif(self, mcdsts=mcdsts):
        s_path = mcdsts.plot_scatter()
        s_opathfile = mcdsts.make_gif(
            path = s_path,
            #interface = 'jpeg',
        )
        assert os.path.exists(s_opathfile) and \
            (s_opathfile == s_path+'cell_cell_type_z0.0_jpeg.gif')
        shutil.rmtree(s_path)

    ## make_movie command ##
    def test_mcdsts_make_movie(self, mcdsts=mcdsts):
        s_path = mcdsts.plot_scatter()
        s_opathfile = mcdsts.make_movie(
            path = s_path,
            #interface = 'jpeg',
            #framerate = 12,
        )
        assert os.path.exists(s_opathfile) and \
            (s_opathfile == s_path+'cell_cell_type_z0.0_jpeg12.mp4')
        shutil.rmtree(s_path)

    ## graph related functions ##
    def test_mcds_get_graph_gml_attached(self, mcdsts=mcdsts):
        s_pathfile_glob = s_path_2d + '/graph_attached_*min.gml'
        s_pathfile = s_path_2d + '/graph_attached_00000720min.gml'
        mcdsts.make_graph_gml(graph_type='attached', node_attr=['cell_type'], edge_attr=True)
        assert os.path.exists(s_pathfile)
        for s_file in glob.glob(s_pathfile_glob):
            os.remove(s_file)

    def test_mcds_get_graph_gml_neighbor(self, mcdsts=mcdsts):
        s_pathfile_glob = s_path_2d + '/graph_neighbor_*min.gml'
        s_pathfile = s_path_2d + '/graph_neighbor_00000720min.gml'
        mcdsts.make_graph_gml(graph_type='neighbor', node_attr=['cell_type'], edge_attr=True)
        assert os.path.exists(s_pathfile)
        for s_file in glob.glob(s_pathfile_glob):
            os.remove(s_file)
