#####
# title: test_snapshot_3d.py
#
# language: python3
# author: Elmar Bucher
# date: 2022-10-15
# license: BSD 3-Clause
#
# description:
#   pytest unit test library for the pcdl library pyMCDS class.
#   focus is only in 3d and speed.
#   + https://docs.pytest.org/
#
#   note:
#   assert actual == expected, message
#   == value equality
#   is reference equality
#   pytest.approx for real values
#####


# load library
import numpy as np
import os
import pathlib
import pcdl
import matplotlib.pyplot as plt


# const
s_path_3d = str(pathlib.Path(pcdl.__file__).parent.resolve()/'data_timeseries_3d')
s_file_3d = 'output00000024.xml'
s_pathfile_3d = f'{s_path_3d}/{s_file_3d}'


## download test dataset ##
if not os.path.exists(s_path_3d):
    pcdl.install_data()


###########
# 3D only #
###########

class TestPyMcds3dOnly(object):
    ''' test for 3D only conditions in pcdl.pyMCDS functions. '''
    mcds = pcdl.pyMCDS(xmlfile=s_pathfile_3d)  # custom_data_type={}, microenv=True, graph=True, physiboss=True, settingxml='PhysiCell_settings.xml', verbose=True

    ## mesh related functions ##
    # bue: check if in 2D
    def test_mcds_get_mesh_spacing(self, mcds=mcds):
        lr_spacing = mcds.get_mesh_spacing()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(lr_spacing)) == "<class 'list'>") and \
              (str(type(lr_spacing[0])) == "<class 'numpy.float64'>") and \
              (str(type(lr_spacing[1])) == "<class 'numpy.float64'>") and \
              (str(type(lr_spacing[-1])) == "<class 'numpy.float64'>") and \
              (lr_spacing == [30.0, 20.0, 10.0])


############################
# test workhorse for speed #
############################

class TestPyMcdsInit(object):
    ''' tests for loading a pcdl.pyMCDS data set. '''
    mcds = pcdl.pyMCDS(xmlfile=s_file_3d, output_path=s_path_3d, custom_data_type={}, microenv=True, graph=True, physiboss=True, settingxml='PhysiCell_settings.xml', verbose=True)
    df_cell = mcds.get_cell_df()
    def test_mcds_init_microenv(self, mcds=mcds, df_cell=df_cell):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_cell.shape == (20460, 118))

    def test_mcds_init_graph(self, mcds=mcds):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(mcds.data['discrete_cells']['graph']['attached_cells'])) == "<class 'dict'>") and \
              (str(type(mcds.data['discrete_cells']['graph']['neighbor_cells'])) == "<class 'dict'>") and \
              (len(mcds.data['discrete_cells']['graph']['attached_cells']) == 20460) and \
              (len(mcds.data['discrete_cells']['graph']['neighbor_cells']) == 20460)

    def test_mcds_init_physiboss(self, mcds=mcds):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (mcds.data['discrete_cells']['physiboss'] == None)

    def test_mcds_init_settingxml(self, mcds=mcds, df_cell=df_cell):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (set(df_cell.columns).issuperset({'cancer_cell_attack_rates'})) and \
              (df_cell.shape == (20460, 118))


class TestPyMcdsInitMicroenvFalse(object):
    ''' tests for loading a pcdl.pyMCDS data set with microenv false. '''
    mcds = pcdl.pyMCDS(xmlfile=s_file_3d, output_path=s_path_3d, custom_data_type={}, microenv=False, graph=True, physiboss=True, settingxml='PhysiCell_settings.xml', verbose=True)
    df_cell = mcds.get_cell_df()
    def test_mcds_init_microenv(self, mcds=mcds, df_cell=df_cell):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_cell.shape == (20460, 112))

    def test_mcds_init_graph(self, mcds=mcds):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(mcds.data['discrete_cells']['graph']['attached_cells'])) == "<class 'dict'>") and \
              (str(type(mcds.data['discrete_cells']['graph']['neighbor_cells'])) == "<class 'dict'>") and \
              (len(mcds.data['discrete_cells']['graph']['attached_cells']) == 20460) and \
              (len(mcds.data['discrete_cells']['graph']['neighbor_cells']) == 20460)

    def test_mcds_init_physiboss(self, mcds=mcds):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (mcds.data['discrete_cells']['physiboss'] == None)

    def test_mcds_init_settingxml(self, mcds=mcds, df_cell=df_cell):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (set(df_cell.columns).issuperset({'cancer_cell_attack_rates'})) and \
              (df_cell.shape == (20460, 112))


class TestPyMcdsInitGraphFalse(object):
    ''' tests for loading a pcdl.pyMCDS data set with graph false. '''
    mcds = pcdl.pyMCDS(xmlfile=s_file_3d, output_path=s_path_3d, custom_data_type={}, microenv=True, graph=False, physiboss=True, settingxml='PhysiCell_settings.xml', verbose=True)
    df_cell = mcds.get_cell_df()
    def test_mcds_init_microenv(self, mcds=mcds, df_cell=df_cell):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_cell.shape == (20460, 118))

    def test_mcds_init_graph(self, mcds=mcds):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(mcds.data['discrete_cells']['graph']['attached_cells'])) == "<class 'dict'>") and \
              (str(type(mcds.data['discrete_cells']['graph']['neighbor_cells'])) == "<class 'dict'>") and \
              (len(mcds.data['discrete_cells']['graph']['attached_cells']) == 0) and \
              (len(mcds.data['discrete_cells']['graph']['neighbor_cells']) == 0)

    def test_mcds_init_physiboss(self, mcds=mcds):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (mcds.data['discrete_cells']['physiboss'] == None)

    def test_mcds_init_settingxml(self, mcds=mcds, df_cell=df_cell):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (set(df_cell.columns).issuperset({'cancer_cell_attack_rates'})) and \
              (df_cell.shape == (20460, 118))


class TestPyMcdsInitPhysibossFalse(object):
    ''' tests for loading a pcdl.pyMCDS data set with physiboss false. '''
    mcds = pcdl.pyMCDS(xmlfile=s_file_3d, output_path=s_path_3d, custom_data_type={}, microenv=True, graph=True, physiboss=False, settingxml='PhysiCell_settings.xml', verbose=True)
    df_cell = mcds.get_cell_df()
    def test_mcds_init_microenv(self, mcds=mcds, df_cell=df_cell):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_cell.shape == (20460, 118))

    def test_mcds_init_graph(self, mcds=mcds):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(mcds.data['discrete_cells']['graph']['attached_cells'])) == "<class 'dict'>") and \
              (str(type(mcds.data['discrete_cells']['graph']['neighbor_cells'])) == "<class 'dict'>") and \
              (len(mcds.data['discrete_cells']['graph']['attached_cells']) == 20460) and \
              (len(mcds.data['discrete_cells']['graph']['neighbor_cells']) == 20460)

    def test_mcds_init_physiboss(self, mcds=mcds):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (mcds.data['discrete_cells']['physiboss'] == None)

    def test_mcds_init_settingxml(self, mcds=mcds, df_cell=df_cell):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (set(df_cell.columns).issuperset({'cancer_cell_attack_rates'})) and \
              (df_cell.shape == (20460, 118))


#class TestPyMcdsInitSettingxmlFalse(object):
#    ''' tests for loading a pcdl.pyMCDS data set with settingxml false. '''
#    mcds = pcdl.pyMCDS(xmlfile=s_file_3d, output_path=s_path_3d, custom_data_type={}, microenv=True, graph=True, physiboss=True, settingxml=False, verbose=True)
#    df_cell = mcds.get_cell_df()
#    def test_mcds_init_microenv(self, mcds=mcds, df_cell=df_cell):
#        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
#              (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
#              (df_cell.shape == (20460, 118))
#
#    def test_mcds_init_graph(self, mcds=mcds):
#        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
#              (str(type(mcds.data['discrete_cells']['graph']['attached_cells'])) == "<class 'dict'>") and \
#              (str(type(mcds.data['discrete_cells']['graph']['neighbor_cells'])) == "<class 'dict'>") and \
#              (len(mcds.data['discrete_cells']['graph']['attached_cells']) == 20460) and \
#              (len(mcds.data['discrete_cells']['graph']['neighbor_cells']) == 20460)
#
#    def test_mcds_init_physiboss(self, mcds=mcds):
#        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
#              (mcds.data['discrete_cells']['physiboss'] == None)
#
#    def test_mcds_init_settingxml(self, mcds=mcds, df_cell=df_cell):
#        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
#              (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
#              (set(df_cell.columns).issuperset({'attack_rates_0'})) and \
#              (df_cell.shape == (20460, 118))


#class TestPyMcdsInitSettingxmlNone(object):
#    ''' tests for loading a pcdl.pyMCDS data set with settingxml none. '''
#    mcds = pcdl.pyMCDS(xmlfile=s_file_3d, output_path=s_path_3d, custom_data_type={}, microenv=True, graph=True, physiboss=True, settingxml=None, verbose=True)
#    df_cell = mcds.get_cell_df()
#    def test_mcds_init_microenv(self, mcds=mcds, df_cell=df_cell):
#        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
#              (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
#              (df_cell.shape == (20460, 118))
#
#    def test_mcds_init_graph(self, mcds=mcds):
#        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
#              (str(type(mcds.data['discrete_cells']['graph']['attached_cells'])) == "<class 'dict'>") and \
#              (str(type(mcds.data['discrete_cells']['graph']['neighbor_cells'])) == "<class 'dict'>") and \
#              (len(mcds.data['discrete_cells']['graph']['attached_cells']) == 20460) and \
#              (len(mcds.data['discrete_cells']['graph']['neighbor_cells']) == 20460)
#
#    def test_mcds_init_physiboss(self, mcds=mcds):
#        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
#              (mcds.data['discrete_cells']['physiboss'] == None)
#
#    def test_mcds_init_settingxml(self, mcds=mcds, df_cell=df_cell):
#        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
#              (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
#              (set(df_cell.columns).issuperset({'attack_rates_0'})) and \
#              (df_cell.shape == (20460, 118))


class TestPyMcdsInitVerboseTrue(object):
    ''' tests for loading a pcdl.pyMCDS data set and set_verbose_false function. '''
    mcds = pcdl.pyMCDS(xmlfile=s_file_3d, output_path=s_path_3d, custom_data_type={}, microenv=True, graph=True, physiboss=True, settingxml='PhysiCell_settings.xml', verbose=True)

    def test_mcds_verbose_true(self, mcds=mcds):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (mcds.verbose)

    def test_mcds_set_verbose_false(self, mcds=mcds):
        mcds.set_verbose_false()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (not mcds.verbose)


class TestPyMcdsInitVerboseFalse(object):
    ''' tests for loading a pcdl.pyMCDS data set and set_verbose_true function. '''
    mcds = pcdl.pyMCDS(xmlfile=s_file_3d, output_path=s_path_3d, custom_data_type={}, microenv=True, graph=True, physiboss=True, settingxml='PhysiCell_settings.xml', verbose=False)

    def test_mcds_verbose_false(self, mcds=mcds):
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (not mcds.verbose)

    def test_mcds_set_verbose_true(self, mcds=mcds):
        mcds.set_verbose_true()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (mcds.verbose)


class TestPyMcds3dSettingWorkhorse(object):
    ''' tests on 3D data set, for speed, for pcdl.pyMCDS unit related workhorse functions. '''
    mcds = pcdl.pyMCDS(xmlfile=s_pathfile_3d)  # custom_data_type={}, microenv=True, graph=True, physiboss=True, settingxml='PhysiCell_settings.xml', verbose=True

    def test_mcds_get_unit_dict(self, mcds=mcds):
        ds_unit = mcds.get_unit_dict()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(ds_unit)) == "<class 'dict'>") and \
              (len(ds_unit) == 105) and \
              (ds_unit['oxygen'] == 'mmHg')


class TestPyMcds3dMicroenvWorkhorse(object):
    ''' tests on 3D data set, for speed, for pcdl.pyMCDS microenvironment related workhorse functions. '''
    mcds = pcdl.pyMCDS(xmlfile=s_pathfile_3d)  # custom_data_type={}, microenv=True, graph=True, physiboss=True, settingxml='PhysiCell_settings.xml', verbose=True

    ## micro environment related functions ##
    def test_mcds_get_conc_df(self, mcds=mcds):
        df_conc = mcds.get_conc_df(z_slice=None, halt=False, values=1, drop=set(), keep=set())
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_conc)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_conc.shape == (1331, 11))

    def test_mcds_get_conc_df_zslice_center(self, mcds=mcds):
        df_conc = mcds.get_conc_df(z_slice=0, halt=False, values=1, drop=set(), keep=set())
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_conc)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_conc.shape == (121, 11))

    def test_mcds_get_conc_df_zslice_outofcenter(self, mcds=mcds):
        df_conc = mcds.get_conc_df(z_slice=-6, halt=False, values=1, drop=set(), keep=set())
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_conc)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_conc.shape == (121, 11))

    def test_mcds_get_conc_df_values(self, mcds=mcds):
        df_conc = mcds.get_conc_df(z_slice=None, halt=False, values=2, drop=set(), keep=set())
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_conc)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_conc.shape == (1331, 11))

    def test_mcds_get_conc_df_drop(self, mcds=mcds):
        df_conc = mcds.get_conc_df(z_slice=None, halt=False, values=1, drop={'oxygen'}, keep=set())
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_conc)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_conc.shape == (1331, 10))

    def test_mcds_get_conc_df_keep(self, mcds=mcds):
        df_conc = mcds.get_conc_df(z_slice=None, halt=False, values=1, drop=set(), keep={'oxygen'})
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_conc)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_conc.shape == (1331, 10))

    def test_mcds_plot_contour(self, mcds=mcds):
        fig, ax = plt.subplots()
        fig = mcds.plot_contour(
            'oxygen',
            z_slice = -3.333,  # test if
            vmin = None,  # test if
            vmax = None,  # test if
            #alpha = 1,  # matplotlib
            fill = False,  # contour case
            #cmap = 'viridis',  matplotlib
            title = 'test_mcds_plot_contour',  # test if
            #grid = False, # matplotlib
            xlim = [-31, 301],  # test if
            ylim = [-21, 201],  # test if
            xyequal = True, # test if
            ax = ax,  # use axis from existing matplotlib figure
            figsizepx = [641, 481],  # test if
            ext = None, # test fig case
            figbgcolor = None,  # not at file
        )
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(fig)) == "<class 'matplotlib.figure.Figure'>")
        plt.close()

    def test_mcds_plot_contourf(self, mcds=mcds):
        s_pathfile = mcds.plot_contour(
            'oxygen',
            z_slice = 0,  # jum over if
            vmin = None,  # test if
            vmax = None,  # test if
            #alpha = 1,  # matplotlib
            fill = True,  # contourf case
            #cmap = 'viridis',  # matplotlib
            title = 'test_mcds_plot_contourf',  # test if
            #grid = True,  # matplotlib
            xlim = None,  # jump over if
            ylim = None,  # jump over if
            xyequal = True,  # test if
            ax = None,  # generate fig ax case
            figsizepx = None,  # test if
            ext = 'tiff', # test file case
            figbgcolor = 'orange',  # jump over if
        )
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.replace('\\','/').endswith('/pcdl/data_timeseries_3d/conc_oxygen_z-5.0/output00000024_oxygen.tiff')) and \
              (os.path.exists(s_pathfile)) and \
              (os.path.getsize(s_pathfile) > 2**10)
        os.remove(s_pathfile)

    def test_mcds_make_conc_vtk(self, mcds=mcds):
        s_pathfile = mcds.make_conc_vtk()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.replace('\\','/').endswith('/pcdl/data_timeseries_3d/output00000024_conc.vtr')) and \
              (os.path.exists(s_pathfile)) and \
              (os.path.getsize(s_pathfile) > 2**10)
        os.remove(s_pathfile)


class TestPyMcds3dCellWorkhorse(object):
    ''' tests on 3D data set, for speed, for pcdl.pyMCDS cell related workhorse functions. '''
    mcds = pcdl.pyMCDS(xmlfile=s_pathfile_3d)  # custom_data_type={}, microenv=True, graph=True, physiboss=True, settingxml='PhysiCell_settings.xml', verbose=True

    ## cell related functions ##
    def test_mcds_get_cell_df(self, mcds=mcds):
        df_cell = mcds.get_cell_df(values=1, drop=set(), keep=set())
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_cell.shape == (20460, 118))

    def test_mcds_get_cell_df_values(self, mcds=mcds):
        df_cell = mcds.get_cell_df(values=2, drop=set(), keep=set())
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_cell.shape == (20460, 33))

    def test_mcds_get_cell_df_drop(self, mcds=mcds):
        df_cell = mcds.get_cell_df(values=1, drop={'oxygen'}, keep=set())
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_cell.shape == (20460, 117))

    def test_mcds_get_cell_df_keep(self, mcds=mcds):
        df_cell = mcds.get_cell_df(values=1, drop=set(), keep={'oxygen'})
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_cell.shape == (20460, 13))

    # scatter categorical
    def test_mcds_plot_scatter_cat_if(self, mcds=mcds):
        fig, ax = plt.subplots()
        fig = mcds.plot_scatter(
            focus='cell_type',  # case categorical
            z_slice = -3.333,   # test if
            z_axis = None,  # test if case categorical
            #alpha = 1,  # matplotlib
            cmap = 'viridis',  # else case es_categorical
            title = 'test_mcds_plot_scatter_cat', # matplotlib
            #grid = True,  # matplotlib
            #legend_loc = 'lower left',  # matplotlib
            xlim = None,  # test if
            ylim = None,  # test if
            xyequal = True,  # test if
            #s = None,  # matplotlib
            ax = ax,  # use axis from existing matplotlib figure
            figsizepx = [701, 501],  # jump over if case ax none
            ext = None,  # test fig case
            figbgcolor = None,  # not a file
        )
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(fig)) == "<class 'matplotlib.figure.Figure'>")
        plt.close()

    def test_mcds_plot_scatter_cat_else1(self, mcds=mcds):
        s_pathfile = mcds.plot_scatter(
            focus='cell_type',  # case categorical
            z_slice = 0,  # jump over if
            z_axis = {'cancer_cell'},  # test else case categorical
            #alpha = 1,  # matplotlib
            cmap = {'cancer_cell': 'maroon'},  # test if case es_categorical
            title ='test_mcds_plot_scatter_else',  # matplotlib
            #grid = True,  # matplotlib
            #legend_loc = 'lower left',  # matplotlib
            xlim = [-31, 301],  # jump over if
            ylim = [-21, 201],  # jump over if
            xyequal = False,  # jump over if
            #s = None,  # matplotlib
            ax = None,  # generate matplotlib figure
            figsizepx = None,  # test if case ax none
            ext = 'tiff',  # test file case
            figbgcolor = 'lime',  # jump over if
        )
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.replace('\\','/').endswith('/pcdl/data_timeseries_3d/cell_cell_type_z-5.0/output00000024_cell_type.tiff')) and \
              (os.path.exists(s_pathfile)) and \
              (os.path.getsize(s_pathfile) > 2**10)
        os.remove(s_pathfile)

    def test_mcds_plot_scatter_cat_else2(self, mcds=mcds):
        fig, ax = plt.subplots()
        mcds.plot_scatter(
            focus='cell_type',  # case categorical
            z_slice = 0,  # jump over if
            z_axis = {'cancer_cell'},  # test else case categorical
            #alpha = 1,  # matplotlib
            cmap = 'viridis',  # test else case es_categorical
            title ='test_mcds_plot_scatter_else2',  # matplotlib
            #grid = True,  # matplotlib
            #legend_loc = 'lower left',  # matplotlib
            xlim = None,  # test if
            ylim = None,  # test if
            xyequal = True,  # test if
            #s = None,  # matplotlib
            ax = ax,  # use axis from existing matplotlib figure
            #figsizepx = None,  # test case ax ax
            ext = None,  # test fig case
            figbgcolor = None,  # not a file
        )
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(fig)) == "<class 'matplotlib.figure.Figure'>")
        plt.close()

    # scatter numerical
    def test_mcds_plot_scatter_num_if(self, mcds=mcds):
        fig = mcds.plot_scatter(
            focus='oxygen',  # case numeric
            z_slice = -3.333,   # test if
            z_axis = None,  # test if numeric
            #alpha = 1,  # matplotlib
            #cmap = 'viridis',  # matplotlib
            #title = None, # matplotlib
            #grid = True,  # matplotlib
            #legend_loc = 'lower left',  # matplotlib
            xlim = None,  # test if
            ylim = None,  # test if
            #xyequal = True,  # test if
            #s = None,  # matplotlib
            #ax = None,  # generate matplotlib figure
            #figsizepx = None,  # test if case
            ext = None,  # test fig case
            figbgcolor = None,  # not a file
        )
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(fig)) == "<class 'matplotlib.figure.Figure'>")
        plt.close()

    def test_mcds_plot_scatter_num_else(self, mcds=mcds):
        fig = mcds.plot_scatter(
            focus='oxygen',  # case numeric
            z_slice = 0,   # jump over if
            z_axis = [0, 38],  # test else numeric
            #alpha = 1,  # matplotlib
            #cmap = 'viridis',  # matplotlib
            #title = None, # matplotlib
            #grid = True,  # matplotlib
            #legend_loc = 'lower left',  # matplotlib
            xlim = None,  # test if
            ylim = None,  # test if
            #xyequal = True,  # test if
            #s = None,  # matplotlib
            #ax = None,  # generate matplotlib figure
            #figsizepx = None,  # if case
            ext = None,  # test fig case
            figbgcolor = None,  # not a file
        )
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(fig)) == "<class 'matplotlib.figure.Figure'>")
        plt.close()

    def test_mcds_make_cell_vtk_attribute_default(self, mcds=mcds):
        s_pathfile = mcds.make_cell_vtk(
            #attribute=['cell_type'],
            visualize=False,
        )
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.replace('\\','/').endswith('/pcdl/data_timeseries_3d/output00000024_cell.vtp')) and \
              (os.path.exists(s_pathfile)) and \
              (os.path.getsize(s_pathfile) > 2**10)
        os.remove(s_pathfile)

    def test_mcds_make_cell_vtk_attribute_zero(self, mcds=mcds):
        s_pathfile = mcds.make_cell_vtk(attribute=[], visualize=False)
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.replace('\\','/').endswith('/pcdl/data_timeseries_3d/output00000024_cell.vtp')) and \
              (os.path.exists(s_pathfile)) and \
              (os.path.getsize(s_pathfile) > 2**10)
        os.remove(s_pathfile)

    def test_mcds_make_cell_vtk_attribute_many(self, mcds=mcds):
        s_pathfile = mcds.make_cell_vtk(
            attribute=['dead', 'cell_count_voxel', 'pressure', 'cell_type'],
            visualize=False,
        )
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.replace('\\','/').endswith('/pcdl/data_timeseries_3d/output00000024_cell.vtp')) and \
              (os.path.exists(s_pathfile)) and \
              (os.path.getsize(s_pathfile) > 2**10)
        os.remove(s_pathfile)


class TestPyMcds3dGraphWorkhorse(object):
    ''' tests on 3D data set, for speed, for pcdl.pyMCDS graph related workhorse functions. '''
    mcds = pcdl.pyMCDS(xmlfile=s_pathfile_3d)  # custom_data_type={}, microenv=True, graph=True, physiboss=True, settingxml='PhysiCell_settings.xml', verbose=True

    ## graph related functions ##
    # attached graph gml files
    def test_mcds_make_graph_gml_attached_defaultattr(self, mcds=mcds):
        s_pathfile = mcds.make_graph_gml(graph_type='attached', edge_attribute=True, node_attribute=[])
        f = open(s_pathfile)
        s_file = f.read()
        f.close()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.replace('\\','/').endswith('pcdl/data_timeseries_3d/output00000024_attached.gml')) and \
              (os.path.exists(s_pathfile)) and \
              (s_file.find('Creator "pcdl_v') > -1) and \
              (s_file.find('graph [\n  id 1440\n  comment "time_min"\n  label "attached_graph"\n  directed 0\n') > -1) and \
              (s_file.find('node [\n    id') > -1) and \
              (s_file.find('edge [\n    source') == -1) and \
              (s_file.find('distance_microns') == -1)
        os.remove(s_pathfile)

    def test_mcds_make_graph_gml_attached_edgeattrfalse(self, mcds=mcds):
        s_pathfile = mcds.make_graph_gml(graph_type='attached', edge_attribute=False, node_attribute=[])
        f = open(s_pathfile)
        s_file = f.read()
        f.close()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.replace('\\','/').endswith('pcdl/data_timeseries_3d/output00000024_attached.gml')) and \
              (os.path.exists(s_pathfile)) and \
              (s_file.find('Creator "pcdl_v') > -1) and \
              (s_file.find('graph [\n  id 1440\n  comment "time_min"\n  label "attached_graph"\n  directed 0\n') > -1) and \
              (s_file.find('node [\n    id') > -1) and \
              (s_file.find('edge [\n    source') == -1) and \
              (s_file.find('distance_microns') == -1)
        os.remove(s_pathfile)

    def test_mcds_make_graph_gml_neighbor_nodeattrtrue(self, mcds=mcds):
        s_pathfile = mcds.make_graph_gml(graph_type='neighbor', edge_attribute=True, node_attribute=['dead','cell_count_voxel','cell_density_micron3','cell_type'])  # bool,int,float,str
        f = open(s_pathfile)
        s_file = f.read()
        f.close()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.replace('\\','/').endswith('pcdl/data_timeseries_3d/output00000024_neighbor.gml')) and \
              (os.path.exists(s_pathfile)) and \
              (s_file.find('Creator "pcdl_v') > -1) and \
              (s_file.find('graph [\n  id 1440\n  comment "time_min"\n  label "neighbor_graph"\n  directed 0\n') > -1) and \
              (s_file.find('node [\n    id') > -1) and \
              (s_file.find('dead') == -1) and \
              (s_file.find('cell_count_voxel') == -1) and \
              (s_file.find('cell_density_micron3') == -1) and \
              (s_file.find('cell_type') == -1) and \
              (s_file.find('edge [\n    source') > -1) and \
              (s_file.find('distance_microns')> -1)
        os.remove(s_pathfile)

    # neighbor graph gml file
    def test_mcds_make_graph_gml_neighbor_defaultattr(self, mcds=mcds):
        s_pathfile = mcds.make_graph_gml(graph_type='neighbor', edge_attribute=True, node_attribute=[])
        f = open(s_pathfile)
        s_file = f.read()
        f.close()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.replace('\\','/').endswith('pcdl/data_timeseries_3d/output00000024_neighbor.gml')) and \
              (os.path.exists(s_pathfile)) and \
              (s_file.find('Creator "pcdl_v') > -1) and \
              (s_file.find('graph [\n  id 1440\n  comment "time_min"\n  label "neighbor_graph"\n  directed 0\n') > -1) and \
              (s_file.find('node [\n    id') > -1) and \
              (s_file.find('edge [\n    source') > -1) and \
              (s_file.find('distance_microns') > -1)
        os.remove(s_pathfile)

    def test_mcds_make_graph_gml_neighbor_edgeattrfalse(self, mcds=mcds):
        s_pathfile = mcds.make_graph_gml(graph_type='neighbor', edge_attribute=False, node_attribute=[])
        f = open(s_pathfile)
        s_file = f.read()
        f.close()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.replace('\\','/').endswith('pcdl/data_timeseries_3d/output00000024_neighbor.gml')) and \
              (os.path.exists(s_pathfile)) and \
              (s_file.find('Creator "pcdl_v') > -1) and \
              (s_file.find('graph [\n  id 1440\n  comment "time_min"\n  label "neighbor_graph"\n  directed 0\n') > -1) and \
              (s_file.find('node [\n    id') > -1) and \
              (s_file.find('edge [\n    source') > -1) and \
              (s_file.find('distance_microns') == -1)
        os.remove(s_pathfile)

    def test_mcds_make_graph_gml_neighbor_nodeattrtrue(self, mcds=mcds):
        s_pathfile = mcds.make_graph_gml(graph_type='neighbor', edge_attribute=True, node_attribute=['dead','cell_count_voxel','cell_density_micron3','cell_type'])  # bool,int,float,str
        f = open(s_pathfile)
        s_file = f.read()
        f.close()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.replace('\\','/').endswith('pcdl/data_timeseries_3d/output00000024_neighbor.gml')) and \
              (os.path.exists(s_pathfile)) and \
              (s_file.find('Creator "pcdl_v') > -1) and \
              (s_file.find('graph [\n  id 1440\n  comment "time_min"\n  label "neighbor_graph"\n  directed 0\n') > -1) and \
              (s_file.find('node [\n    id') > -1) and \
              (s_file.find('dead') > -1) and \
              (s_file.find('cell_count_voxel') > -1) and \
              (s_file.find('cell_density_micron3') > -1) and \
              (s_file.find('cell_type') > -1) and \
              (s_file.find('edge [\n    source') > -1) and \
              (s_file.find('distance_microns') > -1)
        os.remove(s_pathfile)


class TestPyMcds3dOmeTiffWorkhorse(object):
    ''' tests on 3D data set, for speed, for pcdl.pyMCDS ome tiff related workhorse functions. '''
    mcds = pcdl.pyMCDS(xmlfile=s_pathfile_3d)  # custom_data_type={}, microenv=True, graph=True, physiboss=True, settingxml='PhysiCell_settings.xml', verbose=True

    ## ome tiff related functions ##
    def test_mcds_make_ome_tiff_default(self, mcds=mcds):
        s_pathfile = mcds.make_ome_tiff(cell_attribute='ID', conc_cutoff={}, focus=None, file=True)
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.replace('\\','/').endswith('pcdl/data_timeseries_3d/output00000024_oxygen_immunostimulatory_factor_cancer_cell_immune_cell_ID.ome.tiff')) and \
              (os.path.exists(s_pathfile)) and \
              (os.path.getsize(s_pathfile) > 2**10)
        os.remove(s_pathfile)

    def test_mcds_make_ome_tiff_nofile(self, mcds=mcds):
        a_ometiff = mcds.make_ome_tiff(cell_attribute='ID', conc_cutoff={}, focus=None, file=False)
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(a_ometiff)) == "<class 'numpy.ndarray'>") and \
              (a_ometiff.dtype == np.float32) and \
              (a_ometiff.shape == (4, 11, 200, 300))

