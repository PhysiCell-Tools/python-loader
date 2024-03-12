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
import os
import pathlib
import pcdl


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
    mcds = pcdl.pyMCDS(xmlfile=s_pathfile_3d)  # custom_type={}, microenv=True, graph=True, settingxml='PhysiCell_settings.xml', verbose=True

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


##################
# test for speed #
##################

class TestPyMcds3dMicroenvWorkhorse(object):
    ''' tests on 3D data set, for speed, for pcdl.pyMCDS microenvironment related workhorse functions. '''
    mcds = pcdl.pyMCDS(xmlfile=s_pathfile_3d)  # custom_type={}, microenv=True, graph=True, settingxml='PhysiCell_settings.xml', verbose=True

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


class TestPyMcds3dCellWorkhorse(object):
    ''' tests on 3D data set, for speed, for pcdl.pyMCDS cell related workhorse functions. '''
    mcds = pcdl.pyMCDS(xmlfile=s_pathfile_3d)  # custom_type={}, microenv=True, graph=True, settingxml='PhysiCell_settings.xml', verbose=True

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

class TestPyMcds3dGraphWorkhorse(object):
    ''' tests on 3D data set, for speed, for pcdl.pyMCDS graph related workhorse functions. '''
    mcds = pcdl.pyMCDS(xmlfile=s_pathfile_3d)  # custom_type={}, microenv=True, graph=True, settingxml='PhysiCell_settings.xml', verbose=True


    ## graph related functions ##
    # attached graph gml files
    def test_mcds_make_graph_gml_attached_defaultattr(self, mcds=mcds):
        s_pathfile = mcds.make_graph_gml(graph_type='attached', edge_attr=True, node_attr=[])
        f = open(s_pathfile)
        s_file = f.read()
        f.close()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.endswith('pcdl/data_timeseries_3d/output00000024_attached.gml')) and \
              (os.path.exists(s_pathfile)) and \
              (s_file.find('Creator "pcdl_v') > -1) and \
              (s_file.find('graph [\n  id 1440\n  comment "time_min"\n  label "attached_graph"\n  directed 0\n') > -1) and \
              (s_file.find('node [\n    id') > -1) and \
              (s_file.find('edge [\n    source') == -1) and \
              (s_file.find('distance_microns') == -1)
        os.remove(s_pathfile)

    def test_mcds_make_graph_gml_attached_edgeattrfalse(self, mcds=mcds):
        s_pathfile = mcds.make_graph_gml(graph_type='attached', edge_attr=False, node_attr=[])
        f = open(s_pathfile)
        s_file = f.read()
        f.close()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.endswith('pcdl/data_timeseries_3d/output00000024_attached.gml')) and \
              (os.path.exists(s_pathfile)) and \
              (s_file.find('Creator "pcdl_v') > -1) and \
              (s_file.find('graph [\n  id 1440\n  comment "time_min"\n  label "attached_graph"\n  directed 0\n') > -1) and \
              (s_file.find('node [\n    id') > -1) and \
              (s_file.find('edge [\n    source') == -1) and \
              (s_file.find('distance_microns') == -1)
        os.remove(s_pathfile)

    def test_mcds_make_graph_gml_neighbor_nodeattrtrue(self, mcds=mcds):
        s_pathfile = mcds.make_graph_gml(graph_type='neighbor', edge_attr=True, node_attr=['dead','cell_count_voxel','cell_density_micron3','cell_type'])  # bool,int,float,str
        f = open(s_pathfile)
        s_file = f.read()
        f.close()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.endswith('pcdl/data_timeseries_3d/output00000024_neighbor.gml')) and \
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
        s_pathfile = mcds.make_graph_gml(graph_type='neighbor', edge_attr=True, node_attr=[])
        f = open(s_pathfile)
        s_file = f.read()
        f.close()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.endswith('pcdl/data_timeseries_3d/output00000024_neighbor.gml')) and \
              (os.path.exists(s_pathfile)) and \
              (s_file.find('Creator "pcdl_v') > -1) and \
              (s_file.find('graph [\n  id 1440\n  comment "time_min"\n  label "neighbor_graph"\n  directed 0\n') > -1) and \
              (s_file.find('node [\n    id') > -1) and \
              (s_file.find('edge [\n    source') > -1) and \
              (s_file.find('distance_microns') > -1)
        os.remove(s_pathfile)

    def test_mcds_make_graph_gml_neighbor_edgeattrfalse(self, mcds=mcds):
        s_pathfile = mcds.make_graph_gml(graph_type='neighbor', edge_attr=False, node_attr=[])
        f = open(s_pathfile)
        s_file = f.read()
        f.close()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.endswith('pcdl/data_timeseries_3d/output00000024_neighbor.gml')) and \
              (os.path.exists(s_pathfile)) and \
              (s_file.find('Creator "pcdl_v') > -1) and \
              (s_file.find('graph [\n  id 1440\n  comment "time_min"\n  label "neighbor_graph"\n  directed 0\n') > -1) and \
              (s_file.find('node [\n    id') > -1) and \
              (s_file.find('edge [\n    source') > -1) and \
              (s_file.find('distance_microns') == -1)
        os.remove(s_pathfile)

    def test_mcds_make_graph_gml_neighbor_nodeattrtrue(self, mcds=mcds):
        s_pathfile = mcds.make_graph_gml(graph_type='neighbor', edge_attr=True, node_attr=['dead','cell_count_voxel','cell_density_micron3','cell_type'])  # bool,int,float,str
        f = open(s_pathfile)
        s_file = f.read()
        f.close()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (s_pathfile.endswith('pcdl/data_timeseries_3d/output00000024_neighbor.gml')) and \
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


## setting related functions ##
# BUE: TO BE COPIED FROM 2D
class TestPyMcds3dUnitWorkhorse(object):
    ''' tests on 3D data set, for speed, for pcdl.pyMCDS unit related workhorse functions. '''
    mcds = pcdl.pyMCDS(xmlfile=s_pathfile_3d)  # custom_type={}, microenv=True, graph=True, settingxml='PhysiCell_settings.xml', verbose=True

    def test_mcds_get_parameter_dict(self, mcds=mcds):
        d_parameter = mcds.get_parameter_dict()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(d_parameter)) == "<class 'dict'>") and \
              (len(d_parameter) == 82) and \
              (d_parameter['oxygen_initial_condition'] == 38.0)

    def test_mcds_get_rule_df(self, mcds=mcds):
        df_rule = mcds.get_rule_df()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (df_rule is None)

    def test_mcds_get_unit_dict(self, mcds=mcds):
        ds_unit = mcds.get_unit_dict()
        assert(str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>") and \
              (str(type(ds_unit)) == "<class 'dict'>") and \
              (len(ds_unit) == 151) and \
              (ds_unit['oxygen'] == 'mmHg')
