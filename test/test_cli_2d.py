####
# title: test_cli_2d.py
#
# language: python3
# author: Elmar Bucher
# date: 2023-02-25
# license: BSD 3-Clause
#
# description:
#   pytest unit test library for the pcdl library comand line interface functions.
#   + https://docs.pytest.org/
#
#   note:
#   assert actual == expected, message
#   == value equality
#   is reference equality
#   pytest.approx for real values
#####


# load library
import json
import os
import pandas as pd
import pathlib
import pcdl
import pytest
import shutil
import subprocess

# const
s_path_2d = str(pathlib.Path(pcdl.__file__).parent.resolve()/'data_timeseries_2d')
s_file_2d = 'output00000024.xml'
s_pathfile_2d = f'{s_path_2d}/{s_file_2d}'
print("s_path_2d", s_path_2d)
print("s_pathfile_2d", s_pathfile_2d)

# test data
if not os.path.exists(s_path_2d):
    pcdl.install_data()

print(f"process: pcdl pyCLI functions from the command line...")


# load
class TestPyCliCellDf(object):
    ''' test for pcdl.pyMCDS data loader, the complete data set. '''

    # timeseries
    # + collapse (true false) ok
    # timestep and timeseries:
    # + path nop
    # + microenv (true, _false_) ok
    # + settingxml (string, _none_, _false_) ok
    # + verbose (true, _false_) nop
    # + values (int) ok
    # + drop (cell_type oxygen) ok
    # + keep (cell_type oxygen) ok

    def test_pcdl_get_cell_df_timeseries(self):
        result = subprocess.run(['pcdl_get_cell_df', s_path_2d], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell.csv'))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timeseries_collapsed(self):
        result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--collapse', 'false'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        ls_opathfile = result.stderr.decode('UTF8').replace("['","").replace("']\n","").split("', '")
        print(ls_opathfile, len(ls_opathfile))
        assert len(ls_opathfile) == 25 and \
               ls_opathfile[0].endswith('data_timeseries_2d/output00000000_cell.csv') and \
               ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_cell.csv')
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timeseries_microenv(self):
        result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--microenv', 'false'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert not set(df_cell.columns).issuperset({'oxygen'})
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timeseries_settingxmlfalse(self):
        result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--settingxml', 'false'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert set(df_cell.columns).issuperset({'attack_rates_0'})
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timeseries_settingxmlnone(self):
        result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--settingxml', 'none'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert set(df_cell.columns).issuperset({'attack_rates_0'})
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timeseries_value(self):
        result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '2'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert  df_cell.shape == (24758, 41)
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timeseries_drop(self):
        result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--drop', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert not set(df_cell.columns).issuperset({'cell_type', 'oxygen'})
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timeseries_keep(self):
        result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--keep', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert set(df_cell.columns).issuperset({'cell_type', 'oxygen'}) and \
               df_cell.shape == (24758, 15)
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timestep(self):
        result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell.csv'))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timestep_microenv(self):
        result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '--microenv', 'false'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert not set(df_cell.columns).issuperset({'oxygen'})
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timestep_settingxmlfalse(self):
        result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '--settingxml', 'false'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert set(df_cell.columns).issuperset({'attack_rates_0'})
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timestep_settingxmlnone(self):
        result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '--settingxml', 'none'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert set(df_cell.columns).issuperset({'attack_rates_0'})
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timestep_value(self):
        result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '2'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert  df_cell.shape == (1099, 40)
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timestep_drop(self):
        result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '--drop', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert not set(df_cell.columns).issuperset({'cell_type', 'oxygen'})
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timestep_keep(self):
        result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '--keep', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert set(df_cell.columns).issuperset({'cell_type', 'oxygen'}) and \
               df_cell.shape == (1099, 14)
        os.remove(s_opathfile)
"""

"""
class TestPyCliConcDf(object):
    ''' test for pcdl.pyMCDS data loader, the complete data set. '''
    print(f"process: pcdl pyCLI functions from the command line...")

    # timeseries
    # + collapse (true false) ok

    # timestep and timeseries:
    # + verbose (true, _false_) nop
    # + values (int) ok
    # + drop (oxygen) ok
    # + keep (oxygen) ok

    def test_pcdl_get_conc_df_timeseries(self):
        result = subprocess.run(['pcdl_get_conc_df', s_path_2d], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_conc.csv'))
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_timeseries_collapsed(self):
        result = subprocess.run(['pcdl_get_conc_df', s_path_2d, '--collapse', 'false'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        ls_opathfile = result.stderr.decode('UTF8').replace("['","").replace("']\n","").split("', '")
        print(ls_opathfile, len(ls_opathfile))
        assert len(ls_opathfile) == 25 and \
               ls_opathfile[0].endswith('data_timeseries_2d/output00000000_conc.csv') and \
               ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_conc.csv')
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_get_conc_df_timeseries_value(self):
        result = subprocess.run(['pcdl_get_conc_df', s_path_2d, '2'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_conc = pd.read_csv(s_opathfile, index_col=0)
        assert  df_conc.shape == (3025, 9)
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_timeseries_drop(self):
        result = subprocess.run(['pcdl_get_conc_df', s_path_2d, '--drop', 'oxygen'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_conc = pd.read_csv(s_opathfile, index_col=0)
        assert not set(df_conc.columns).issuperset({'conc_type', 'oxygen'})
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_timeseries_keep(self):
        result = subprocess.run(['pcdl_get_conc_df', s_path_2d, '--keep', 'oxygen'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_conc = pd.read_csv(s_opathfile, index_col=0)
        assert set(df_conc.columns).issuperset({'oxygen'}) and \
               df_conc.shape == (3025, 9)
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_timestep(self):
        result = subprocess.run(['pcdl_get_conc_df', s_pathfile_2d], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_conc.csv'))
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_timestep_value(self):
        result = subprocess.run(['pcdl_get_conc_df', s_pathfile_2d, '2'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_conc = pd.read_csv(s_opathfile, index_col=0)
        assert  df_conc.shape == (121, 9)
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_timestep_drop(self):
        result = subprocess.run(['pcdl_get_conc_df', s_pathfile_2d, '--drop', 'oxygen'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_conc = pd.read_csv(s_opathfile, index_col=0)
        assert not set(df_conc.columns).issuperset({'oxygen'})
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_timestep_keep(self):
        result = subprocess.run(['pcdl_get_conc_df', s_pathfile_2d, '--keep', 'oxygen'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_conc = pd.read_csv(s_opathfile, index_col=0)
        assert set(df_conc.columns).issuperset({'oxygen'}) and \
               df_conc.shape == (121, 9)
        os.remove(s_opathfile)
"""

"""
class TestPyCliCellDfFeature(object):
    ''' tests for one pcdl.pyCli function. '''

    # timeseries collapsed:
    # + path (str) nop
    # + microenv (true, _false_) ok
    # + settingxml (string, _none_, _false_) ok
    # + verbose (true, _false_) nop
    # + values (int) ok
    # + drop (cell_type oxygen) ok
    # + keep (cell_type oxygen) ok
    # + allvalues (false _true_) ok

    def test_pcdl_get_cell_df_features_timeseries(self):
        result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_features_minmax.json'))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_features_timeseries_allvalues(self):
        result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d, '--allvalues', 'true'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_features_all.json'))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_features_timeseries_microenv(self):
        result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d, '--microenv', 'false'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert not set(d_feature.keys()).issuperset({'oxygen'}) and \
               (len(d_feature) == 80)
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_features_timeseries_settingxmlfalse(self):
        result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d, '--settingxml', 'false'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert set(d_feature.keys()).issuperset({'attack_rates_0'}) and \
               (len(d_feature) == 83)
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_features_timeseries_settingxmlnone(self):
        result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d, '--settingxml', 'none'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert set(d_feature.keys()).issuperset({'attack_rates_0'}) and \
               (len(d_feature) == 83)
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_features_timeseries_value(self):
        result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d, '2'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert (len(d_feature) == 28)
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_features_timeseries_drop(self):
        result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d, '--drop', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert not set(d_feature.keys()).issuperset({'cell_type', 'oxygen'}) and \
               (len(d_feature) == 81)
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_features_timeseries_keep(self):
        result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d, '--keep', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert set(d_feature.keys()).issuperset({'cell_type', 'oxygen'}) and \
               (len(d_feature) == 2)
        os.remove(s_opathfile)
"""

"""
class TestPyCliConcDfFeature(object):
    ''' tests for one pcdl.pyCli function. '''

    # timeseries collapsed:
    # + path (str) nop
    # + verbose (true, _false_) nop
    # + values (int) ok
    # + drop (oxygen) ok
    # + keep (oxygen) ok
    # + allvalues (false _true_) ok

    def test_pcdl_get_conc_df_features_timeseries(self):
        result = subprocess.run(['pcdl_get_conc_df_features', s_path_2d], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_conc_features_minmax.json'))
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_features_timeseries_allvalues(self):
        result = subprocess.run(['pcdl_get_conc_df_features', s_path_2d, '--allvalues', 'true'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_conc_features_all.json'))
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_features_timeseries_value(self):
        result = subprocess.run(['pcdl_get_conc_df_features', s_path_2d, '2'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert (len(d_feature) == 1)
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_features_timeseries_drop(self):
        result = subprocess.run(['pcdl_get_conc_df_features', s_path_2d, '--drop', 'conc_type', 'oxygen'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert not set(d_feature.keys()).issuperset({'oxygen'}) and \
               (len(d_feature) == 0)
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_features_timeseries_keep(self):
        result = subprocess.run(['pcdl_get_conc_df_features', s_path_2d, '--keep', 'conc_type', 'oxygen'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert set(d_feature.keys()).issuperset({'oxygen'}) and \
               (len(d_feature) == 1)
        os.remove(s_opathfile)


class TestPyCliGraphGml(object):
    ''' test from pcdl.pyCLI get_graph_gml. '''

    # timestep and timeseries:
    # + path nop
    # + microenv (true, false) ok
    # + settingxml (string, _none_, _false_) ok
    # + verbose (true, _false_) nop
    # + graph_type (neighbor, _attached_) ok
    # + edge_attr (true, _false_) ok
    # + node_attr (cell_type oxygen) ok

    def test_pcdl_get_graph_gml_timeseries_default(self):
        result = subprocess.run(['pcdl_get_graph_gml', s_path_2d], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        ls_opathfile = result.stderr.decode('UTF8').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert len(ls_opathfile) == 25 and \
               ls_opathfile[0].endswith('data_timeseries_2d/output00000000_neighbor.gml') and \
               ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_neighbor.gml')
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_get_graph_gml_timeseries_microenv(self):
        result = subprocess.run(['pcdl_get_graph_gml', s_path_2d, '--microenv', 'false'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        ls_opathfile = result.stderr.decode('UTF8').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert len(ls_opathfile) == 25 and \
               ls_opathfile[0].endswith('data_timeseries_2d/output00000000_neighbor.gml') and \
               ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_neighbor.gml')
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_get_graph_gml_timeseries_settingxmlfalse_nodeattr_one(self):
        result = subprocess.run(['pcdl_get_graph_gml', s_path_2d, '--settingxml', 'false', '--node_attr', 'attack_rates_0'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        ls_opathfile = result.stderr.decode('UTF8').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert len(ls_opathfile) == 25 and \
               ls_opathfile[0].endswith('data_timeseries_2d/output00000000_neighbor.gml') and \
               ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_neighbor.gml')
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_get_graph_gml_timeseries_settingxmlnone_nodeattr_one(self):
        result = subprocess.run(['pcdl_get_graph_gml', s_path_2d, '--settingxml', 'none', '--node_attr', 'attack_rates_0'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        ls_opathfile = result.stderr.decode('UTF8').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert len(ls_opathfile) == 25 and \
               ls_opathfile[0].endswith('data_timeseries_2d/output00000000_neighbor.gml') and \
               ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_neighbor.gml')
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_get_graph_gml_timeseries_graph_type(self):
        result = subprocess.run(['pcdl_get_graph_gml', s_path_2d, 'attached'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        ls_opathfile = result.stderr.decode('UTF8').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert len(ls_opathfile) == 25 and \
               ls_opathfile[0].endswith('data_timeseries_2d/output00000000_attached.gml') and \
               ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_attached.gml')
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_get_graph_gml_timeseries_edge_attr(self):
        result = subprocess.run(['pcdl_get_graph_gml', s_path_2d, '--edge_attr', 'false'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        ls_opathfile = result.stderr.decode('UTF8').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert len(ls_opathfile) == 25 and \
               ls_opathfile[0].endswith('data_timeseries_2d/output00000000_neighbor.gml') and \
               ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_neighbor.gml')
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_get_graph_gml_timeseries_nodeattr_one(self):
        result = subprocess.run(['pcdl_get_graph_gml', s_path_2d, '--node_attr', 'cell_type'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        ls_opathfile = result.stderr.decode('UTF8').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert len(ls_opathfile) == 25 and \
               ls_opathfile[0].endswith('data_timeseries_2d/output00000000_neighbor.gml') and \
               ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_neighbor.gml')
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_get_graph_gml_timeseries_nodeattr_many(self):
        result = subprocess.run(['pcdl_get_graph_gml', s_path_2d, '--node_attr', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        ls_opathfile = result.stderr.decode('UTF8').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert len(ls_opathfile) == 25 and \
               ls_opathfile[0].endswith('data_timeseries_2d/output00000000_neighbor.gml') and \
               ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_neighbor.gml')
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_get_graph_gml_timestep_default(self):
        result = subprocess.run(['pcdl_get_graph_gml', s_pathfile_2d], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_neighbor.gml'))
        os.remove(s_opathfile)

    def test_pcdl_get_graph_gml_timestep_microenv(self):
        result = subprocess.run(['pcdl_get_graph_gml', s_pathfile_2d, '--microenv', 'false'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_neighbor.gml'))
        os.remove(s_opathfile)

    def test_pcdl_get_graph_gml_timestep_settingxmlfalse_nodeattr_one(self):
        result = subprocess.run(['pcdl_get_graph_gml', s_pathfile_2d, '--settingxml', 'false', '--node_attr', 'attack_rates_0'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_neighbor.gml'))
        os.remove(s_opathfile)

    def test_pcdl_get_graph_gml_timestep_settingxmlnone_nodeattr_one(self):
        result = subprocess.run(['pcdl_get_graph_gml', s_pathfile_2d, '--settingxml', 'none', '--node_attr', 'attack_rates_0'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_neighbor.gml'))
        os.remove(s_opathfile)

    def test_pcdl_get_graph_gml_timestep_graph_type(self):
        result = subprocess.run(['pcdl_get_graph_gml', s_pathfile_2d, 'attached'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_attached.gml'))
        os.remove(s_opathfile)

    def test_pcdl_get_graph_gml_timestep_edge_attr(self):
        result = subprocess.run(['pcdl_get_graph_gml', s_pathfile_2d, '--edge_attr', 'false'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_neighbor.gml'))
        os.remove(s_opathfile)

    def test_pcdl_get_graph_gml_timestep_node_attr_one(self):
        result = subprocess.run(['pcdl_get_graph_gml', s_pathfile_2d, '--node_attr', 'cell_type'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_neighbor.gml'))
        os.remove(s_opathfile)

    def test_pcdl_get_graph_gml_timestep_node_attr_many(self):
        result = subprocess.run(['pcdl_get_graph_gml', s_pathfile_2d, '--node_attr', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_neighbor.gml'))
        os.remove(s_opathfile)


class TestPyCliUnitSe(object):
    ''' test from pcdl.pyCLI get_graph_gml. '''

    # timestep:
    # + path (pathfile, path) ok
    # + microenv (true, _false_) ok
    # + settingxml (string, _none_, _false_) ok
    # + verbose (true, _false_) nop

    def test_pcdl_get_unit_se_timeseries(self):
        result = subprocess.run(['pcdl_get_unit_se', s_path_2d], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_unit.csv'))
        os.remove(s_opathfile)

    def test_pcdl_get_unit_se_timestep(self):
        result = subprocess.run(['pcdl_get_unit_se', s_pathfile_2d], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_unit.csv'))
        os.remove(s_opathfile)

    def test_pcdl_get_unit_se_timestep_microenv(self):
        result = subprocess.run(['pcdl_get_unit_se', s_pathfile_2d, '--microenv', 'false'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert not set(df_cell.index).issuperset({'oxygen'})
        os.remove(s_opathfile)

    def test_pcdl_get_unit_se_timestep_settingxmlfalse(self):
        result = subprocess.run(['pcdl_get_unit_se', s_pathfile_2d, '--settingxml', 'false'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert set(df_cell.index).issuperset({'attack_rates_0'})
        os.remove(s_opathfile)

    def test_pcdl_get_unit_se_timestep_settingxmlnone(self):
        result = subprocess.run(['pcdl_get_unit_se', s_pathfile_2d, '--settingxml', 'none'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert set(df_cell.index).issuperset({'attack_rates_0'})
        os.remove(s_opathfile)


class TestPyCliVersion(object):
    # timestep:
    # + path (pathfile, path) ok
    # + verbose (true, _false_) nop

    def test_pcdl_get_version_timeseries(self):
        result = subprocess.run(['pcdl_get_version', s_path_2d], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_version = result.stderr.decode('UTF8')
        assert s_version.startswith('version:\nPhysiCell_')

    def test_pcdl_get_version_timestep(self):
        result = subprocess.run(['pcdl_get_version', s_pathfile_2d], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_version = result.stderr.decode('UTF8')
        assert s_version.startswith('version:\nPhysiCell_')

class TestPyCliMakeGif(object):
    # time series
    # + path nop
    # + interface (default, 'tiff')

    def test_pcdl_get_version_timeseries_default(self):
        s_path = subprocess.run(['pcdl_plot_scatter', s_path_2d], check=False, capture_output=True)
        s_path = s_path.stderr.decode('UTF8').replace('\n','')
        result = subprocess.run(['pcdl_make_gif', s_path], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert s_opathfile.endswith('data_timeseries_2d/cell_cell_type_z0.0/cell_cell_type_z0.0_jpeg.gif')
        shutil.rmtree(s_path)

    def test_pcdl_get_version_timeseries_interface(self):
        s_path = subprocess.run(['pcdl_plot_contour', s_path_2d, 'oxygen', '--ext', 'tiff'], check=False, capture_output=True)
        s_path = s_path.stderr.decode('UTF8').replace('\n','')
        result = subprocess.run(['pcdl_make_gif', s_path, 'tiff'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert s_opathfile.endswith('data_timeseries_2d/substrate_oxygen_z0.0/substrate_oxygen_z0.0_tiff.gif')
        shutil.rmtree(s_path)


class TestPyCliMakeMove(object):
    # time series
    # + path nop
    # + interface (default, 'tiff')
    # + farme (default, 'tiff')

    def test_pcdl_get_version_timeseries_default(self):
        s_path = subprocess.run(['pcdl_plot_scatter', s_path_2d], check=False, capture_output=True)
        s_path = s_path.stderr.decode('UTF8').replace('\n','')
        result = subprocess.run(['pcdl_make_movie', s_path], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert s_opathfile.endswith('data_timeseries_2d/cell_cell_type_z0.0/cell_cell_type_z0.0_jpeg12.mp4')
        shutil.rmtree(s_path)

    def test_pcdl_get_version_timeseries_interface(self):
        s_path = subprocess.run(['pcdl_plot_contour', s_path_2d, 'oxygen', '--ext', 'tiff'], check=False, capture_output=True)
        s_path = s_path.stderr.decode('UTF8').replace('\n','')
        result = subprocess.run(['pcdl_make_movie', s_path, 'tiff'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert s_opathfile.endswith('data_timeseries_2d/substrate_oxygen_z0.0/substrate_oxygen_z0.0_tiff12.mp4')
        shutil.rmtree(s_path)

    def test_pcdl_get_version_timeseries_farme(self):
        s_path = subprocess.run(['pcdl_plot_contour', s_path_2d, 'oxygen', '--ext', 'jpeg'], check=False, capture_output=True)
        s_path = s_path.stderr.decode('UTF8').replace('\n','')
        result = subprocess.run(['pcdl_make_movie', s_path, '--framerate', '9'], check=False, capture_output=True)
        #print(f'result: {result}\n')
        #print(f'result.stderr: {result.stderr}\n')
        s_opathfile = result.stderr.decode('UTF8').replace('\n','')
        assert s_opathfile.endswith('data_timeseries_2d/substrate_oxygen_z0.0/substrate_oxygen_z0.0_jpeg9.mp4')
        shutil.rmtree(s_path)
