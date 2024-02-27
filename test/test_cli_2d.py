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
#import capsys
import os
import pandas as pd
import pathlib
import pcdl
import pytest
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

# load
class TestPyCliCellDf(object):
    ''' test for pcdl.pyMCDS data loader, the complete data set. '''
    print(f"process: pcdl pyCLI functions from the command line...")

    # timeseries
    # + collapse (true false) ok

    # timestep and timeseries:
    # + microenv (true, _false_) ok
    # + settingxml (string, _none_, _false_) ok
    # + verbose (true, _false_) nop
    # + values (ok) ok
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


class TestPyCliConcDf(object):
    ''' test for pcdl.pyMCDS data loader, the complete data set. '''
    print(f"process: pcdl pyCLI functions from the command line...")

    # timeseries
    # + collapse (true false) ok

    # timestep and timeseries:
    # + verbose (true, _false_) nop
    # + values (ok) ok
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

