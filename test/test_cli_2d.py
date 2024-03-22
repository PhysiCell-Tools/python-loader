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
import anndata as ad
import json
import os
import pandas as pd
import pathlib
import pcdl
import shutil
import subprocess

# const
s_path_2d = str(pathlib.Path(pcdl.__file__).parent.resolve()/'data_timeseries_2d')
s_file_2d = 'output00000024.xml'
s_pathfile_2d = f'{s_path_2d}/{s_file_2d}'
print("s_path_2d", s_path_2d)
print("s_pathfile_2d", s_pathfile_2d)

# probe for  data
if not os.path.exists(s_path_2d):
    pcdl.install_data()

print(f"process: pcdl pyCLI functions from the command line...")


# test code
class TestPyCliAnndata(object):
    ''' tests for one pcdl.pyCli function. '''

    # timeseries
    # + collapse (true false) ok

    # timestep and timeseries:
    # + path nop
    # + customtype ([], _oncoprotein:str_) ok
    # + microenv (true, _false_) ok
    # + graph (true, _false_) ok
    # + physiboss (true, _false_) ok
    # + settingxml (string, _none_, _false_) ok
    # + verbose (true, _false_) nop
    # + values (int) ok
    # + drop (cell_type oxygen) ok
    # + keep (cell_type oxygen) ok
    # + scale (maxabs, _std_) ok

    # timeseries
    def test_pcdl_get_anndata_timeseries(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_path_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_maxabs.h5ad')) and \
               (ann.shape == (24758, 79)) and \
               (ann.obs.shape == (24758, 7)) and \
               (len(ann.obsp) == 0) and \
               (ann.var.shape == (79, 0)) and \
               (len(ann.uns) == 0)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timeseries_collapsed(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--collapse', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        ls_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace("['","").replace("']\n","").split("', '")
        assert len(ls_opathfile) == 25 and \
               ls_opathfile[0].endswith('data_timeseries_2d/output00000000_cell_maxabs.h5ad') and \
               ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_cell_maxabs.h5ad') and \
               os.path.exists(ls_opathfile[12])
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_get_anndata_timeseries_customtype(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--custom_type', 'oncoprotein:str'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_maxabs.h5ad')) and \
               (set(ann.obs.columns).issuperset({'oncoprotein'})) and \
               (ann.shape == (24758, 78)) and \
               (ann.obs.shape == (24758, 8)) and \
               (len(ann.obsp) == 0) and \
               (ann.var.shape == (78, 0)) and \
               (len(ann.uns) == 0)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timeseries_microenv(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--microenv', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_maxabs.h5ad')) and \
               (not set(ann.var_names).issuperset({'oxygen'})) and \
               (ann.shape == (24758, 76)) and \
               (ann.obs.shape == (24758, 7)) and \
               (len(ann.obsp) == 0) and \
               (ann.var.shape == (76, 0)) and \
               (len(ann.uns) == 0)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timeseries_graph(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--graph', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_maxabs.h5ad')) and \
               (ann.shape == (24758, 79)) and \
               (ann.obs.shape == (24758, 7)) and \
               (len(ann.obsp) == 0) and \
               (ann.var.shape == (79, 0)) and \
               (len(ann.uns) == 0)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timeseries_physiboss(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--physiboss', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_maxabs.h5ad')) and \
               (ann.shape == (24758, 79)) and \
               (ann.obs.shape == (24758, 7)) and \
               (len(ann.obsp) == 0) and \
               (ann.var.shape == (79, 0)) and \
               (len(ann.uns) == 0)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timeseries_settingxmlfalse(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--settingxml', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_maxabs.h5ad')) and \
               (set(ann.var_names).issuperset({'attack_rates_0'})) and \
               (ann.shape == (24758, 79)) and \
               (ann.obs.shape == (24758, 7)) and \
               (len(ann.obsp) == 0) and \
               (ann.var.shape == (79, 0)) and \
               (len(ann.uns) == 0)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timeseries_settingxmlnone(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--settingxml', 'none'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_maxabs.h5ad')) and \
               (set(ann.var_names).issuperset({'attack_rates_0'})) and \
               (ann.shape == (24758, 79)) and \
               (ann.obs.shape == (24758, 7)) and \
               (len(ann.obsp) == 0) and \
               (ann.var.shape == (79, 0)) and \
               (len(ann.uns) == 0)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timeseries_value(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '2'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_maxabs.h5ad')) and \
               (ann.shape == (24758, 26)) and \
               (ann.obs.shape == (24758, 5)) and \
               (len(ann.obsp) == 0) and \
               (ann.var.shape == (26, 0)) and \
               (len(ann.uns) == 0)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timeseries_drop(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--drop', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_maxabs.h5ad')) and \
               (not set(ann.var_names).issuperset({'cell_type'})) and \
               (not set(ann.obs_keys()).issuperset({'oxygen'})) and \
               (ann.shape == (24758, 78)) and \
               (ann.obs.shape == (24758, 6)) and \
               (len(ann.obsp) == 0) and \
               (ann.var.shape == (78, 0)) and \
               (len(ann.uns) == 0)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timeseries_keep(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--keep', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_maxabs.h5ad')) and \
               (set(ann.var_names).issuperset({'oxygen'})) and \
               (set(ann.obs_keys()).issuperset({'cell_type'})) and \
               (ann.shape == (24758, 1)) and \
               (ann.obs.shape == (24758, 4)) and \
               (len(ann.obsp) == 0) and \
               (ann.var.shape == (1, 0)) and \
               (len(ann.uns) == 0)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timeseries_scale(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--scale', 'std'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_std.h5ad')) and \
               (ann.shape == (24758, 79)) and \
               (ann.obs.shape == (24758, 7)) and \
               (len(ann.obsp) == 0) and \
               (ann.var.shape == (79, 0)) and \
               (len(ann.uns) == 0)
        os.remove(s_opathfile)

    # timestep
    def test_pcdl_get_anndata_timestep(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell_maxabs.h5ad')) and \
               (ann.shape == (1099, 79)) and \
               (ann.obs.shape == (1099, 6)) and \
               (len(ann.obsp) == 2) and \
               (ann.var.shape == (79, 0)) and \
               (len(ann.uns) == 1)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timestep_microenv(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '--microenv', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell_maxabs.h5ad')) and \
               (not set(ann.var_names).issuperset({'oxygen'})) and \
               (ann.shape == (1099, 76)) and \
               (ann.obs.shape == (1099, 6)) and \
               (len(ann.obsp) == 2) and \
               (ann.var.shape == (76, 0)) and \
               (len(ann.uns) == 1)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timestep_graph(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '--graph', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell_maxabs.h5ad')) and \
               (ann.shape == (1099, 79)) and \
               (ann.obs.shape == (1099, 6)) and \
               (len(ann.obsp) == 0) and \
               (ann.var.shape == (79, 0)) and \
               (len(ann.uns) == 0)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timestep_physiboss(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '--physiboss', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell_maxabs.h5ad')) and \
               (ann.shape == (1099, 79)) and \
               (ann.obs.shape == (1099, 6)) and \
               (len(ann.obsp) == 2) and \
               (ann.var.shape == (79, 0)) and \
               (len(ann.uns) == 1)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timestep_settingxmlfalse(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '--settingxml', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell_maxabs.h5ad')) and \
               (set(ann.var_names).issuperset({'attack_rates_0'})) and \
               (ann.shape == (1099, 79)) and \
               (ann.obs.shape == (1099, 6)) and \
               (len(ann.obsp) == 2) and \
               (ann.var.shape == (79, 0)) and \
               (len(ann.uns) == 1)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timestep_settingxmlnone(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '--settingxml', 'none'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell_maxabs.h5ad')) and \
               (set(ann.var_names).issuperset({'attack_rates_0'})) and \
               (ann.shape == (1099, 79)) and \
               (ann.obs.shape == (1099, 6)) and \
               (len(ann.obsp) == 2) and \
               (ann.var.shape == (79, 0)) and \
               (len(ann.uns) == 1)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timestep_value(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '2'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell_maxabs.h5ad')) and \
               (ann.shape == (1099, 26)) and \
               (ann.obs.shape == (1099, 4)) and \
               (len(ann.obsp) == 2) and \
               (ann.var.shape == (26, 0)) and \
               (len(ann.uns) == 1)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timestep_drop(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '--drop', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell_maxabs.h5ad')) and \
               (not set(ann.var_names).issuperset({'cell_type'})) and \
               (not set(ann.obs_keys()).issuperset({'oxygen'})) and \
               (ann.shape == (1099, 78)) and \
               (ann.obs.shape == (1099, 5)) and \
               (len(ann.obsp) == 2) and \
               (ann.var.shape == (78, 0)) and \
               (len(ann.uns) == 1)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timestep_keep(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '--keep', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell_maxabs.h5ad')) and \
               (set(ann.var_names).issuperset({'oxygen'})) and \
               (set(ann.obs_keys()).issuperset({'cell_type'})) and \
               (ann.shape == (1099, 1)) and \
               (ann.obs.shape == (1099, 3)) and \
               (len(ann.obsp) == 2) and \
               (ann.var.shape == (1, 0)) and \
               (len(ann.uns) == 1)
        os.remove(s_opathfile)

    def test_pcdl_get_anndata_timestep_scale(self):
        s_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '--scale', 'std'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        ann = ad.read_h5ad(s_opathfile)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell_std.h5ad')) and \
               (ann.shape == (1099, 79)) and \
               (ann.obs.shape == (1099, 6)) and \
               (len(ann.obsp) == 2) and \
               (ann.var.shape == (79, 0)) and \
               (len(ann.uns) == 1)
        os.remove(s_opathfile)


class TestPyCliCellDf(object):
    ''' tests for one pcdl.pyCli function. '''

    # timeseries
    # + collapse (true false) ok

    # timestep and timeseries:
    # + path nop
    # + customtype nop (because the dataframe is straightaway saved as csv)
    # + microenv (true, _false_) ok
    # + physiboss (true, _false_) ok
    # + settingxml (string, _none_, _false_) ok
    # + verbose (true, _false_) nop
    # + values (int) ok
    # + drop (cell_type oxygen) ok
    # + keep (cell_type oxygen) ok

    def test_pcdl_get_cell_df_timeseries(self):
        s_result = subprocess.run(['pcdl_get_cell_df', s_path_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell.csv')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timeseries_collapsed(self):
        s_result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--collapse', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        ls_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, len(ls_opathfile))
        assert len(ls_opathfile) == 25 and \
               ls_opathfile[0].endswith('data_timeseries_2d/output00000000_cell.csv') and \
               ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_cell.csv') and \
               os.path.exists(ls_opathfile[12])
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timeseries_microenv(self):
        s_result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--microenv', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell.csv')) and \
               (not set(df_cell.columns).issuperset({'oxygen'}))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timeseries_physiboss(self):
        s_result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--physiboss', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell.csv')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timeseries_settingxmlfalse(self):
        s_result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--settingxml', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell.csv')) and \
               (set(df_cell.columns).issuperset({'attack_rates_0'}))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timeseries_settingxmlnone(self):
        s_result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--settingxml', 'none'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell.csv')) and \
               (set(df_cell.columns).issuperset({'attack_rates_0'}))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timeseries_value(self):
        s_result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '2'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell.csv')) and \
               (df_cell.shape == (24758, 41))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timeseries_drop(self):
        s_result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--drop', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell.csv')) and \
               (not set(df_cell.columns).issuperset({'cell_type', 'oxygen'})) and \
               (df_cell.shape == (24758, 94))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timeseries_keep(self):
        s_result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--keep', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell.csv')) and \
               (set(df_cell.columns).issuperset({'cell_type', 'oxygen'})) and \
               (df_cell.shape == (24758, 15))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timestep(self):
        s_result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell.csv'))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timestep_microenv(self):
        s_result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '--microenv', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell.csv')) and \
               (not set(df_cell.columns).issuperset({'oxygen'}))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timestep_physiboss(self):
        s_result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '--physiboss', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell.csv'))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timestep_settingxmlfalse(self):
        s_result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '--settingxml', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell.csv')) and \
               (set(df_cell.columns).issuperset({'attack_rates_0'}))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timestep_settingxmlnone(self):
        s_result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '--settingxml', 'none'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell.csv')) and \
               (set(df_cell.columns).issuperset({'attack_rates_0'}))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timestep_value(self):
        s_result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '2'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell.csv')) and \
               (df_cell.shape == (1099, 40))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timestep_drop(self):
        s_result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '--drop', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell.csv')) and \
               (not set(df_cell.columns).issuperset({'cell_type', 'oxygen'})) and \
               (df_cell.shape == (1099, 93))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_timestep_keep(self):
        s_result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '--keep', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_cell.csv')) and \
               (set(df_cell.columns).issuperset({'cell_type', 'oxygen'})) and \
               (df_cell.shape == (1099, 14))
        os.remove(s_opathfile)


class TestPyCliConcDf(object):
    ''' tests for one pcdl.pyCli function. '''

    # timeseries
    # + collapse (true false) ok

    # timestep and timeseries:
    # + verbose (true, _false_) nop
    # + values (int) ok
    # + drop (oxygen) ok
    # + keep (oxygen) ok

    def test_pcdl_get_conc_df_timeseries(self):
        s_result = subprocess.run(['pcdl_get_conc_df', s_path_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_conc.csv')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_timeseries_collapsed(self):
        s_result = subprocess.run(['pcdl_get_conc_df', s_path_2d, '--collapse', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        ls_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, len(ls_opathfile))
        assert len(ls_opathfile) == 25 and \
               ls_opathfile[0].endswith('data_timeseries_2d/output00000000_conc.csv') and \
               ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_conc.csv') and \
               os.path.exists(ls_opathfile[12])
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_get_conc_df_timeseries_value(self):
        s_result = subprocess.run(['pcdl_get_conc_df', s_path_2d, '2'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_conc = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_conc.csv')) and \
               (df_conc.shape == (3025, 9))
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_timeseries_drop(self):
        s_result = subprocess.run(['pcdl_get_conc_df', s_path_2d, '--drop', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_conc = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_conc.csv')) and \
               (not set(df_conc.columns).issuperset({'oxygen'})) and \
               (df_conc.shape == (3025, 8))
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_timeseries_keep(self):
        s_result = subprocess.run(['pcdl_get_conc_df', s_path_2d, '--keep', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_conc = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_conc.csv')) and \
               (set(df_conc.columns).issuperset({'oxygen'})) and \
               (df_conc.shape == (3025, 9))
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_timestep(self):
        s_result = subprocess.run(['pcdl_get_conc_df', s_pathfile_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_conc.csv')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_timestep_value(self):
        s_result = subprocess.run(['pcdl_get_conc_df', s_pathfile_2d, '2'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_conc = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_conc.csv')) and \
               (df_conc.shape == (121, 9))
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_timestep_drop(self):
        s_result = subprocess.run(['pcdl_get_conc_df', s_pathfile_2d, '--drop', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_conc = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_conc.csv')) and \
               (not set(df_conc.columns).issuperset({'oxygen'})) and \
               (df_conc.shape == (121, 8))
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_timestep_keep(self):
        s_result = subprocess.run(['pcdl_get_conc_df', s_pathfile_2d, '--keep', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_conc = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_conc.csv')) and \
               (set(df_conc.columns).issuperset({'oxygen'})) and \
               (df_conc.shape == (121, 9))
        os.remove(s_opathfile)


class TestPyCliCellDfFeature(object):
    ''' tests for one pcdl.pyCli function. '''

    # timeseries collapsed:
    # + path (str) nop
    # + customtype ([], _oncoprotein:str_) ok
    # + microenv (true, _false_) ok
    # + physiboss (true, _false_)
    # + settingxml (string, _none_, _false_) ok
    # + verbose (true, _false_) nop
    # + values (int) ok
    # + drop (cell_type oxygen) ok
    # + keep (cell_type oxygen) ok
    # + allvalues (false _true_) ok

    def test_pcdl_get_cell_df_features_timeseries(self):
        s_result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_features_minmax.json')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_features_timeseries_customtype(self):
        s_result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d, '--custom_type', 'oncoprotein:str'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        d_cell = json.load(open(s_opathfile))
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_features_minmax.json')) and \
               (d_cell['oncoprotein'] == ['0','1','2'])
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_features_timeseries_microenv(self):
        s_result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d, '--microenv', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_features_minmax.json')) and \
               (not set(d_feature.keys()).issuperset({'oxygen'})) and \
               (len(d_feature) == 80)
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_features_timeseries_physiboss(self):
        s_result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d, '--physiboss', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_features_minmax.json')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_features_timeseries_settingxmlfalse(self):
        s_result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d, '--settingxml', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_features_minmax.json')) and \
               (set(d_feature.keys()).issuperset({'attack_rates_0'})) and \
               (len(d_feature) == 83)
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_features_timeseries_settingxmlnone(self):
        s_result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d, '--settingxml', 'none'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_features_minmax.json')) and \
               (set(d_feature.keys()).issuperset({'attack_rates_0'})) and \
               (len(d_feature) == 83)
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_features_timeseries_value(self):
        s_result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d, '2'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_features_minmax.json')) and \
               (len(d_feature) == 28)
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_features_timeseries_drop(self):
        s_result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d, '--drop', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_features_minmax.json')) and \
               (not set(d_feature.keys()).issuperset({'cell_type', 'oxygen'})) and \
               (len(d_feature) == 81)
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_features_timeseries_keep(self):
        s_result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d, '--keep', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_features_minmax.json')) and \
               (set(d_feature.keys()).issuperset({'cell_type', 'oxygen'})) and \
               (len(d_feature) == 2)
        os.remove(s_opathfile)

    def test_pcdl_get_cell_df_features_timeseries_allvalues(self):
        s_result = subprocess.run(['pcdl_get_cell_df_features', s_path_2d, '--allvalues', 'true'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_cell_features_all.json')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)


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
        s_result = subprocess.run(['pcdl_get_conc_df_features', s_path_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_conc_features_minmax.json')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_features_timeseries_value(self):
        s_result = subprocess.run(['pcdl_get_conc_df_features', s_path_2d, '2'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_conc_features_minmax.json')) and \
               (len(d_feature) == 1)
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_features_timeseries_drop(self):
        s_result = subprocess.run(['pcdl_get_conc_df_features', s_path_2d, '--drop', 'conc_type', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_conc_features_minmax.json')) and \
               (not set(d_feature.keys()).issuperset({'oxygen'})) and \
               (len(d_feature) == 0)
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_features_timeseries_keep(self):
        s_result = subprocess.run(['pcdl_get_conc_df_features', s_path_2d, '--keep', 'conc_type', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        d_feature = json.load(open(s_opathfile))
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_conc_features_minmax.json')) and \
               (set(d_feature.keys()).issuperset({'oxygen'})) and \
               (len(d_feature) == 1)
        os.remove(s_opathfile)

    def test_pcdl_get_conc_df_features_timeseries_allvalues(self):
        s_result = subprocess.run(['pcdl_get_conc_df_features', s_path_2d, '--allvalues', 'true'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_conc_features_all.json')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)


class TestPyCliGraphGml(object):
    ''' tests for one pcdl.pyCli function. '''

    # timestep and timeseries:
    # + path nop
    # + customtype ([], _oncoprotein:str_) ok
    # + microenv (true, false) ok
    # + physiboss (true, _false_)
    # + settingxml (string, _none_, _false_) ok
    # + verbose (true, _false_) nop
    # + graph_type (neighbor, _attached_) ok
    # + edge_attr (true, _false_) ok
    # + node_attr (cell_type oxygen) ok

    def test_pcdl_make_graph_gml_timeseries_default(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        ls_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert (len(ls_opathfile) == 25) and \
               (ls_opathfile[0].endswith('data_timeseries_2d/output00000000_neighbor.gml')) and \
               (ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(ls_opathfile[12]))
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_make_grah_gml_timeseries_customtype_nodeattr_one(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, '--custom_type', 'oncoprotein:str', '--node_attr', 'oncoprotein'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        ls_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert (len(ls_opathfile) == 25) and \
               (ls_opathfile[0].endswith('data_timeseries_2d/output00000000_neighbor.gml')) and \
               (ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(ls_opathfile[12]))
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timeseries_microenv(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, '--microenv', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        ls_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert (len(ls_opathfile) == 25) and \
               (ls_opathfile[0].endswith('data_timeseries_2d/output00000000_neighbor.gml')) and \
               (ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(ls_opathfile[12]))
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timeseries_physiboss(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, '--physiboss', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        ls_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert (len(ls_opathfile) == 25) and \
               (ls_opathfile[0].endswith('data_timeseries_2d/output00000000_neighbor.gml')) and \
               (ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(ls_opathfile[12]))
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timeseries_settingxmlfalse_nodeattr_one(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, '--settingxml', 'false', '--node_attr', 'attack_rates_0'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        ls_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert (len(ls_opathfile) == 25) and \
               (ls_opathfile[0].endswith('data_timeseries_2d/output00000000_neighbor.gml')) and \
               (ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(ls_opathfile[12]))
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timeseries_settingxmlnone_nodeattr_one(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, '--settingxml', 'none', '--node_attr', 'attack_rates_0'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        ls_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert (len(ls_opathfile) == 25) and \
               (ls_opathfile[0].endswith('data_timeseries_2d/output00000000_neighbor.gml')) and \
               (ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(ls_opathfile[12]))
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timeseries_graph_type(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, 'attached'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        ls_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert (len(ls_opathfile) == 25) and \
               (ls_opathfile[0].endswith('data_timeseries_2d/output00000000_attached.gml')) and \
               (ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_attached.gml')) and \
               (os.path.exists(ls_opathfile[12]))
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timeseries_edge_attr(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, '--edge_attr', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        ls_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert (len(ls_opathfile) == 25) and \
               (ls_opathfile[0].endswith('data_timeseries_2d/output00000000_neighbor.gml')) and \
               (ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(ls_opathfile[12]))
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timeseries_nodeattr_one(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, '--node_attr', 'cell_type'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        ls_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert (len(ls_opathfile) == 25) and \
               (ls_opathfile[0].endswith('data_timeseries_2d/output00000000_neighbor.gml')) and \
               (ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(ls_opathfile[12]))
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timeseries_nodeattr_many(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, '--node_attr', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        ls_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace("['","").replace("']\n","").split("', '")
        #print(ls_opathfile, ls_opathfile)
        assert (len(ls_opathfile) == 25) and \
               (ls_opathfile[0].endswith('data_timeseries_2d/output00000000_neighbor.gml')) and \
               (ls_opathfile[-1].endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(ls_opathfile[12]))
        for s_opathfile in ls_opathfile:
            os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timestep_default(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timestep_customtype_nodeattr_one(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, '--custom_type', 'oncoprotein:str', '--node_attr', 'oncoprotein'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timestep_microenv(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, '--microenv', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timestep_physiboss(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, '--physiboss', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timestep_settingxmlfalse_nodeattr_one(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, '--settingxml', 'false', '--node_attr', 'attack_rates_0'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timestep_settingxmlnone_nodeattr_one(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, '--settingxml', 'none', '--node_attr', 'attack_rates_0'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timestep_graph_type(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, 'attached'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_attached.gml')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timestep_edge_attr(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, '--edge_attr', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timestep_node_attr_one(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, '--node_attr', 'cell_type'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_make_graph_gml_timestep_node_attr_many(self):
        s_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, '--node_attr', 'cell_type', 'oxygen'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/output00000024_neighbor.gml')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)


class TestPyCliParameterDict(object):
    ''' tests for one pcdl.pyCli function. '''

    # timestep:
    # + path (pathfile, path) ok
    # + microenv (true, _false_) ok
    # + verbose (true, _false_) nop

    def test_pcdl_get_parameter_dict_timeseries(self):
        s_result = subprocess.run(['pcdl_get_parameter_dict', s_path_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        d_parameter = json.load(open(s_opathfile))
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_parameter.json')) and \
               (len(d_parameter) == 63)
        os.remove(s_opathfile)

    def test_pcdl_get_parameter_dict_timestep(self):
        s_result = subprocess.run(['pcdl_get_parameter_dict', s_pathfile_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        d_parameter = json.load(open(s_opathfile))
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_parameter.json')) and \
               (len(d_parameter) == 63)
        os.remove(s_opathfile)

    def test_pcdl_get_parameter_dict_timestep_microenv(self):
        s_result = subprocess.run(['pcdl_get_parameter_dict', s_pathfile_2d, '--microenv', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        d_parameter = json.load(open(s_opathfile))
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_parameter.json')) and \
               (len(d_parameter) == 63)
        os.remove(s_opathfile)


class TestPyCliRuleDf(object):
    ''' tests for one pcdl.pyCli function. '''

    # timestep:
    # + path (pathfile, path) ok
    # + verbose (true, _false_) nop

    def test_pcdl_get_rule_df_timeseries(self):
        s_result = subprocess.run(['pcdl_get_rule_df', s_path_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('None'))

    def test_pcdl_get_rule_df_timestep(self):
        s_result = subprocess.run(['pcdl_get_rule_df', s_pathfile_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('None'))


class TestPyCliUnitDict(object):
    ''' tests for one pcdl.pyCli function. '''

    # timestep:
    # + path (pathfile, path) ok
    # + microenv (true, _false_) ok
    # + settingxml (string, _none_, _false_) ok
    # + verbose (true, _false_) nop

    def test_pcdl_get_unit_dict_timeseries(self):
        s_result = subprocess.run(['pcdl_get_unit_dict', s_path_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_unit.csv')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_get_unit_dict_timestep(self):
        s_result = subprocess.run(['pcdl_get_unit_dict', s_pathfile_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_unit.csv')) and \
               (os.path.exists(s_opathfile))
        os.remove(s_opathfile)

    def test_pcdl_get_unit_dict_timestep_microenv(self):
        s_result = subprocess.run(['pcdl_get_unit_dict', s_pathfile_2d, '--microenv', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_unit.csv')) and \
               (not set(df_cell.index).issuperset({'oxygen'}))
        os.remove(s_opathfile)

    def test_pcdl_get_unit_dict_timestep_settingxmlfalse(self):
        s_result = subprocess.run(['pcdl_get_unit_dict', s_pathfile_2d, '--settingxml', 'false'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_unit.csv')) and \
               (set(df_cell.index).issuperset({'attack_rates_0'}))
        os.remove(s_opathfile)

    def test_pcdl_get_unit_dict_timestep_settingxmlnone(self):
        s_result = subprocess.run(['pcdl_get_unit_dict', s_pathfile_2d, '--settingxml', 'none'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        df_cell = pd.read_csv(s_opathfile, index_col=0)
        assert (s_opathfile.endswith('data_timeseries_2d/timeseries_unit.csv')) and \
               (set(df_cell.index).issuperset({'attack_rates_0'}))
        os.remove(s_opathfile)


class TestPyCliVersion(object):
    ''' tests for one pcdl.pyCli function. '''

    # timestep:
    # + path (pathfile, path) ok
    # + verbose (true, _false_) nop

    def test_pcdl_get_version_timeseries(self):
        s_result = subprocess.run(['pcdl_get_version', s_path_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_version = s_result.stderr.decode('UTF8').replace('\r','')
        assert (s_version.startswith('version:\nPhysiCell_'))

    def test_pcdl_get_version_timestep(self):
        s_result = subprocess.run(['pcdl_get_version', s_pathfile_2d], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_version = s_result.stderr.decode('UTF8').replace('\r','')
        assert (s_version.startswith('version:\nPhysiCell_'))


class TestPyCliMakeGif(object):
    ''' tests for one pcdl.pyCli function. '''

    # time series
    # + path nop
    # + interface (default, 'tiff')

    def test_pcdl_make_gif_timeseries_default(self):
        s_path = subprocess.run(['pcdl_plot_scatter', s_path_2d], check=False, capture_output=True)
        s_path = s_path.stderr.decode('UTF8').replace('\r','').replace('\n','')
        s_result = subprocess.run(['pcdl_make_gif', s_path], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').split('\n')[-2]
        assert (s_opathfile.endswith('data_timeseries_2d/cell_cell_type_z0.0/cell_cell_type_z0.0_jpeg.gif')) and \
               (os.path.exists(s_opathfile))
        shutil.rmtree(s_path)

    def test_pcdl_make_gif_timeseries_interface(self):
        s_path = subprocess.run(['pcdl_plot_contour', s_path_2d, 'oxygen', '--ext', 'tiff'], check=False, capture_output=True)
        s_path = s_path.stderr.decode('UTF8').replace('\r','').replace('\n','')
        s_result = subprocess.run(['pcdl_make_gif', s_path, 'tiff'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').split('\n')[-2]
        assert (s_opathfile.endswith('data_timeseries_2d/conc_oxygen_z0.0/conc_oxygen_z0.0_tiff.gif')) and \
               (os.path.exists(s_opathfile))
        shutil.rmtree(s_path)


class TestPyCliMakeMove(object):
    ''' tests for one pcdl.pyCli function. '''

    # time series
    # + path nop
    # + interface (default, 'tiff')
    # + farme (default, 'tiff')

    def test_pcdl_make_movie_timeseries_default(self):
        s_path = subprocess.run(['pcdl_plot_scatter', s_path_2d], check=False, capture_output=True)
        s_path = s_path.stderr.decode('UTF8').replace('\r','').replace('\n','')
        s_result = subprocess.run(['pcdl_make_movie', s_path], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').split('\n')[-2]
        assert (s_opathfile.endswith('data_timeseries_2d/cell_cell_type_z0.0/cell_cell_type_z0.0_jpeg12.mp4')) and \
               (os.path.exists(s_opathfile))
        shutil.rmtree(s_path)

    def test_pcdl_make_movie_timeseries_interface(self):
        s_path = subprocess.run(['pcdl_plot_contour', s_path_2d, 'oxygen', '--ext', 'tiff'], check=False, capture_output=True)
        s_path = s_path.stderr.decode('UTF8').replace('\r','').replace('\n','')
        s_result = subprocess.run(['pcdl_make_movie', s_path, 'tiff'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').split('\n')[-2]
        assert (s_opathfile.endswith('data_timeseries_2d/conc_oxygen_z0.0/conc_oxygen_z0.0_tiff12.mp4')) and \
               (os.path.exists(s_opathfile))
        shutil.rmtree(s_path)

    def test_pcdl_make_movie_timeseries_farme(self):
        s_path = subprocess.run(['pcdl_plot_contour', s_path_2d, 'oxygen', '--ext', 'jpeg'], check=False, capture_output=True)
        s_path = s_path.stderr.decode('UTF8').replace('\r','').replace('\n','')
        s_result = subprocess.run(['pcdl_make_movie', s_path, '--framerate', '9'], check=False, capture_output=True)
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').split('\n')[-2]
        assert (s_opathfile.endswith('data_timeseries_2d/conc_oxygen_z0.0/conc_oxygen_z0.0_jpeg9.mp4')) and \
               (os.path.exists(s_opathfile))
        shutil.rmtree(s_path)


class TestPyCliPlotContour(object):
    ''' tests for one pcdl.pyCli function. '''

    # time series and time steps.
    # + path nop
    # + verbose (true, _false_) nop
    # + focus (_oxygen_) ok
    # + z_slize (0.0, _1.1_) ok
    # + extrema (none, _0.0 40.0_) ok
    # + alpha (1.0, _0.2_) ok
    # + fill (true, _false_)
    # + cmp (viridis, _magma_)
    # + title (, _abc_)
    # + grid (true, _false_)
    # + xlim (none, _-40 400.0_)
    # + ylim (none, _-30 300_)
    # + xyequal (true, _false_)
    # + figsizepx (none, _[641, 481]_)
    # + ext (jpeg, _tiff_)
    # + figbgcolor (none, _yellow_)

    def test_pcdl_plot_contour_default(self):
        s_result = subprocess.run([
            'pcdl_plot_contour', s_pathfile_2d, 'oxygen'
        ], check=False, capture_output=True)
        #print(f'\ns_result: {s_result}\n')
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_stdout = s_result.stdout.decode('UTF8').replace('\r','')
        s_opath = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (os.path.exists(s_opath + 'output00000024_oxygen.jpeg')) and \
           (s_opath.endswith('pcdl/data_timeseries_2d/conc_oxygen_z0.0/')) and \
           (s_stdout.find("focus='oxygen'") > -1) and \
           (s_stdout.find("z_slice=0.0") > -1) and \
           (s_stdout.find("extrema=['none']") > -1) and \
           (s_stdout.find("alpha=1.0") > -1) and \
           (s_stdout.find("fill='true'") > -1) and \
           (s_stdout.find("cmap='viridis'") > -1) and \
           (s_stdout.find("title=''") > -1) and \
           (s_stdout.find("grid='true'") > -1) and \
           (s_stdout.find("xlim=['none']") > -1) and \
           (s_stdout.find("ylim=['none']") > -1) and \
           (s_stdout.find("xyequal='true'") > -1) and \
           (s_stdout.find("figsizepx=['none']") > -1) and \
           (s_stdout.find("ext='jpeg'") > -1) and \
           (s_stdout.find("figbgcolor='none'") > -1)
        shutil.rmtree(s_opath)

    def test_pcdl_plot_contour_set(self):
        s_result = subprocess.run([
            'pcdl_plot_contour', s_pathfile_2d, 'oxygen',
            '--z_slice', '1.1',
            '--extrema', '0.0', '40.0',
            '--alpha', '0.5',
            '--fill', 'false',
            '--cmap', 'magma',
            '--title', 'abc',
            '--grid', 'false',
            '--xlim', '-40', '400',
            '--ylim', '-30', '300',
            '--xyequal', 'false',
            '--figsizepx', '842', '531',
            '--ext', 'tiff',
            '--figbgcolor', 'yellow',
        ], check=False, capture_output=True)
        #print(f'\ns_result: {s_result}\n')
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_stdout = s_result.stdout.decode('UTF8').replace('\r','')
        s_opath = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (os.path.exists(s_opath + 'output00000024_oxygen.tiff')) and \
           (s_opath.endswith('pcdl/data_timeseries_2d/conc_oxygen_z0.0/')) and \
           (s_stdout.find("focus='oxygen'") > -1) and \
           (s_stdout.find("z_slice=1.1") > -1) and \
           (s_stdout.find("extrema=['0.0', '40.0']") > -1) and \
           (s_stdout.find("alpha=0.5") > -1) and \
           (s_stdout.find("fill='false'") > -1) and \
           (s_stdout.find("cmap='magma'") > -1) and \
           (s_stdout.find("title='abc'") > -1) and \
           (s_stdout.find("grid='false'") > -1) and \
           (s_stdout.find("xlim=['-40', '400']") > -1) and \
           (s_stdout.find("ylim=['-30', '300']") > -1) and \
           (s_stdout.find("xyequal='false'") > -1) and \
           (s_stdout.find("figsizepx=['842', '531']") > -1) and \
           (s_stdout.find("ext='tiff'") > -1) and \
           (s_stdout.find("figbgcolor='yellow'") > -1)
        shutil.rmtree(s_opath)


class TestPyCliPlotScatter(object):
    ''' tests for one pcdl.pyCli function. '''

    # time series and time steps.
    # + path nop
    # + customtype ([], _oncoprotein:str_) ok
    # + microenv (true, _false_)
    # + physiboss (true, _false_)
    # + settingxml (PhysiCell_settings.xml, _false_)
    # + verbose (true, _false_) nop
    # + focus (_oxygen_) ok
    # + z_slize (0.0, _1.1_) ok
    # + z_axis (none, _0.0_40.0_) ok
    # + alpha (1.0, _0.5_) ok
    # + cmap (viridis, _magma_)
    # + title (, _abc_)
    # + grid (true, _false_)
    # + legend_loc ('lower left', _'upper right'_)
    # + xlim (none, _-40_400_)
    # + ylim (none, _-30_300_)
    # + xyequal (true, _false_)
    # + s ('none', '74')
    # + figsizepx (none, _[641, 481]_)
    # + ext (jpeg, _tiff_)
    # + figbgcolor (none, _yellow_)

    def test_pcdl_plot_scatter_default(self):
        s_result = subprocess.run([
            'pcdl_plot_scatter', s_pathfile_2d,
        ], check=False, capture_output=True)
        #print(f'\ns_result: {s_result}\n')
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_stdout = s_result.stdout.decode('UTF8').replace('\r','')
        s_opath = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (os.path.exists(s_opath + 'output00000024_cell_type.jpeg')) and \
           (s_opath.endswith('pcdl/data_timeseries_2d/cell_cell_type_z0.0/')) and \
           (s_stdout.find("custom_type=[]") > -1) and \
           (s_stdout.find("microenv='true'") > -1) and \
           (s_stdout.find("physiboss='true'") > -1) and \
           (s_stdout.find("settingxml='PhysiCell_settings.xml'") > -1) and \
           (s_stdout.find("focus='cell_type'") > -1) and \
           (s_stdout.find("z_slice=0.0") > -1) and \
           (s_stdout.find("z_axis=['none']") > -1) and \
           (s_stdout.find("alpha=1.0") > -1) and \
           (s_stdout.find("cmap='viridis'") > -1) and \
           (s_stdout.find("title=''") > -1) and \
           (s_stdout.find("grid='true'") > -1) and \
           (s_stdout.find("legend_loc='lower left'") > -1) and \
           (s_stdout.find("xlim=['none']") > -1) and \
           (s_stdout.find("ylim=['none']") > -1) and \
           (s_stdout.find("xyequal='true'") > -1) and \
           (s_stdout.find("s='none'") > -1) and \
           (s_stdout.find("figsizepx=['none']") > -1) and \
           (s_stdout.find("ext='jpeg'") > -1) and \
           (s_stdout.find("figbgcolor='none'") > -1)
        shutil.rmtree(s_opath)

    def test_pcdl_plot_scatter_set(self):
        s_result = subprocess.run([
            'pcdl_plot_scatter', s_pathfile_2d, 'oxygen',
            '--custom_type', 'oncoprotein:str',
            '--microenv', 'True',
            '--physiboss', 'false',
            '--settingxml', 'false',
            '--z_slice', '1.1',
            '--z_axis', '0.0', '40.0',
            '--alpha', '0.5',
            '--cmap', 'magma',
            '--title', 'abc',
            '--grid', 'false',
            '--legend_loc', 'upper right',
            '--xlim', '-40', '400',
            '--ylim', '-30', '300',
            '--xyequal', 'false',
            '--s', '74',
            '--figsizepx', '842', '531',
            '--ext', 'tiff',
            '--figbgcolor', 'yellow',
        ], check=False, capture_output=True)
        #print(f'\ns_result: {s_result}\n')
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_stdout = s_result.stdout.decode('UTF8').replace('\r','')
        s_opath = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (os.path.exists(s_opath + 'output00000024_oxygen.tiff')) and \
           (s_opath.endswith('pcdl/data_timeseries_2d/cell_oxygen_z0.0/')) and \
           (s_stdout.find("custom_type=['oncoprotein:str']") > -1) and \
           (s_stdout.find("microenv='True'") > -1) and \
           (s_stdout.find("physiboss='false'") > -1) and \
           (s_stdout.find("settingxml='false'") > -1) and \
           (s_stdout.find("focus='oxygen'") > -1) and \
           (s_stdout.find("z_slice=1.1") > -1) and \
           (s_stdout.find("z_axis=['0.0', '40.0']") > -1) and \
           (s_stdout.find("alpha=0.5") > -1) and \
           (s_stdout.find("cmap='magma'") > -1) and \
           (s_stdout.find("title='abc'") > -1) and \
           (s_stdout.find("grid='false'") > -1) and \
           (s_stdout.find("legend_loc='upper right'") > -1) and \
           (s_stdout.find("xlim=['-40', '400']") > -1) and \
           (s_stdout.find("ylim=['-30', '300']") > -1) and \
           (s_stdout.find("xyequal='false'") > -1) and \
           (s_stdout.find("s='74'") > -1) and \
           (s_stdout.find("figsizepx=['842', '531']") > -1) and \
           (s_stdout.find("ext='tiff'") > -1) and \
           (s_stdout.find("figbgcolor='yellow'") > -1)
        shutil.rmtree(s_opath)


class TestPyCliPlotTimeSeries(object):
    ''' tests for one pcdl.pyCli function. '''

    # time series.
    # + path nop
    # + customtype ([], _oncoprotein:str_) ok
    # + microenv (true, _false_)
    # + physiboss (true, _false_)
    # + settingxml (PhysiCell_settings.xml, _false_)
    # + verbose (true, _false_) nop
    # + focus_cat ('none', _cell_type_) ok
    # + focus_num ('none', _oxygen_) ok
    # + aggregate_num ('mean', 'entropy') ok
    # + frame ('cell', 'conc')
    # + z_slice ('none', _1.1_) ok
    # + logy (false, _true_)
    # + ylim (['none'], ['', ''])
    # + secondary_y (false, _true_)
    # + subplots (false, _true_)
    # + sharex (false, _true_)
    # + sharey (false, _true_)
    # + linestyle ('-', '-.')
    # + linewidth ('none', 9)
    # + cmap ('none', 'magma')
    # + color ('none', 'maroon')
    # + grid (true, _false_)
    # + legend (true, _false_)
    # + yunit ('none', 'myunit')
    # + title ('none', 'my title')
    # + figsizepx (none, _[641, 481]_)
    # + ext (jpeg, _tiff_)
    # + figbgcolor (none, _yellow_)

    def test_pcdl_plot_timeseries_default(self):
        s_result = subprocess.run([
            'pcdl_plot_timeseries', s_path_2d, '-v', 'false'
        ], check=False, capture_output=True)
        #print(f'\ns_result: {s_result}\n')
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_stdout = s_result.stdout.decode('UTF8').replace('\r','')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (os.path.exists(s_opathfile)) and \
           (s_opathfile.endswith('pcdl/data_timeseries_2d/timeseries_cell_total_count.jpeg')) and \
           (s_stdout.find("custom_type=[]") > -1) and \
           (s_stdout.find("microenv='true'") > -1) and \
           (s_stdout.find("physiboss='true'") > -1) and \
           (s_stdout.find("settingxml='PhysiCell_settings.xml'") > -1) and \
           (s_stdout.find("focus_cat='none'") > -1) and \
           (s_stdout.find("focus_num='none'") > -1) and \
           (s_stdout.find("aggregate_num='mean'") > -1) and \
           (s_stdout.find("frame='cell'") > -1) and \
           (s_stdout.find("z_slice='none'") > -1) and \
           (s_stdout.find("logy='false'") > -1) and \
           (s_stdout.find("ylim=['none']") > -1) and \
           (s_stdout.find("secondary_y=['false']") > -1) and \
           (s_stdout.find("subplots='false'") > -1) and \
           (s_stdout.find("sharex='false'") > -1) and \
           (s_stdout.find("sharey='false'") > -1) and \
           (s_stdout.find("linestyle='-'") > -1) and \
           (s_stdout.find("linewidth='none'") > -1) and \
           (s_stdout.find("cmap='none'") > -1) and \
           (s_stdout.find("color=['none']") > -1) and \
           (s_stdout.find("grid='true'") > -1) and \
           (s_stdout.find("legend='true'") > -1) and \
           (s_stdout.find("yunit='none'") > -1) and \
           (s_stdout.find("title='none'") > -1) and \
           (s_stdout.find("figsizepx=['640', '480']") > -1) and \
           (s_stdout.find("ext='jpeg'") > -1) and \
           (s_stdout.find("figbgcolor='none'") > -1)
        os.remove(s_opathfile)

    def test_pcdl_plot_timeseries_set(self):
        s_result = subprocess.run([
            'pcdl_plot_timeseries', s_path_2d, 'None', 'oxygen', 'entropy', '-v', 'false',
            '--custom_type', 'oncoprotein:str',
            '--microenv', 'True',
            '--physiboss', 'false',
            '--settingxml', 'false',
            '--frame', 'conc',
            '--z_slice', '1.1',
            '--logy', 'true',
            '--ylim', '6.0', '7.0',
            '--secondary_y', 'true', 'abc', 'def',
            '--subplots', 'true',
            '--sharex', 'true',
            '--sharey', 'true',
            '--linestyle', ':',
            '--linewidth', '9',
            '--cmap', 'magma',
            '--color', 'maroon', 'orange', 'yellow',
            '--grid', 'false',
            '--legend', 'reverse',
            '--yunit', 'myunit',
            '--title', 'my title',
            '--figsizepx', '641', '481',
            '--ext', 'tiff',
            '--figbgcolor', 'cyan',
        ], check=False, capture_output=True)
        #print(f'\ns_result: {s_result}\n')
        #print(f'\ns_result.stdout: {s_result.stdout}\n')
        #print(f'\ns_result.stderr: {s_result.stderr}\n')
        s_stdout = s_result.stdout.decode('UTF8').replace('\r','')
        s_opathfile = s_result.stderr.decode('UTF8').replace('\r','').replace('\n','')
        assert (os.path.exists(s_opathfile)) and \
           (s_opathfile.endswith('pcdl/data_timeseries_2d/timeseries_conc_total_oxygen_entropy.tiff')) and \
           (s_stdout.find("custom_type=['oncoprotein:str']") > -1) and \
           (s_stdout.find("microenv='True'") > -1) and \
           (s_stdout.find("physiboss='false'") > -1) and \
           (s_stdout.find("settingxml='false'") > -1) and \
           (s_stdout.find("focus_cat='None'") > -1) and \
           (s_stdout.find("focus_num='oxygen'") > -1) and \
           (s_stdout.find("aggregate_num='entropy'") > -1) and \
           (s_stdout.find("frame='conc'") > -1) and \
           (s_stdout.find("z_slice='1.1'") > -1) and \
           (s_stdout.find("logy='true'") > -1) and \
           (s_stdout.find("ylim=['6.0', '7.0']") > -1) and \
           (s_stdout.find("secondary_y=['true', 'abc', 'def']") > -1) and \
           (s_stdout.find("subplots='true'") > -1) and \
           (s_stdout.find("sharex='true'") > -1) and \
           (s_stdout.find("sharey='true'") > -1) and \
           (s_stdout.find("linestyle=':'") > -1) and \
           (s_stdout.find("linewidth='9'") > -1) and \
           (s_stdout.find("cmap='magma'") > -1) and \
           (s_stdout.find("color=['maroon', 'orange', 'yellow']") > -1) and \
           (s_stdout.find("grid='false'") > -1) and \
           (s_stdout.find("legend='reverse'") > -1) and \
           (s_stdout.find("yunit='myunit'") > -1) and \
           (s_stdout.find("title='my title'") > -1) and \
           (s_stdout.find("figsizepx=['641', '481']") > -1) and \
           (s_stdout.find("ext='tiff'") > -1) and \
           (s_stdout.find("figbgcolor='cyan'") > -1)
        os.remove(s_opathfile)

