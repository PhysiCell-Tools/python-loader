###
# title: test_cli_2d.py
#
# language: python3
# author: Elmar Bucher
# date: 2023-02-25
# license: BSD 3-Clause
#
# description:
#   pytest unit test library for the pcdl library command line interface functions.
#   + https://docs.pytest.org/
#
# note:
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
s_path_2d = str(pathlib.Path(pcdl.__file__).parent.resolve()/'output_2d')
s_file_2d = 'output00000024.xml'
s_pathfile_2d = f'{s_path_2d}/{s_file_2d}'
print("s_path_2d", s_path_2d)
print("s_pathfile_2d", s_pathfile_2d)

# probe for  data
if not os.path.exists(s_path_2d):
    pcdl.install_data()

print(f"process: pcdl pyCLI functions from the command line...")


##############################
# metadata realted test code #
##############################

class TestPyCliVersion(object):
    ''' tests for one pcdl.pyCli function. '''

    # timestep:
    # + path (pathfile, path) ok
    # + verbose (true, _false_) nop

    def test_pcdl_get_version_timeseries(self):
        o_result = subprocess.run(['pcdl_get_version', s_path_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        assert o_result.returncode == 0

    def test_pcdl_get_version_timestep(self):
        o_result = subprocess.run(['pcdl_get_version', s_pathfile_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        assert o_result.returncode == 0


class TestPyCliUnitDict(object):
    ''' tests for one pcdl.pyCli function. '''

    # timestep:
    # + path (pathfile, path) ok
    # + microenv (true, _false_) ok
    # + settingxml (string, _none_, _false_) ok
    # + verbose (true, _false_) nop

    def test_pcdl_get_unit_dict_timeseries(self):
        o_result = subprocess.run(['pcdl_get_unit_dict', s_path_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_unit.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_unit_dict_timestep(self):
        o_result = subprocess.run(['pcdl_get_unit_dict', s_pathfile_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_unit.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_unit_dict_timestep_microenv(self):
        o_result = subprocess.run(['pcdl_get_unit_dict', s_pathfile_2d, '--microenv', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_unit.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_unit_dict_timestep_settingxmlfalse(self):
        o_result = subprocess.run(['pcdl_get_unit_dict', s_pathfile_2d, '--settingxml', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_unit.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_unit_dict_timestep_settingxmlnone(self):
        o_result = subprocess.run(['pcdl_get_unit_dict', s_pathfile_2d, '--settingxml', 'none'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_unit.csv')
        assert o_result.returncode == 0


##############################
# substarte relatd test code #
##############################

class TestPyCliSubstrateList(object):
    ''' tests for one pcdl.pyCli function. '''

    # timestep:
    # + path (pathfile, path) ok
    # + verbose (true, _false_) nop

    def test_pcdl_get_substrate_list_timeseries(self):
        o_result = subprocess.run(['pcdl_get_substrate_list', s_path_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        assert o_result.returncode == 0

    def test_pcdl_get_substrate_list_timestep(self):
        o_result = subprocess.run(['pcdl_get_substrate_list', s_pathfile_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        assert o_result.returncode == 0


class TestPyCliConcDfAttribute(object):
    ''' tests for one pcdl.pyCli function. '''

    # timeseries collapsed:
    # + path (str) nop
    # + verbose (true, _false_) nop
    # + values (int) ok
    # + drop (oxygen) ok
    # + keep (oxygen) ok
    # + allvalues (false _true_) ok

    def test_pcdl_get_conc_attribute_timeseries(self):
        o_result = subprocess.run(['pcdl_get_conc_attribute', s_path_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_conc_attribute_minmax.json')
        assert o_result.returncode == 0

    def test_pcdl_get_conc_attribute_timeseries_value(self):
        o_result = subprocess.run(['pcdl_get_conc_attribute', s_path_2d, '2'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_conc_attribute_minmax.json')
        assert o_result.returncode == 0

    def test_pcdl_get_conc_attribute_timeseries_drop(self):
        o_result = subprocess.run(['pcdl_get_conc_attribute', s_path_2d, '--drop', 'conc_type', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_conc_attribute_minmax.json')
        assert o_result.returncode == 0

    def test_pcdl_get_conc_attribute_timeseries_keep(self):
        o_result = subprocess.run(['pcdl_get_conc_attribute', s_path_2d, '--keep', 'conc_type', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_conc_attribute_minmax.json')
        assert o_result.returncode == 0

    def test_pcdl_get_conc_attribute_timeseries_allvalues(self):
        o_result = subprocess.run(['pcdl_get_conc_attribute', s_path_2d, '--allvalues', 'true'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_conc_attribute_all.json')
        assert o_result.returncode == 0


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
        o_result = subprocess.run(['pcdl_get_conc_df', s_path_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_conc.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_conc_df_timeseries_collapsed(self):
        o_result = subprocess.run(['pcdl_get_conc_df', s_path_2d, '--collapse', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_conc.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_conc_df_timeseries_value(self):
        o_result = subprocess.run(['pcdl_get_conc_df', s_path_2d, '2'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_conc.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_conc_df_timeseries_drop(self):
        o_result = subprocess.run(['pcdl_get_conc_df', s_path_2d, '--drop', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_conc.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_conc_df_timeseries_keep(self):
        o_result = subprocess.run(['pcdl_get_conc_df', s_path_2d, '--keep', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_conc.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_conc_df_timestep(self):
        o_result = subprocess.run(['pcdl_get_conc_df', s_pathfile_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_conc.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_conc_df_timestep_value(self):
        o_result = subprocess.run(['pcdl_get_conc_df', s_pathfile_2d, '2'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_conc.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_conc_df_timestep_drop(self):
        o_result = subprocess.run(['pcdl_get_conc_df', s_pathfile_2d, '--drop', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_conc.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_conc_df_timestep_keep(self):
        o_result = subprocess.run(['pcdl_get_conc_df', s_pathfile_2d, '--keep', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_conc.csv')
        assert o_result.returncode == 0


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
        o_result = subprocess.run([
            'pcdl_plot_contour', s_pathfile_2d, 'oxygen'
        ], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        shutil.rmtree(f'{s_path_2d}/conc_oxygen_z0.0/')
        assert o_result.returncode == 0

    def test_pcdl_plot_contour_set(self):
        o_result = subprocess.run([
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
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        shutil.rmtree(f'{s_path_2d}/conc_oxygen_z0.0/')
        assert o_result.returncode == 0


class TestPyCliConcVtk(object):
    ''' tests for one pcdl.pyCli function. '''

    # timestep and timeseries:
    # + path nop
    # + verbose (true, _false_) nop

    def test_pcdl_make_conc_vtk_timeseries_default(self):
        o_result = subprocess.run(['pcdl_make_conc_vtk', s_path_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_conc.vtr')
        assert o_result.returncode == 0

    def test_pcdl_make_conc_vtk_timestep_default(self):
        o_result = subprocess.run(['pcdl_make_conc_vtk', s_pathfile_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_conc.vtr')
        assert o_result.returncode == 0


################################
# cell agent realted test code #
################################

class TestPyCliCelltypeList(object):
    ''' tests for one pcdl.pyCli function. '''

    # timestep:
    # + path (pathfile, path) ok
    # + verbose (true, _false_) nop

    def test_pcdl_get_celltype_list_timeseries(self):
        o_result = subprocess.run(['pcdl_get_celltype_list', s_path_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        assert o_result.returncode == 0

    def test_pcdl_get_celltype_list_timestep(self):
        o_result = subprocess.run(['pcdl_get_celltype_list', s_pathfile_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        assert o_result.returncode == 0


class TestPyCliCellDfAttribute(object):
    ''' tests for one pcdl.pyCli function. '''

    # timeseries collapsed:
    # + path (str) nop
    # + customtype ([], _sample:bool_) ok
    # + microenv (true, _false_) ok
    # + physiboss (true, _false_)
    # + settingxml (string, _none_, _false_) ok
    # + verbose (true, _false_) nop
    # + values (int) ok
    # + drop (cell_type oxygen) ok
    # + keep (cell_type oxygen) ok
    # + allvalues (false _true_) ok

    def test_pcdl_get_cell_attribute_timeseries(self):
        o_result = subprocess.run(['pcdl_get_cell_attribute', s_path_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_attribute_minmax.json')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_attribute_timeseries_customtype(self):
        o_result = subprocess.run(['pcdl_get_cell_attribute', s_path_2d, '--custom_data_type', 'sample:bool'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_attribute_minmax.json')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_attribute_timeseries_microenv(self):
        o_result = subprocess.run(['pcdl_get_cell_attribute', s_path_2d, '--microenv', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_attribute_minmax.json')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_attribute_timeseries_physiboss(self):
        o_result = subprocess.run(['pcdl_get_cell_attribute', s_path_2d, '--physiboss', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_attribute_minmax.json')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_attribute_timeseries_settingxmlfalse(self):
        o_result = subprocess.run(['pcdl_get_cell_attribute', s_path_2d, '--settingxml', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_attribute_minmax.json')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_attribute_timeseries_settingxmlnone(self):
        o_result = subprocess.run(['pcdl_get_cell_attribute', s_path_2d, '--settingxml', 'none'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_attribute_minmax.json')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_attribute_timeseries_value(self):
        o_result = subprocess.run(['pcdl_get_cell_attribute', s_path_2d, '2'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_attribute_minmax.json')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_attribute_timeseries_drop(self):
        o_result = subprocess.run(['pcdl_get_cell_attribute', s_path_2d, '--drop', 'cell_type', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_attribute_minmax.json')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_attribute_timeseries_keep(self):
        o_result = subprocess.run(['pcdl_get_cell_attribute', s_path_2d, '--keep', 'cell_type', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_attribute_minmax.json')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_attribute_timeseries_allvalues(self):
        o_result = subprocess.run(['pcdl_get_cell_attribute', s_path_2d, '--allvalues', 'true'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_attribute_all.json')
        assert o_result.returncode == 0


class TestPyCliCellAttributeList(object):
    ''' tests for one  pcdl command line interface  function. '''

    # timeseries collapsed:
    # + path (str) nop
    # + microenv (true, _false_) ok
    # + physiboss (true, _false_)
    # + settingxml (string, _none_, _false_) ok
    # + verbose (true, _false_) nop

    def test_pcdl_get_cell_attribute_list_timeseries(self):
        o_result = subprocess.run(['pcdl_get_cell_attribute_list', s_path_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_attribute_list_timeseries_microenv(self):
        o_result = subprocess.run(['pcdl_get_cell_attribute_list', s_path_2d, '--microenv', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_attribute_list_timeseries_physiboss(self):
        o_result = subprocess.run(['pcdl_get_cell_attribute_list', s_path_2d, '--physiboss', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_attribute_list_timeseries_settingxmlfalse(self):
        o_result = subprocess.run(['pcdl_get_cell_attribute_list', s_path_2d, '--settingxml', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_attribute_list_timeseries_settingxmlnone(self):
        o_result = subprocess.run(['pcdl_get_cell_attribute_list', s_path_2d, '--settingxml', 'none'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        assert o_result.returncode == 0


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
        o_result = subprocess.run(['pcdl_get_cell_df', s_path_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_df_timeseries_collapsed(self):
        o_result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--collapse', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_cell.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_df_timeseries_microenv(self):
        o_result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--microenv', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_df_timeseries_physiboss(self):
        o_result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--physiboss', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_df_timeseries_settingxmlfalse(self):
        o_result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--settingxml', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_df_timeseries_settingxmlnone(self):
        o_result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--settingxml', 'none'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_df_timeseries_value(self):
        o_result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '2'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_df_timeseries_drop(self):
        o_result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--drop', 'cell_type', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_df_timeseries_keep(self):
        o_result = subprocess.run(['pcdl_get_cell_df', s_path_2d, '--keep', 'cell_type', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_df_timestep(self):
        o_result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_df_timestep_microenv(self):
        o_result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '--microenv', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_df_timestep_physiboss(self):
        o_result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '--physiboss', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_df_timestep_settingxmlfalse(self):
        o_result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '--settingxml', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_df_timestep_settingxmlnone(self):
        o_result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '--settingxml', 'none'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_df_timestep_value(self):
        o_result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '2'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_df_timestep_drop(self):
        o_result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '--drop', 'cell_type', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell.csv')
        assert o_result.returncode == 0

    def test_pcdl_get_cell_df_timestep_keep(self):
        o_result = subprocess.run(['pcdl_get_cell_df', s_pathfile_2d, '--keep', 'cell_type', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell.csv')
        assert o_result.returncode == 0


class TestPyCliAnndata(object):
    ''' tests for one pcdl.pyCli function. '''

    # timeseries
    # + collapse (true false) ok

    # timestep and timeseries:
    # + path nop
    # + customtype ([], _sample:bool_) ok
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
        o_result = subprocess.run(['pcdl_get_anndata', s_path_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timeseries_collapsed(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--collapse', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timeseries_customtype(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--custom_data_type', 'sample:bool'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timeseries_microenv(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--microenv', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timeseries_graph(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--graph', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timeseries_physiboss(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--physiboss', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timeseries_settingxmlfalse(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--settingxml', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timeseries_settingxmlnone(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--settingxml', 'none'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timeseries_value(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '2'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timeseries_drop(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--drop', 'cell_type', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timeseries_keep(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--keep', 'cell_type', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timeseries_scale(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_path_2d, '--scale', 'std'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_std.h5ad')
        assert o_result.returncode == 0

    # timestep
    def test_pcdl_get_anndata_timestep(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timestep_microenv(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '--microenv', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timestep_graph(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '--graph', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timestep_physiboss(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '--physiboss', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timestep_settingxmlfalse(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '--settingxml', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timestep_settingxmlnone(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '--settingxml', 'none'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timestep_value(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '2'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timestep_drop(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '--drop', 'cell_type', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timestep_keep(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '--keep', 'cell_type', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell_maxabs.h5ad')
        assert o_result.returncode == 0

    def test_pcdl_get_anndata_timestep_scale(self):
        o_result = subprocess.run(['pcdl_get_anndata', s_pathfile_2d, '--scale', 'std'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell_std.h5ad')
        assert o_result.returncode == 0


class TestPyCliGraphGml(object):
    ''' tests for one pcdl.pyCli function. '''

    # timestep and timeseries:
    # + path nop
    # + customtype ([], _sample:bool_) ok
    # + microenv (true, false) ok
    # + physiboss (true, _false_)
    # + settingxml (string, _none_, _false_) ok
    # + verbose (true, _false_) nop
    # + graph_type (neighbor, _attached_) ok
    # + edge_attribute (true, _false_) ok
    # + node_attribute (cell_type oxygen) ok

    def test_pcdl_make_graph_gml_timeseries_default(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, 'neighbor'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_neighbor.gml')
        assert o_result.returncode == 0


    def test_pcdl_make_graph_gml_timeseries_customtype_nodeattribute_one(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, 'neighbor', '--custom_data_type', 'sample:bool', '--node_attribute', 'sample'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_neighbor.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timeseries_microenv(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, 'neighbor', '--microenv', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_neighbor.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timeseries_physiboss(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, 'neighbor', '--physiboss', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_neighbor.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timeseries_settingxmlfalse_nodeattribute_one(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, 'neighbor', '--settingxml', 'false', '--node_attribute', 'default_fusion_rates'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_neighbor.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timeseries_settingxmlnone_nodeattribute_one(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, 'neighbor', '--settingxml', 'none', '--node_attribute', 'default_fusion_rates'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_neighbor.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timeseries_graph_type(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, 'attached'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_attached.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timeseries_edge_attribute(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, 'neighbor', '--edge_attribute', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_neighbor.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timeseries_nodeattribute_one(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, 'neighbor', '--node_attribute', 'cell_type'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_neighbor.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timeseries_nodeattribute_many(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_path_2d, 'neighbor', '--node_attribute', 'cell_type', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_neighbor.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timestep_default(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, 'neighbor'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_neighbor.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timestep_customtype_nodeattribute_one(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, 'neighbor', '--custom_data_type', 'sample:bool', '--node_attribute', 'sample'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_neighbor.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timestep_microenv(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, 'neighbor', '--microenv', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_neighbor.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timestep_physiboss(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, 'neighbor', '--physiboss', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_neighbor.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timestep_settingxmlfalse_nodeattribute_one(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, 'neighbor', '--settingxml', 'false', '--node_attribute', 'default_fusion_rates'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_neighbor.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timestep_settingxmlnone_nodeattribute_one(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, 'neighbor', '--settingxml', 'none', '--node_attribute', 'default_fusion_rates'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_neighbor.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timestep_graph_type(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, 'attached'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_attached.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timestep_edge_attribute(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, 'neighbor', '--edge_attribute', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_neighbor.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timestep_node_attribute_one(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, 'neighbor', '--node_attribute', 'cell_type'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_neighbor.gml')
        assert o_result.returncode == 0

    def test_pcdl_make_graph_gml_timestep_node_attribute_many(self):
        o_result = subprocess.run(['pcdl_make_graph_gml', s_pathfile_2d, 'neighbor', '--node_attribute', 'cell_type', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_neighbor.gml')
        assert o_result.returncode == 0


class TestPyCliPlotScatter(object):
    ''' tests for one pcdl.pyCli function. '''

    # time series and time steps.
    # + path nop
    # + customtype ([], _sample:bool_) ok
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
        o_result = subprocess.run([
            'pcdl_plot_scatter', s_pathfile_2d,
        ], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        shutil.rmtree(f'{s_path_2d}/cell_cell_type_z0.0/')
        assert o_result.returncode == 0

    def test_pcdl_plot_scatter_set(self):
        o_result = subprocess.run([
            'pcdl_plot_scatter', s_pathfile_2d, 'oxygen',
            '--custom_data_type', 'sample:bool',
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
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        shutil.rmtree(f'{s_path_2d}/cell_oxygen_z0.0/')
        assert o_result.returncode == 0


class TestPyCliCellVtk(object):
    ''' tests for one pcdl.pyCli function. '''

    # timestep and timeseries:
    # + path nop
    # + customtype ([], _sample:bool_) ok
    # + microenv (true, _false) ok
    # + physiboss (true, _false_)
    # + settingxml (string, _none_, _false_) ok
    # + verbose (true, _false_) nop
    # + attribute (cell_type oxygen, _empty_) ok

    def test_pcdl_make_cell_vtk_timeseries_default(self):
        o_result = subprocess.run(['pcdl_make_cell_vtk', s_path_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_cell.vtp')
        assert o_result.returncode == 0

    def test_pcdl_make_cell_vtk_timeseries_customtype_attribute_one(self):
        o_result = subprocess.run(['pcdl_make_cell_vtk', s_path_2d, 'sample', '--custom_data_type', 'sample:bool'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_cell.vtp')
        assert o_result.returncode == 0

    def test_pcdl_make_cell_vtk_timeseries_microenv(self):
        o_result = subprocess.run(['pcdl_make_cell_vtk', s_path_2d, '--microenv', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_cell.vtp')
        assert o_result.returncode == 0

    def test_pcdl_make_cell_vtk_timeseries_physiboss(self):
        o_result = subprocess.run(['pcdl_make_cell_vtk', s_path_2d, '--physiboss', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_cell.vtp')
        assert o_result.returncode == 0

    def test_pcdl_make_cell_vtk_timeseries_settingxmlfalse_attribute_one(self):
        o_result = subprocess.run(['pcdl_make_cell_vtk', s_path_2d, 'default_fusion_rates', '--settingxml', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_cell.vtp')
        assert o_result.returncode == 0

    def test_pcdl_make_cell_vtk_timeseries_settingxmlnone_attribute_one(self):
        o_result = subprocess.run(['pcdl_make_cell_vtk', s_path_2d, 'default_fusion_rates', '--settingxml', 'none'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_cell.vtp')
        assert o_result.returncode == 0

    def test_pcdl_make_cell_vtk_timeseries_attribute_many(self):
        o_result = subprocess.run(['pcdl_make_cell_vtk', s_path_2d, 'cell_type', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        for i_step in range(25):
            os.remove(f'{s_path_2d}/output000000{str(i_step).zfill(2)}_cell.vtp')
        assert o_result.returncode == 0

    def test_pcdl_make_cell_vtk_timestep_default(self):
        o_result = subprocess.run(['pcdl_make_cell_vtk', s_pathfile_2d], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell.vtp')
        assert o_result.returncode == 0

    def test_pcdl_make_cell_vtk_timestep_customtype_attribute_one(self):
        o_result = subprocess.run(['pcdl_make_cell_vtk', s_pathfile_2d, 'sample', '--custom_data_type', 'sample:bool'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell.vtp')
        assert o_result.returncode == 0

    def test_pcdl_make_cell_vtk_timestep_microenv(self):
        o_result = subprocess.run(['pcdl_make_cell_vtk', s_pathfile_2d, '--microenv', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell.vtp')
        assert o_result.returncode == 0

    def test_pcdl_make_cell_vtk_timestep_physiboss(self):
        o_result = subprocess.run(['pcdl_make_cell_vtk', s_pathfile_2d, '--physiboss', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell.vtp')
        assert o_result.returncode == 0

    def test_pcdl_make_cell_vtk_timestep_settingxmlfalse_attribute_one(self):
        o_result = subprocess.run(['pcdl_make_cell_vtk', s_pathfile_2d, 'default_fusion_rates', '--settingxml', 'false'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell.vtp')
        assert o_result.returncode == 0

    def test_pcdl_make_cell_vtk_timestep_settingxmlnone_attribute_one(self):
        o_result = subprocess.run(['pcdl_make_cell_vtk', s_pathfile_2d, 'default_fusion_rates', '--settingxml', 'none'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell.vtp')
        assert o_result.returncode == 0

    def test_pcdl_make_cell_vtk_timestep_attribute_many(self):
        o_result = subprocess.run(['pcdl_make_cell_vtk', s_pathfile_2d, 'cell_type', 'oxygen'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/output00000024_cell.vtp')
        assert o_result.returncode == 0


#######################################
# substrate and cell agenat test code #
#######################################

class TestPyCliPlotTimeSeries(object):
    ''' tests for one pcdl.pyCli function. '''

    # time series.
    # + path nop
    # + customtype ([], _sample:bool_) ok
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
        o_result = subprocess.run([
            'pcdl_plot_timeseries', s_path_2d, '-v', 'false'
        ], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_cell_total_count.jpeg')
        assert o_result.returncode == 0

    def test_pcdl_plot_timeseries_set(self):
        o_result = subprocess.run([
            'pcdl_plot_timeseries', s_path_2d, 'None', 'oxygen', 'entropy', '-v', 'false',
            '--custom_data_type', 'sample:bool',
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
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        os.remove(f'{s_path_2d}/timeseries_conc_total_oxygen_entropy.tiff')
        assert o_result.returncode == 0


###########################
# making movies test code #
###########################

class TestPyCliMakeGif(object):
    ''' tests for one pcdl.pyCli function. '''

    # time series
    # + path nop
    # + interface (default, 'tiff')

    def test_pcdl_make_gif_timeseries_default(self):
        o_path = subprocess.run(['pcdl_plot_scatter', s_path_2d], check=False, capture_output=True)
        print(f'o_path: {o_path}\n')
        s_path = f'{s_path_2d}/cell_cell_type_z0.0/'
        o_result = subprocess.run(['pcdl_make_gif', s_path], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        shutil.rmtree(s_path)
        assert o_result.returncode == 0

    def test_pcdl_make_gif_timeseries_interface(self):
        o_path = subprocess.run(['pcdl_plot_contour', s_path_2d, 'oxygen', '--ext', 'tiff'], check=False, capture_output=True)
        print(f'o_path: {o_path}\n')
        s_path = f'{s_path_2d}/conc_oxygen_z0.0/'
        o_result = subprocess.run(['pcdl_make_gif', s_path, 'tiff'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        shutil.rmtree(s_path)
        assert o_result.returncode == 0


class TestPyCliMakeMovie(object):
    ''' tests for one pcdl.pyCli function. '''

    # time series
    # + path nop
    # + interface (default, 'tiff')
    # + farme (default, 'tiff')

    def test_pcdl_make_movie_timeseries_default(self):
        o_path = subprocess.run(['pcdl_plot_scatter', s_path_2d], check=False, capture_output=True)
        print(f'o_path: {o_path}\n')
        s_path = f'{s_path_2d}/cell_cell_type_z0.0/'
        o_result = subprocess.run(['pcdl_make_movie', s_path], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        shutil.rmtree(s_path)
        assert o_result.returncode == 0

    def test_pcdl_make_movie_timeseries_interface(self):
        o_path = subprocess.run(['pcdl_plot_contour', s_path_2d, 'oxygen', '--ext', 'tiff'], check=False, capture_output=True)
        print(f'o_path: {o_path}\n')
        s_path = f'{s_path_2d}/conc_oxygen_z0.0/'
        o_result = subprocess.run(['pcdl_make_movie', s_path, 'tiff'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        shutil.rmtree(s_path)
        assert o_result.returncode == 0

    def test_pcdl_make_movie_timeseries_farme(self):
        o_path = subprocess.run(['pcdl_plot_contour', s_path_2d, 'oxygen', '--ext', 'jpeg'], check=False, capture_output=True)
        print(f'o_path: {o_path}\n')
        s_path = f'{s_path_2d}/conc_oxygen_z0.0/'
        o_result = subprocess.run(['pcdl_make_movie', s_path, '--framerate', '9'], check=False, capture_output=True)
        print(f'o_result: {o_result}\n')
        print(f'o_result.returncode: {o_result.returncode}\n')
        print(f'o_result.stdout: {o_result.stdout}\n')
        print(f'o_result.stderr: {o_result.stderr}\n')
        shutil.rmtree(s_path)
        assert o_result.returncode == 0
