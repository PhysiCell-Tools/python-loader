####
# title: test_anndata.py
#
# language: python3
# author: Elmar Bucher
# date: 2023-06-24
# license: BSD 3-Clause
#
# description:
#   pytest unit test library for the pcdl library TimeStep and TimeSeries class.
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
import pandas as pd
import pathlib
import pcdl
import pytest


# const
s_path_3d = str(pathlib.Path(pcdl.__file__).parent.resolve()/'data_timeseries_3d')
s_file_3d = 'output00000024.xml'
s_pathfile_3d = f'{s_path_3d}/{s_file_3d}'


# test data
if not os.path.exists(s_path_3d):
    pcdl.install_data()


# load physicell data time step
class TestTimeStep(object):
    ''' test for pcdl.TimeStep class. '''
    mcds = pcdl.TimeStep(s_pathfile_3d, verbose=False)

    ## get_anndata command ##
    @pytest.mark.filterwarnings("ignore:invalid value encountered in divide")
    def test_mcds_get_anndata(self, mcds=mcds):
        ann = mcds.get_anndata(states=1, drop=set(), keep=set(), scale='maxabs')
        assert (str(type(ann)) == "<class 'anndata._core.anndata.AnnData'>") and \
               (ann.X.shape == (20460, 102)) and \
               (ann.obs.shape == (20460, 6)) and \
               (ann.obsm['spatial'].shape == (20460, 3)) and \
               (ann.var.shape == (102, 0)) and \
               (len(ann.uns) == 0)


# load physicell data time series
class TestTimeSeries(object):
    ''' test for pcdl.TestSeries class. '''
    mcdsts = pcdl.TimeSeries(s_path_3d, verbose=False)

    ## get_anndata command ##
    @pytest.mark.filterwarnings("ignore:invalid value encountered in divide")
    def test_mcdsts_get_anndata_collapse_mcds(self, mcdsts=mcdsts):
        ann = mcdsts.get_anndata(states=1, drop=set(), keep=set(), scale='maxabs', collapse=True, keep_mcds=True)
        assert (len(mcdsts.l_mcds) == 25) and \
               (str(type(ann)) == "<class 'anndata._core.anndata.AnnData'>") and \
               (ann.X.shape == (481651, 102)) and \
               (ann.obs.shape == (481651, 7)) and \
               (ann.obsm['spatial'].shape == (481651, 3)) and \
               (ann.var.shape == (102, 0)) and \
               (len(ann.uns) == 0)

    @pytest.mark.filterwarnings("ignore:invalid value encountered in divide")
    def test_mcdsts_get_anndata_noncollapse_nonmcds(self, mcdsts=mcdsts):
        l_ann = mcdsts.get_anndata(states=1, drop=set(), keep=set(), scale='maxabs', collapse=False, keep_mcds=False)
        assert (len(mcdsts.l_mcds) == 0) and \
               (str(type(l_ann)) == "<class 'list'>") and \
               (len(l_ann) == 25) and \
               (all([str(type(ann)) == "<class 'anndata._core.anndata.AnnData'>" for ann in l_ann])) and \
               (l_ann[24].X.shape == (20460, 102)) and \
               (l_ann[24].obs.shape == (20460, 6)) and \
               (l_ann[24].obsm['spatial'].shape == (20460, 3)) and \
               (l_ann[24].var.shape == (102, 0)) and \
               (len(l_ann[24].uns) == 0)

