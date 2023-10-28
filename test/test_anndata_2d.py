####
# title: test_anndata_2d.py
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
s_path_2d = str(pathlib.Path(pcdl.__file__).parent.resolve()/'data_timeseries_2d')
s_file_2d = 'output00000024.xml'
s_pathfile_2d = f'{s_path_2d}/{s_file_2d}'


# test data
if not os.path.exists(s_path_2d):
    pcdl.install_data()


# helper function
class TestScaler(object):
    ''' test for pcdl.scaler function '''
    a_x = np.array([[ 1.,-1., 2., 0.],[ 2., 0., 0.,0.],[ 0., 1.,-1.,0.]])
    df_x = pd.DataFrame(a_x, columns=['a','b','c','d'])

    def test_scaler_none(self, df_x=df_x):
        df_scaled = pcdl.pyAnnData.scaler(df_x=df_x, scale=None)
        assert all(df_scaled == df_x)

    @pytest.mark.filterwarnings("ignore:invalid value encountered in divide")
    def test_scaler_minabs(self, df_x=df_x):
        df_scaled = pcdl.pyAnnData.scaler(df_x=df_x, scale='maxabs')
        assert (df_scaled.values.sum().round(3) == 2.0) and \
               (df_scaled.values.min().round(3) == -1.0) and \
               (df_scaled.values.max().round(3) == 1.0)

    @pytest.mark.filterwarnings("ignore:invalid value encountered in divide")
    def test_scaler_minmax(self, df_x=df_x):
        df_scaled = pcdl.pyAnnData.scaler(df_x=df_x, scale='minmax')
        assert (df_scaled.values.sum().round(3) == 4.333) and \
               (df_scaled.values.min().round(3) == 0.0) and \
               (df_scaled.values.max().round(3) == 1.0)

    @pytest.mark.filterwarnings("ignore:invalid value encountered in divide")
    def test_scaler_std(self, df_x=df_x):
        df_scaled = pcdl.pyAnnData.scaler(df_x=df_x, scale='std')
        assert (df_scaled.values.sum().round(3) == 0.0) and \
               (df_scaled.values.min().round(3) == -1.0) and \
               (df_scaled.values.max().round(3) == 1.091)


# load physicell data time step
class TestTimeStep(object):
    ''' test for pcdl.TimeStep class. '''
    mcds = pcdl.TimeStep(s_pathfile_2d, verbose=False)

    ## get_anndata command ##
    @pytest.mark.filterwarnings("ignore:invalid value encountered in divide")
    def test_mcds_get_anndata(self, mcds=mcds):
        ann = mcds.get_anndata(values=1, drop=set(), keep=set(), scale='maxabs')
        assert (str(type(ann)) == "<class 'anndata._core.anndata.AnnData'>") and \
               (ann.X.shape == (1099, 79)) and \
               (ann.obs.shape == (1099, 6)) and \
               (ann.obsm['spatial'].shape == (1099, 2)) and \
               (ann.var.shape == (79, 0)) and \
               (len(ann.uns) == 0)


# load physicell data time series
class TestTimeSeries(object):
    ''' test for pcdl.TestSeries class. '''
    mcdsts = pcdl.TimeSeries(s_path_2d, verbose=False)

    ## get_anndata command ##
    @pytest.mark.filterwarnings("ignore:invalid value encountered in divide")
    def test_mcdsts_get_anndata_collapse_mcds(self, mcdsts=mcdsts):
        ann = mcdsts.get_anndata(values=1, drop=set(), keep=set(), scale='maxabs', collapse=True, keep_mcds=True)
        assert (len(mcdsts.l_mcds) == 25) and \
               (str(type(ann)) == "<class 'anndata._core.anndata.AnnData'>") and \
               (ann.X.shape == (24758, 79)) and \
               (ann.obs.shape == (24758, 7)) and \
               (ann.obsm['spatial'].shape == (24758, 2)) and \
               (ann.var.shape == (79, 0)) and \
               (len(ann.uns) == 0)

    @pytest.mark.filterwarnings("ignore:invalid value encountered in divide")
    def test_mcdsts_get_anndata_noncollapse_nonmcds(self, mcdsts=mcdsts):
        l_ann = mcdsts.get_anndata(values=1, drop=set(), keep=set(), scale='maxabs', collapse=False, keep_mcds=False)
        assert (len(mcdsts.l_mcds) == 0) and \
               (str(type(l_ann)) == "<class 'list'>") and \
               (len(l_ann) == 25) and \
               (all([str(type(ann)) == "<class 'anndata._core.anndata.AnnData'>" for ann in l_ann])) and \
               (l_ann[24].X.shape == (1099, 79)) and \
               (l_ann[24].obs.shape == (1099, 6)) and \
               (l_ann[24].obsm['spatial'].shape == (1099, 2)) and \
               (l_ann[24].var.shape == (79, 0)) and \
               (len(l_ann[24].uns) == 0)

    ## get_annmcds_list command ##
    def test_mcdsts_get_annmcds_list(self, mcdsts=mcdsts):
        l_annmcds = mcdsts.get_annmcds_list()
        assert l_annmcds == mcdsts.l_annmcds

