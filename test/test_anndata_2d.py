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


# const
s_path_2d = str(pathlib.Path(pcdl.__file__).parent.resolve()/'output_2d')
s_file_2d = 'output00000024.xml'
s_pathfile_2d = f'{s_path_2d}/{s_file_2d}'


# test data
if not os.path.exists(s_path_2d):
    pcdl.install_data()


## helper function ##
class TestPyAnndataScaler(object):
    ''' test for pcdl.scaler function '''
    a_x = np.array([[ 1.,-1., 2., 0.],[ 2., 0., 0.,0.],[ 0., 1.,-1.,0.]])
    df_x = pd.DataFrame(a_x, columns=['a','b','c','d'])

    def test_scaler_none(self, df_x=df_x):
        df_scaled = pcdl.pyAnnData.scaler(df_x=df_x, scale=None)
        assert(str(type(df_scaled)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (all(df_scaled == df_x))

    def test_scaler_minabs(self, df_x=df_x):
        df_scaled = pcdl.pyAnnData.scaler(df_x=df_x, scale='maxabs')
        assert(str(type(df_scaled)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_scaled.values.sum().round(3) == 2.0) and \
              (df_scaled.values.min().round(3) == -1.0) and \
              (df_scaled.values.max().round(3) == 1.0)

    def test_scaler_minmax(self, df_x=df_x):
        df_scaled = pcdl.pyAnnData.scaler(df_x=df_x, scale='minmax')
        assert(str(type(df_scaled)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_scaled.values.sum().round(3) == 4.333) and \
              (df_scaled.values.min().round(3) == 0.0) and \
              (df_scaled.values.max().round(3) == 1.0)

    def test_scaler_std(self, df_x=df_x):
        df_scaled = pcdl.pyAnnData.scaler(df_x=df_x, scale='std')
        assert(str(type(df_scaled)) == "<class 'pandas.core.frame.DataFrame'>") and \
              (df_scaled.values.sum().round(3) == 0.0) and \
              (df_scaled.values.min().round(3) == -1.0) and \
              (df_scaled.values.max().round(3) == 1.091)


## load physicell data time step  ##
class TestPyAnndataTimeStep(object):
    ''' test for pcdl.TimeStep class. '''

    ## get_anndata command ##
    def test_mcds_get_anndata(self):
        mcds = pcdl.TimeStep(s_pathfile_2d, verbose=False)
        ann = mcds.get_anndata(values=1, drop=set(), keep=set(), scale='maxabs')
        assert(str(type(mcds)) == "<class 'pcdl.pyAnnData.TimeStep'>") and \
              (str(type(ann)) == "<class 'anndata._core.anndata.AnnData'>") and \
              (ann.X.shape[0] > 9) and \
              (ann.X.shape[1] == 105) and \
              (ann.obs.shape[0] > 9) and \
              (ann.obs.shape[1] == 7) and \
              (ann.obsm['spatial'].shape[0] > 9) and \
              (ann.obsm['spatial'].shape[1] == 2) and \
              (len(ann.obsp) == 2) and \
              (ann.var.shape == (105, 0)) and \
              (len(ann.uns) == 1)


## load physicell data time series ##
class TestPyAnndataTimeSeries(object):
    ''' test for pcdl.TestSeries class. '''

    # get_anndata
    # get_annmcds_list {integrated}
    # value {1, _2_}
    # collaps {True, _False_}
    # keep_mcds {True, _False_}

    ## get_anndata command ##
    def test_mcdsts_get_anndata(self):
        mcdsts = pcdl.TimeSeries(s_path_2d, verbose=True)
        ann = mcdsts.get_anndata(values=1, drop=set(), keep=set(), scale='maxabs', collapse=True, keep_mcds=True)
        l_annmcds = mcdsts.get_annmcds_list()
        assert(str(type(mcdsts)) == "<class 'pcdl.pyAnnData.TimeSeries'>") and \
              (len(mcdsts.l_mcds) == 25) and \
              (l_annmcds == mcdsts.l_annmcds) and \
              (mcdsts.l_annmcds is None) and \
              (str(type(ann)) == "<class 'anndata._core.anndata.AnnData'>") and \
              (ann.X.shape[0] > 9) and \
              (ann.X.shape[1] == 105) and \
              (ann.obs.shape[0] > 9) and \
              (ann.obs.shape[1] == 8) and \
              (ann.obsm['spatial'].shape[0] > 9) and \
              (ann.obsm['spatial'].shape[1] == 2) and \
              (len(ann.obsp) == 0) and \
              (ann.var.shape == (105, 0)) and \
              (len(ann.uns) == 0)

    def test_mcdsts_get_anndata_value(self):
        mcdsts = pcdl.TimeSeries(s_path_2d, verbose=True)
        ann = mcdsts.get_anndata(values=2, drop=set(), keep=set(), scale='maxabs', collapse=True, keep_mcds=True)
        l_annmcds = mcdsts.get_annmcds_list()
        assert(str(type(mcdsts)) == "<class 'pcdl.pyAnnData.TimeSeries'>") and \
              (len(mcdsts.l_mcds) == 25) and \
              (l_annmcds == mcdsts.l_annmcds) and \
              (mcdsts.l_annmcds is None) and \
              (str(type(ann)) == "<class 'anndata._core.anndata.AnnData'>") and \
              (ann.X.shape[0] > 9) and \
              (ann.X.shape[1] == 50) and \
              (ann.obs.shape[0] > 9) and \
              (ann.obs.shape[1] == 7) and \
              (ann.obsm['spatial'].shape[0] > 9) and \
              (ann.obsm['spatial'].shape[1] == 2) and \
              (len(ann.obsp) == 0) and \
              (ann.var.shape == (50, 0)) and \
              (len(ann.uns) == 0)

    def test_mcdsts_get_anndata_collapsefalse(self):
        mcdsts = pcdl.TimeSeries(s_path_2d, verbose=True)
        ann = mcdsts.get_anndata(values=1, drop=set(), keep=set(), scale='maxabs', collapse=False, keep_mcds=True)
        l_annmcds = mcdsts.get_annmcds_list()
        assert(str(type(mcdsts)) == "<class 'pcdl.pyAnnData.TimeSeries'>") and \
              (len(mcdsts.l_mcds) == 25) and \
              (l_annmcds == mcdsts.l_annmcds) and \
              (str(type(mcdsts.l_annmcds)) == "<class 'list'>") and \
              (len(mcdsts.l_annmcds) == 25) and \
              (all([str(type(ann)) == "<class 'anndata._core.anndata.AnnData'>" for ann in mcdsts.l_annmcds])) and \
              (mcdsts.l_annmcds[24].X.shape[0] > 9) and \
              (mcdsts.l_annmcds[24].X.shape[1] == 105) and \
              (mcdsts.l_annmcds[24].obs.shape[0] > 9) and \
              (mcdsts.l_annmcds[24].obs.shape[1] == 7) and \
              (mcdsts.l_annmcds[24].obsm['spatial'].shape[0] > 9) and \
              (mcdsts.l_annmcds[24].obsm['spatial'].shape[1] == 2) and \
              (len(mcdsts.l_annmcds[24].obsp) == 4) and \
              (mcdsts.l_annmcds[24].var.shape == (105, 0)) and \
              (len(mcdsts.l_annmcds[24].uns) == 2)

    def test_mcdsts_get_anndata_keepmcdsfalse(self):
        mcdsts = pcdl.TimeSeries(s_path_2d, verbose=True)
        ann = mcdsts.get_anndata(values=1, drop=set(), keep=set(), scale='maxabs', collapse=True, keep_mcds=False)
        l_annmcds = mcdsts.get_annmcds_list()
        assert(str(type(mcdsts)) == "<class 'pcdl.pyAnnData.TimeSeries'>") and \
              (len(mcdsts.l_mcds) == 0) and \
              (l_annmcds == mcdsts.l_annmcds) and \
              (mcdsts.l_annmcds is None) and \
              (str(type(ann)) == "<class 'anndata._core.anndata.AnnData'>") and \
              (ann.X.shape[0] > 9) and \
              (ann.X.shape[1] == 105) and \
              (ann.obs.shape[0] > 9) and \
              (ann.obs.shape[1] == 8) and \
              (ann.obsm['spatial'].shape[0] > 9) and \
              (ann.obsm['spatial'].shape[1] == 2) and \
              (len(ann.obsp) == 0) and \
              (ann.var.shape == (105, 0)) and \
              (len(ann.uns) == 0)

