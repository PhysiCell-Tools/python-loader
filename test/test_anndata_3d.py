####
# title: test_anndata_3d.py
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
s_path_3d = str(pathlib.Path(pcdl.__file__).parent.resolve()/'output_3d')
s_file_3d = 'output00000024.xml'
s_pathfile_3d = f'{s_path_3d}/{s_file_3d}'


# test data
if not os.path.exists(s_path_3d):
    pcdl.install_data()


###########
# 3D only #
###########

## load physicell data time step  ##
class TestPyAnndata3DTimeStep(object):
    ''' test for pcdl.TimeStep class. '''

    ## get_anndata command ##
    def test_mcds_get_anndata(self):
        mcds = pcdl.TimeStep(s_pathfile_3d, verbose=False)
        ann = mcds.get_anndata(values=1, drop=set(), keep=set(), scale='maxabs')
        assert(str(type(mcds)) == "<class 'pcdl.pyAnnData.TimeStep'>") and \
              (str(type(ann)) == "<class 'anndata._core.anndata.AnnData'>") and \
              (ann.X.shape[0] > 9) and \
              (ann.X.shape[1] == 105) and \
              (ann.obs.shape[0] > 9) and \
              (ann.obs.shape[1] == 7) and \
              (ann.obsm['spatial'].shape[0] > 9) and \
              (ann.obsm['spatial'].shape[1] == 3) and \
              (len(ann.obsp) == 2) and \
              (ann.var.shape == (105, 0)) and \
              (len(ann.uns) == 1)


##################
# test for speed #
##################
# BUE: test functions are mirrored test_anndata_2d.py

## load physicell data time series ##
class TestPyAnndata3DTimeSeries(object):
    ''' test for pcdl.TestSeries class. '''

    # get_anndata
    # get_annmcds_list {integrated}
    # value {1, _2_}
    # collaps {True, _False_}
    # keep_mcds {True, _False_}

    ## get_anndata command ##
    def test_mcdsts_get_anndata(self):
        mcdsts = pcdl.TimeSeries(s_path_3d, verbose=True)
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
              (ann.obsm['spatial'].shape[1] == 3) and \
              (len(ann.obsp) == 0) and \
              (ann.var.shape == (105, 0)) and \
              (len(ann.uns) == 0)

    def test_mcdsts_get_anndata_value(self):
        mcdsts = pcdl.TimeSeries(s_path_3d, verbose=True)
        ann = mcdsts.get_anndata(values=2, drop=set(), keep=set(), scale='maxabs', collapse=True, keep_mcds=True)
        l_annmcds = mcdsts.get_annmcds_list()
        assert(str(type(mcdsts)) == "<class 'pcdl.pyAnnData.TimeSeries'>") and \
              (len(mcdsts.l_mcds) == 25) and \
              (l_annmcds == mcdsts.l_annmcds) and \
              (mcdsts.l_annmcds is None) and \
              (str(type(ann)) == "<class 'anndata._core.anndata.AnnData'>") and \
              (ann.X.shape[0] > 9) and \
              (ann.X.shape[1] == 56) and \
              (ann.obs.shape[0] > 9) and \
              (ann.obs.shape[1] == 7) and \
              (ann.obsm['spatial'].shape[0] > 9) and \
              (ann.obsm['spatial'].shape[1] == 3) and \
              (len(ann.obsp) == 0) and \
              (ann.var.shape == (56, 0)) and \
              (len(ann.uns) == 0)

    def test_mcdsts_get_anndata_collapsefalse(self):
        mcdsts = pcdl.TimeSeries(s_path_3d, verbose=True)
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
              (mcdsts.l_annmcds[24].obsm['spatial'].shape[1] == 3) and \
              (len(mcdsts.l_annmcds[24].obsp) == 4) and \
              (mcdsts.l_annmcds[24].var.shape == (105, 0)) and \
              (len(mcdsts.l_annmcds[24].uns) == 2)

    def test_mcdsts_get_anndata_keepmcdsfalse(self):
        mcdsts = pcdl.TimeSeries(s_path_3d, verbose=True)
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
              (ann.obsm['spatial'].shape[1] == 3) and \
              (len(ann.obsp) == 0) and \
              (ann.var.shape == (105, 0)) and \
              (len(ann.uns) == 0)

