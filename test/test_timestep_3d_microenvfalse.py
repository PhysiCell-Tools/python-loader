#####
# title: test_snapshot_3d_microenvfalse.py
#
# language: python3
# author: Elmar Bucher
# date: 2022-10-15
# license: BSD 3-Clause
#
# description:
#   pytest unit test library for the pcdl library pyMCDS class.
#   + https://docs.pytest.org/
#
#   note:
#   assert actual == expected, message
#   == value equality
#   is reference equality
#   pytest.approx for real values
#####


# load library
import matplotlib.pyplot as plt
import os
import pathlib
import pcdl


# const
s_path_3d = str(pathlib.Path(pcdl.__file__).parent.resolve()/'data_timeseries_3d')
s_file_3d = 'output00000024.xml'
s_pathfile_3d = f'{s_path_3d}/{s_file_3d}'


# test data
if not os.path.exists(s_path_3d):
    pcdl.install_data()


# test function
class TestyMcdsMicroenvFalse3D(object):
    ''' test for pcdl.pyMCDS data loader with microenvironment, graph, and settingxml set to False '''
    mcds = pcdl.pyMCDS(xmlfile=s_file_3d, output_path=s_path_3d, custom_type={}, microenv=False, graph=False, settingxml=False, verbose=True)

    def test_pyMCDS(self, mcds=mcds):
        # load physicell data
        print(f"process: pcdl.pyMCDS(xmlfile={s_file_3d}, output_path={s_path_3d}, custom_type={{}}, microenv=False, graph=False, settingxml=False, verbose=True) ...")
        assert str(type(mcds)) == "<class 'pcdl.pyMCDS.pyMCDS'>"

    ## metadata related functions
    # nop

    ## mesh related functions
    def test_mcds_get_voxel_ijk_range(self, mcds=mcds):
        ltr_range = mcds.get_voxel_ijk_range()
        assert ltr_range == [(0, 10), (0, 10), (0, 10)]

    def test_mcds_get_mesh_mnp_range(self, mcds=mcds):
        ltr_range = mcds.get_mesh_mnp_range()
        assert ltr_range == [(-15, 285), (-10, 190), (-5, 95)]

    def test_mcds_get_xyz_range(self, mcds=mcds):
        ltr_range = mcds.get_xyz_range()
        assert ltr_range == [(-30, 300), (-20, 200), (-10, 100)]

    def test_mcds_get_voxel_ijk_axis(self, mcds=mcds):
        lar_axis = mcds.get_voxel_ijk_axis()
        assert (str(type(lar_axis)) == "<class 'list'>") and \
               (len(lar_axis) == 3) and \
               (str(type(lar_axis[0])) == "<class 'numpy.ndarray'>") and \
               (str(lar_axis[0].dtype) in {"int64", "int32"}) and \
               (lar_axis[0].shape == (11,)) and \
               (str(type(lar_axis[1])) == "<class 'numpy.ndarray'>") and \
               (str(lar_axis[1].dtype) in {"int64", "int32"}) and \
               (lar_axis[1].shape == (11,)) and \
               (str(type(lar_axis[2])) == "<class 'numpy.ndarray'>") and \
               (str(lar_axis[2].dtype) in {"int64", "int32"}) and \
               (lar_axis[2].shape == (11,))

    def test_mcds_get_mesh_mnp_axis(self, mcds=mcds):
        lar_axis = mcds.get_mesh_mnp_axis()
        assert (str(type(lar_axis)) == "<class 'list'>") and \
               (len(lar_axis) == 3) and \
               (str(type(lar_axis[0])) == "<class 'numpy.ndarray'>") and \
               (str(lar_axis[0].dtype) == "float64") and \
               (lar_axis[0].shape == (11,)) and \
               (str(type(lar_axis[1])) == "<class 'numpy.ndarray'>") and \
               (str(lar_axis[1].dtype) == "float64") and \
               (lar_axis[1].shape == (11,)) and \
               (str(type(lar_axis[2])) == "<class 'numpy.ndarray'>") and \
               (str(lar_axis[2].dtype) == "float64") and \
               (lar_axis[2].shape == (11,))

    def test_mcds_get_mesh_flat_false(self, mcds=mcds):
        aar_mesh = mcds.get_mesh(flat=False)
        assert (str(type(aar_mesh)) == "<class 'numpy.ndarray'>") and \
               (len(aar_mesh) == 3) and \
               (str(type(aar_mesh[0])) == "<class 'numpy.ndarray'>") and \
               (str(aar_mesh[0].dtype) == "float64") and \
               (aar_mesh[0].shape == (11, 11, 11)) and \
               (str(type(aar_mesh[1])) == "<class 'numpy.ndarray'>") and \
               (str(aar_mesh[1].dtype) == "float64") and \
               (aar_mesh[1].shape == (11, 11, 11)) and \
               (str(type(aar_mesh[2])) == "<class 'numpy.ndarray'>") and \
               (str(aar_mesh[2].dtype) == "float64") and \
               (aar_mesh[2].shape == (11, 11, 11))

    def test_mcds_get_mesh_flat_true(self, mcds=mcds):
        aar_mesh = mcds.get_mesh(flat=True)
        assert (str(type(aar_mesh)) == "<class 'numpy.ndarray'>") and \
               (len(aar_mesh) == 2) and \
               (str(type(aar_mesh[0])) == "<class 'numpy.ndarray'>") and \
               (str(aar_mesh[0].dtype) == "float64") and \
               (aar_mesh[0].shape == (11, 11)) and \
               (str(type(aar_mesh[1])) == "<class 'numpy.ndarray'>") and \
               (str(aar_mesh[1].dtype) == "float64") and \
               (aar_mesh[1].shape == (11, 11))

    def test_mcds_get_mesh_2d(self, mcds=mcds):
        aar_mesh_flat = mcds.get_mesh(flat=True)
        aar_mesh_2d = mcds.get_mesh_2D()
        assert (str(type(aar_mesh_2d)) == "<class 'numpy.ndarray'>") and \
               (len(aar_mesh_2d) == 2) and \
               (aar_mesh_2d[0] == aar_mesh_flat[0]).all() and \
               (aar_mesh_2d[1] == aar_mesh_flat[1]).all()

    def test_mcds_get_mesh_coordinate(self, mcds=mcds):
        # cube coordinates
        ar_m_cube, ar_n_cube, ar_p_cube = mcds.get_mesh(flat=False)
        er_m_cube = set(ar_m_cube.flatten())
        er_n_cube = set(ar_n_cube.flatten())
        er_p_cube = set(ar_p_cube.flatten())
        # linear coordinates
        aar_voxel = mcds.get_mesh_coordinate()
        assert (str(type(aar_voxel)) == "<class 'numpy.ndarray'>") and \
               (len(aar_voxel) == 3) and \
               (str(type(aar_voxel[0])) == "<class 'numpy.ndarray'>") and \
               (str(aar_voxel[0].dtype) == "float64") and \
               (set(aar_voxel[0]) == er_m_cube) and \
               (aar_voxel[0].shape == (1331,)) and \
               (str(type(aar_voxel[1])) == "<class 'numpy.ndarray'>") and \
               (str(aar_voxel[1].dtype) == "float64") and \
               (set(aar_voxel[1]) == er_n_cube) and \
               (aar_voxel[1].shape == (1331,)) and \
               (str(type(aar_voxel[2])) == "<class 'numpy.ndarray'>") and \
               (str(aar_voxel[2].dtype) == "float64") and \
               (set(aar_voxel[2]) == er_p_cube) and \
               (aar_voxel[2].shape == (1331,))

    def test_mcds_get_voxel_volume(self, mcds=mcds):
        r_volume = mcds.get_voxel_volume()
        assert r_volume == 6000.0

    def test_mcds_get_mesh_spacing(self, mcds=mcds):
        lr_spacing = mcds.get_mesh_spacing()
        assert lr_spacing == [30.0, 20.0, 10.0]

    def test_mcds_get_voxel_spacing(self, mcds=mcds):
        lr_spacing = mcds.get_voxel_spacing()
        assert lr_spacing == [30.0, 20.0, 10.0]

    def test_mcds_is_in_mesh(self, mcds=mcds):
        assert mcds.is_in_mesh(x=42, y=42, z=42, halt=False) and \
               not mcds.is_in_mesh(x=-42, y=-42, z=-42, halt=False)

    def test_mcds_get_voxel_ijk(self, mcds=mcds):
        li_voxel_0 = mcds.get_voxel_ijk(x=0, y=0, z=0, is_in_mesh=True)
        li_voxel_1 = mcds.get_voxel_ijk(x=15, y=10, z=5, is_in_mesh=True)
        li_voxel_2 = mcds.get_voxel_ijk(x=30, y=20, z=10, is_in_mesh=True)
        li_voxel_3 = mcds.get_voxel_ijk(x=333, y=222, z=111, is_in_mesh=True)
        assert (li_voxel_0 == [0, 0, 0]) and \
               (li_voxel_1 == [1, 1, 1]) and \
               (li_voxel_2 == [2, 2, 2]) and \
               (li_voxel_3 is None)

    ## micro environment related functions
    def test_mcds_get_substrate_names(self, mcds=mcds):
        ls_substrate = mcds.get_substrate_names()
        assert (str(type(ls_substrate)) == "<class 'list'>") and \
               (len(ls_substrate) == 0)

    def test_mcds_get_substrate_dict(self, mcds=mcds):
        ds_substrate = mcds.get_substrate_dict()
        assert (str(type(ds_substrate)) == "<class 'dict'>") and \
               (len(ds_substrate) == 0)  # settingxml=False

    ## cell related functions
    def test_mcds_get_cell_variables(self, mcds=mcds):
        ls_variable = mcds.get_cell_variables()
        assert (str(type(ls_variable)) == "<class 'list'>") and \
               (len(ls_variable) == 97) and \
               (ls_variable[0] == 'ID')

    def test_mcds_get_celltype_dict(self, mcds=mcds):
        ds_celltype = mcds.get_celltype_dict()
        assert (str(type(ds_celltype)) == "<class 'dict'>") and \
               (len(ds_celltype) == 0)  # settingxml=False

    def test_mcds_get_cell_df(self, mcds=mcds):
        df_cell = mcds.get_cell_df(values=0, drop=set(), keep=set())
        assert (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_cell.shape == (20460, 112))

    def test_mcds_get_cell_df_features(self, mcds=mcds):
        df_cell = mcds.get_cell_df(values=2, drop=set(), keep=set())
        assert (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_cell.shape == (20460, 31))

    def test_mcds_get_cell_df_keep(self, mcds=mcds):
        df_cell = mcds.get_cell_df(values=0, drop=set(), keep={'total_volume'})
        assert (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_cell.shape == (20460, 13))

    def test_mcds_get_cell_df_at(self, mcds=mcds):
        df_cell = mcds.get_cell_df_at(x=0, y=0, z=0, values=1, drop=set(), keep=set())
        assert (str(type(df_cell)) == "<class 'pandas.core.frame.DataFrame'>") and \
               (df_cell.shape == (5, 112))

    def test_mcds_plot_scatter_cat(self, mcds=mcds):
        fig = mcds.plot_scatter(
            focus='cell_type',  # case categorical
            z_slice = -3.333,   # test if
            z_axis = None,  # test if categorical
            #alpha = 1,  # matplotlib 
            #cmap = 'viridis',  # matplotlib
            #title = None, # matplotlib
            #grid = True,  # matplotlib
            #legend_loc = 'lower left',  # matplotlib
            xlim = None,  # test if
            ylim = None,  # test if
            xyequal = True,  # test if
            s = None,  # test if
            figsize = (6.4, 4.8),  # test if case
            ax = None,  # generate matplotlib figure
        )
        assert (str(type(fig)) == "<class 'matplotlib.figure.Figure'>")

    def test_mcds_plot_scatter_cat_cmap(self, mcds=mcds):
        fig, ax = plt.subplots()
        mcds.plot_scatter(
            focus='cell_type',  # case categorical
            z_slice = 0,  # jump over if
            z_axis = None,  # test if categorical
            #alpha = 1,  # matplotlib
            cmap = {'0': 'maroon'},  # matplotlib and label '0' because settingxml=False.
            title ='scatter cat',  # matplotlib
            #grid = True,  # matplotlib
            #legend_loc = 'lower left',  # matplotlib
            xlim = None,  # test if
            ylim = None,  # test if
            xyequal = True,  # test if
            s = None,  # test if
            #figsize = None,  # test if ax case
            ax = ax,  # use axis from existing matplotlib figure
        )
        assert (str(type(fig)) == "<class 'matplotlib.figure.Figure'>")

    def test_mcds_plot_scatter_num(self, mcds=mcds):
        fig = mcds.plot_scatter(
            focus='pressure',  # case numeric
            z_slice = -3.333,   # test if
            z_axis = None,  # test if numeric
            #alpha = 1,  # matplotlib
            #cmap = 'viridis',  # matplotlib
            #title = None, # matplotlib
            #grid = True,  # matplotlib
            #legend_loc = 'lower left',  # matplotlib
            xlim = None,  # test if
            ylim = None,  # test if
            xyequal = True,  # test if
            s = None,  # test if
            figsize = None,  # else case
            ax = None,  # generate matplotlib figure
        )
        assert (str(type(fig)) == "<class 'matplotlib.figure.Figure'>")

    ## unit related functions
    def test_mcds_get_unit_se(self, mcds=mcds):
        se_unit = mcds.get_unit_se()
        assert (str(type(se_unit)) == "<class 'pandas.core.series.Series'>") and \
               (se_unit.shape == (99,))

