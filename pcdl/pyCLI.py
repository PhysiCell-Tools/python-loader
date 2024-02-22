###
# title: pyCLI.py
#
# language: python3
# date: 2024-02-21
# license: BSD-3-Clause
# author: Elmar Bucher
#
# description:
#     pyCLI.py provides command line interface commands for appropriate pcdl functions.
#     all clis mirror the related python function interface as close as possible.
#     i like to thank Miguel Ponce-de-Leon for making me aware of the
#     entry point implementation technic which makes all of this possible.
#
# + https://hatch.pypa.io/latest/config/metadata/#entry-points
# + https://setuptools.pypa.io/en/latest/userguide/entry_point.html
####


# library
import argparse
import json
import numpy as np
import os
import pcdl
from scipy import stats
import sys

# functions
def entropy(pk):
    return stats.entropy(pk=pk, qk=None, base=2, nan_policy='omit', axis=0)[0]


# command line functions aplphabetically ordered.
def get_cell_df():
    # argv
    parser = argparse.ArgumentParser(
        prog = 'pcdl_get_cell_df',
        description = 'this function extracts dataframes with a cell centric view of the simulation and saves them as csv files.',
        epilog = 'homepage: https://github.com/elmbeech/physicelldataloader',
    )

    # TimeSeries path
    parser.add_argument(
        'path',
        nargs = '?',
        default = '.',
        help = 'path to the PhysiCell output directory or a outputnnnnnnnn.xml file. default is . .'
    )
    # TimeSeries output_path '.'
    # TimeSeries microenv
    parser.add_argument(
        '--microenv',
        nargs = '?',
        default = 'true',
        help = 'should the microenvironment be extracted? setting microenv to False will use less memory and speed up processing, similar to the original pyMCDS_cells.py script. default is True.'
    )
    # TimeSeries graph False
    # TimeSeries settingxml
    parser.add_argument(
        '--settingxml',
        nargs = '?',
        default = 'PhysiCell_settings.xml',
        help = 'from which settings.xml should the substrate and cell type ID label mapping be extracted? set to None or False if the xml file is missing! default is PhysiCell_settings.xml.',
    )
    # TimeSeries verbose
    parser.add_argument(
        '-v', '--verbose',
        nargs = '?',
        default = 'true',
        help = 'setting verbose to False for less text output, while processing. default is True.'
    )

    # get_cell_df values
    parser.add_argument(
        'values',
        nargs = '?',
        default = 1,
        type = int,
        help = 'minimal number of values a variable has to have in any of the mcds time steps to be outputted. variables that have only 1 state carry no information. None is a state too. default is 1.'
    )
    # get_cell_df drop
    parser.add_argument(
        '--drop',
        nargs = '?',
        default = '',
        help = "set of column labels to be dropped for the dataframe. don't worry: essential columns like ID, coordinates and time will never be dropped. Attention: when the keep parameter is given, then the drop parameter has to be an empty string! default is an empty string."
    )
    # get_cell_df keep
    parser.add_argument(
        '--keep',
        nargs = '?',
        default = '',
        help = "set of column labels to be kept in the dataframe. set values=1 to be sure that all variables are kept. don't worry: essential columns like ID, coordinates and time will always be kept. default is an empty string."
    )

    # parse arguments
    args = parser.parse_args()
    print(args)

    # process arguments
    s_pathfile = args.path
    if not s_pathfile.endswith('.xml'):
        s_pathfile = s_pathfile + '/initial.xml'
    if not os.path.exists(s_pathfile):
        sys.exit(f"Error @ pcdl_plot_timeseries : {args.path} path does not look like a outputnnnnnnnn.xml file or physicell output directory ({args.path}/initial.xml is missing).")
    if os.path.isfile(args.path):
        mcds = pcdl.pyMCDS(
            xmlfile = args.path,
            output_path = '.',
            microenv = False if args.microenv.lower().startswith('f') else True,
            graph = False,
            settingxml = args.settingxml,
            verbose = False if args.verbose.lower().startswith('f') else True
        )
        l_mcds = [mcds]
        ls_opathfile = [args.path.replace('.xml','_cell.csv')]
    else:
        mcdsts = pcdl.pyMCDSts(
            output_path = args.path,
            #custom_type,
            load = True,
            microenv = False if args.microenv.lower().startswith('f') else True,
            graph = False,
            settingxml = args.settingxml,
            verbose = False if args.verbose.lower().startswith('f') else True,
        )
        l_mcds = mcdsts.l_mcds
        ls_opathfile = [s_xmlfile.replace('.xml','_cell.csv') for s_xmlfile in mcdsts.get_xmlfile_list()]
    # run
    for i, mcds in enumerate(l_mcds):
        df_cell = mcds.get_cell_df(
            values = args.values,
            drop = set(args.drop.split()),
            keep = set(args.keep.split()),
        )
        df_cell.to_csv(ls_opathfile[i])
    # going home
    return ls_opathfile


def get_cell_df_features():
    # argv
    parser = argparse.ArgumentParser(
        prog = 'pcdl_get_cell_df_features',
        description = 'function to detect informative variables in a time series. this function detects even variables which have less than the minimal state count in each time step, but different values from time step to time step. the output is a json file with an entry of all non-coordinate column names that at least in one of the time steps or in between time steps, reach the given minimal value count. key is the column name, mapped is a list of all values (bool, str, and, if allvalues is True, int and float) or a list with minimum and maximum values (int, float).',
        epilog = 'homepage: https://github.com/elmbeech/physicelldataloader',
    )

    # TimeSeries path
    parser.add_argument(
        'path',
        nargs = '?',
        default = '.',
        help = 'path to the PhysiCell output directory. default is . .',
    )
    # TimeSeries output_path '.'
    # TimeSeries microenv
    parser.add_argument(
        '--microenv',
        nargs = '?',
        default = 'true',
        help = 'should the microenvironment be extracted? setting microenv to False will use less memory and speed up processing, similar to the original pyMCDS_cells.py script. default is True.',
    )
    # TimeSeries graph False
    # TimeSeries settingxml
    parser.add_argument(
        '--settingxml',
        nargs = '?',
        default = 'PhysiCell_settings.xml',
        help = 'from which settings.xml should the substrate and cell type ID label mapping be extracted? set to None or False if the xml file is missing! default is PhysiCell_settings.xml.',
    )
    # TimeSeries verbose
    parser.add_argument(
        '-v', '--verbose',
        nargs = '?',
        default = 'true',
        help = 'setting verbose to False for less text output, while processing. default is True.',
    )

    # get_cell_df_features values
    parser.add_argument(
        'values',
        nargs = '?',
        default = 1,
        type = int,
        help = 'minimal number of values a variable has to have in any of the mcds time steps to be outputted. variables that have only 1 state carry no information. None is a state too. default is 1.',
    )
    # get_cell_df_features drop
    parser.add_argument(
        '--drop',
        nargs = '?',
        default = '',
        help = "set of column labels to be dropped for the dataframe. don't worry: essential columns like ID, coordinates and time will never be dropped. Attention: when the keep parameter is given, then the drop parameter has to be an empty string! default is an empty string.",
    )
    # get_cell_df_features keep
    parser.add_argument(
        '--keep',
        nargs = '?',
        default = '',
        help = "set of column labels to be kept in the dataframe. set values=1 to be sure that all variables are kept. don't worry: essential columns like ID, coordinates and time will always be kept. default is an empty string.",
    )
    # get_cell_df_features allvalues
    parser.add_argument(
        'allvalues',
        nargs = '?',
        default = 'false',
        help = 'for numeric data, should only the min and max values or all values be returned? default is false.',
    )

    # parse arguments
    args = parser.parse_args()
    print(args)

    # process arguments
    if not os.path.exists(args.path + '/initial.xml'):
        sys.exit(f"Error @ pcdl_plot_timeseries : {args.path} path does not look like a physicell output directory ({args.path}/initial.xml is missing).")
    mcdsts = pcdl.pyMCDSts(
        output_path = args.path,
        #custom_type,
        load = True,
        microenv = False if args.microenv.lower().startswith('f') else True,
        graph = False,
        settingxml = args.settingxml,
        verbose = False if args.verbose.lower().startswith('f') else True,
    )
    # run
    s_values = 'minmax'
    b_allvalues = True if args.allvalues.lower().startswith('t') else False
    if b_allvalues:
        s_values = 'all'
    dl_variable = mcdsts.get_cell_df_features(
        values = args.values,
        drop = set(args.drop.split()),
        keep = set(args.keep.split()),
        allvalues = b_allvalues,
    )
    s_opathfile = f'{args.path}/cell_feature_{s_values}.json'
    json.dump(dl_variable, open(s_opathfile, 'w'), sort_keys=True)
    # going home
    return s_opathfile


def get_conc_df():
    # argv
    parser = argparse.ArgumentParser(
        prog = 'pcdl_get_conc_df',
        description = 'this function extracts dataframes with concentration values for all chemical species in all voxels and saves them as csv files. additionally, this dataframe lists voxel and mesh center coordinates.',
        epilog = 'homepage: https://github.com/elmbeech/physicelldataloader',
    )

    # TimeSeries path
    parser.add_argument(
        'path',
        nargs = '?',
        default = '.',
        help = 'path to the PhysiCell output directory or a outputnnnnnnnn.xml file. default is . .',
    )
    # TimeSeries output_path '.'
    # TimeSeries microenv True
    # TimeSeries graph False
    # TimeSeries settingxml
    parser.add_argument(
        '--settingxml',
        nargs = '?',
        default = 'PhysiCell_settings.xml',
        help = 'from which settings.xml should the substrate and cell type ID label mapping be extracted? set to None or False if the xml file is missing! default is PhysiCell_settings.xml.',
    )
    # TimeSeries verbose
    parser.add_argument(
        '-v', '--verbose',
        nargs = '?',
        default = 'true',
        help = 'setting verbose to False for less text output, while processing. default is True.',
    )

    # get_conc_df values
    parser.add_argument(
        'values',
        nargs = '?',
        default = 1,
        type = int,
        help = 'minimal number of values a variable has to have in any of the mcds time steps to be outputted. variables that have only 1 state carry no information. None is a state too. default is 1.',
    )
    # get_conc_df drop
    parser.add_argument(
        '--drop',
        nargs = '?',
        default = '',
        help = "set of column labels to be dropped for the dataframe. don't worry: essential columns like ID, coordinates and time will never be dropped. Attention: when the keep parameter is given, then the drop parameter has to be an empty string! default is an empty string.",
    )
    # get_conc_df keep
    parser.add_argument(
        '--keep',
        nargs = '?',
        default = '',
        help = "set of column labels to be kept in the dataframe. set values=1 to be sure that all variables are kept. don't worry: essential columns like ID, coordinates and time will always be kept. default is an empty string.",
    )

    # parse arguments
    args = parser.parse_args()
    print(args)

    # process arguments
    s_pathfile = args.path
    if not s_pathfile.endswith('.xml'):
        s_pathfile = s_pathfile + '/initial.xml'
    if not os.path.exists(s_pathfile):
        sys.exit(f"Error @ pcdl_plot_timeseries : {args.path} path does not look like a outputnnnnnnnn.xml file or physicell output directory ({args.path}/initial.xml is missing).")
    if os.path.isfile(args.path):
        mcds = pcdl.pyMCDS(
            xmlfile = args.path,
            output_path = '.',
            microenv = True,
            graph = False,
            settingxml = args.settingxml,
            verbose = False if args.verbose.lower().startswith('f') else True
        )
        l_mcds = [mcds]
        ls_opathfile = [args.path.replace('.xml','_conc.csv')]
    else:
        mcdsts = pcdl.pyMCDSts(
            output_path = args.path,
            #custom_type,
            load = True,
            microenv = True,
            graph = False,
            settingxml = args.settingxml,
            verbose = False if args.verbose.lower().startswith('f') else True,
        )
        l_mcds = mcdsts.l_mcds
        ls_opathfile = [s_xmlfile.replace('.xml','_conc.csv') for s_xmlfile in mcdsts.get_xmlfile_list()]
    # run
    for i, mcds in enumerate(l_mcds):
        df_conc = mcds.get_conc_df(
            values = args.values,
            drop = set(args.drop.split()),
            keep = set(args.keep.split()),
        )
        df_conc.to_csv(ls_opathfile[i], index=False)
    # going home
    return ls_opathfile


def get_conc_df_features():
    pass


def get_graph_gml():
    #parser.add_argument(
    #    'graph',
    #    nargs = '?',
    #    default = 'true',
    #    help = 'should the graphs be extracted? setting graph to False will use less memory and speed up processing.'
    #)
    pass


def get_version():
    # argv
    parser = argparse.ArgumentParser(
        prog = 'pcdl_get_version',
        description = 'this function is extracting PhysiCell and MultiCellDS version from the dataset and the installed pcdl module version.',
        epilog = 'homepage: https://github.com/elmbeech/physicelldataloader',
    )
    # TimeSeries path
    parser.add_argument(
        'path',
        nargs = '?',
        default = '.',
        help = 'path to the PhysiCell output directory or a outputnnnnnnnn.xml file. default is . .',
    )
    # TimeSeries output_path '.'
    # TimeSeries microenv False
    # TimeSeries graph False
    # TimeSeries settingxml None
    # TimeSeries verbose
    parser.add_argument(
        '-v', '--verbose',
        nargs = '?',
        default = 'false',
        help = 'setting verbose to True for more text output, while processing. default is False.',
    )

    # parse arguments
    args = parser.parse_args()
    print(args)

    # process arguments
    # run
    s_pathfile = args.path
    if not s_pathfile.endswith('.xml'):
        s_pathfile = s_pathfile + '/initial.xml'
    if not os.path.exists(s_pathfile):
        sys.exit(f"Error @ pcdl_plot_timeseries : {args.path} path does not look like a outputnnnnnnnn.xml file or physicell output directory ({args.path}/initial.xml is missing).")
    mcds = pcdl.pyMCDS(
        xmlfile = s_pathfile,
        output_path = '.',
        microenv = False,
        graph = False,
        settingxml = None,
        verbose = True if args.verbose.lower().startswith('t') else False
    )
    s_version = f'version:\n{mcds.get_physicell_version()}\n{mcds.get_multicellds_version()}\npcdl_{pcdl.__version__}'
    # going home
    return s_version


def make_gif():
    pass


def make_movie():
    pass


def plot_contour():
    pass


def plot_scatter():
    pass


def plot_timeseries():
    # argv
    parser = argparse.ArgumentParser(
        prog = 'pcdl_plot_timeseries',
        description = 'this function to generate a timeseries plot and either returns a matplotlib figure or an image file (jpeg, png, tiff).',
        epilog = 'homepage: https://github.com/elmbeech/physicelldataloader',
    )

    # TimeSeries path
    parser.add_argument(
        'path',
        nargs = '?',
        default = '.',
        help = 'path to the PhysiCell output directory. default is . .',
    )
    # TimeSeries custom_type
    # nop
    # TimeSeries microenv
    parser.add_argument(
        '--microenv',
        nargs = '?',
        default = 'true',
        help = 'should the microenvironment be extracted? setting microenv to False will use less memory and speed up processing, similar to the original pyMCDS_cells.py script. default is True.',
    )
    # TimeSeries graph
    # nop
    # TimeSeries settingxml
    parser.add_argument(
        '--settingxml',
        nargs = '?',
        default = 'PhysiCell_settings.xml',
        help = 'from which settings.xml should the substrate and cell type ID label mapping be extracted? set to None or False if the xml file is missing! default is PhysiCell_settings.xml.',
    )
    # TimeSeries verbose
    parser.add_argument(
        '-v', '--verbose',
        nargs = '?',
        default = 'true',
        help = 'setting verbose to False for less text output, while processing. default is True.',
    )

    # plot_timeseries focus_cat
    parser.add_argument(
        'focus_cat',
        nargs = '?',
        default = 'none',
        help = 'categorical or boolean data column within dataframe specified under frame. default is None, which is total, which is all agents or voxels, no categories. default is None.',
    )
    # plot_timeseries focus_num
    parser.add_argument(
        'focus_num',
        nargs = '?',
        default = 'none',
        help = 'numerical data column within dataframe specified under frame. default is None, which is count, agent or voxel count. default is None.',
    )
    # plot_timeseries aggregate_num
    parser.add_argument(
        'aggregate_num',
        nargs = '?',
        default = 'mean',
        help = 'aggregation function {max, mean, median, min, std, var} for focus_num data. default is mean.',
    )
    # plot_timeseries frame
    parser.add_argument(
        '--frame',
        nargs = '?',
        default = 'cell',
        help = 'to specifies the data dataframe. cell: dataframe will be retrieved through the mcds.get_cell_df function. conc: dataframe will be retrieved through the mcds.get_conc_df function. default is cell.',
    )
    # plot_timeseries z_slice
    parser.add_argument(
        '--z_slice',
        nargs = '?',
        default = 'none',
        help = 'z-axis position to slice a 2D xy-plain out of the 3D mesh. if z_slice position numeric but not an exact mesh center coordinate, then z_slice will be adjusted to the nearest mesh center value, the smaller one, if the coordinate lies on a saddle point. if set to None, the whole domain is taken. default is None.',
    )
    # plot_timeseries logy
    parser.add_argument(
        '--logy',
        nargs = '?',
        default = 'false',
        help = 'if True, then y axis is natural log scaled. default is False.',
    )
    # plot_timeseries ylim
    parser.add_argument(
        '--ylim',
        nargs = '?',
        default = 'none',
        help = 'two floats. y axis min and max value. default is None, which automatically detects min and max value. default is None.',
    )
    # plot_timeseries secondary_y
    parser.add_argument(
        '--secondary_y',
        nargs = '?',
        default = 'false',
        help = 'whether to plot on the secondary y-axis. if a listing of string, which columns to plot on the secondary y-axis. default is False.',
    )
    # plot_timeseries subplots
    # partly nop
    parser.add_argument(
        '--subplots',
        nargs = '?',
        default = 'false',
        help = 'whether to split the plot into subplots, one per column. default is False.',
    )
    # plot_timeseries sharex
    parser.add_argument(
        '--sharex',
        nargs = '?',
        default = 'false',
        help = 'in case subplots is True, share x-axis by setting some x-axis labels to invisible. default is False.',
    )
    # plot_timeseries sharey
    parser.add_argument(
        '--sharey',
        nargs = '?',
        default = 'false',
        help = 'in case subplots is True, share y-axis range and possibly setting some y-axis labels to invisible. default is False.',
    )
    # plot_timeseries linestyle
    parser.add_argument(
        '--linestyle',
        nargs = '?',
        default = '-',
        help = 'matplotlib line style {-, --, -., :} string. default is - .',
    )
    # plot_timeseries linewidth
    parser.add_argument(
        '--linewidth',
        nargs = '?',
        default = 'none',
        help = 'line width in points, integer. default is None.',
    )
    # plot_timeseries cmap
    parser.add_argument(
        '--cmap',
        nargs = '?',
        default = 'none',
        help = 'matplotlib colormap string from https://matplotlib.org/stable/tutorials/colors/colormaps.html . default is None.',
    )
    # plot_timeseries color
    parser.add_argument(
        '--color',
        nargs = '+',
        default = 'none',
        help = 'color string or listing of color string referred to by name, RGB or RGBA code. default is None.',
    )
    # plot_timeseries grid
    parser.add_argument(
        '--grid',
        nargs = '?',
        default = 'true',
        help = 'plot axis grid lines. default is True.',
    )
    # plot_timeseries legend
    parser.add_argument(
        '--legend',
        nargs = '?',
        default = 'true',
        help = 'if True or reverse, place legend on axis subplots. default is True.',
    )
    # plot_timeseries yunit
    parser.add_argument(
        '--yunit',
        nargs = '?',
        default = 'none',
        help = 'string to specify y-axis unit. None will not print a unit on the y-axis. default is None.',
    )
    # plot_timeseries title
    # nop partly
    parser.add_argument(
        '--title',
        nargs = '?',
        default = 'none',
        help = 'title to use for the plot. None will print no title. default is None.',
    )
    # plot_timeseries ax
    # nop
    # plot_timeseries figsizepx
    parser.add_argument(
        '--figsizepx',
        nargs = '?',
        default = '640 480',
        help = 'size of the figure in pixels, x y. the given x and y will be rounded to the nearest even number, to be able to generate movies from the images. default is 640 480.',
    )
    # plot_timeseries ext
    # nop partly
    parser.add_argument(
        '--ext',
        nargs = '?',
        default = 'jpeg',
        help = 'output image format. possible formats are jpeg, png, and tiff. default is jpeg.',
    )
    # plot_timeseries figbgcolor
    parser.add_argument(
        '--figbgcolor',
        nargs = '?',
        default = 'none',
        help = 'figure background color. None is transparent (png) or white (jpeg, tiff). default is None.',
    )

    # parse arguments
    args = parser.parse_args()
    print(args)

    # process arguments
    # aggregate_num
    if (args.aggregate_num == 'entropy'): o_aggregate_num = entropy
    elif (args.aggregate_num == 'max'): o_aggregate_num = np.nanmax
    elif (args.aggregate_num == 'mean'): o_aggregate_num = np.nanmean
    elif (args.aggregate_num == 'median'): o_aggregate_num = np.nanmedian
    elif (args.aggregate_num == 'min'): o_aggregate_num = np.nanmin
    elif (args.aggregate_num == 'std'): o_aggregate_num = np.nanstd
    elif (args.aggregate_num == 'var'): o_aggregate_num = np.nanvar
    else: sys.exit(f"Error @ pcdl_plot_timeseries : unknowen aggregate_num {args.aggregate_num}. knowen are entropy, max, mean, median, min, std, var.")
    # secondary_y
    if (args.secondary_y.lower() == 'false'): ls_secondary_y = False
    elif (args.secondary_y.lower() == 'true'): ls_secondary_y = True
    else: ls_secondary_y = args.secondary_y.split()
    # legend
    if (args.legend.lower() == 'reverse'): b_legend = 'reverse'
    elif args.legend.lower().startswith('f'): b_legend = False
    else: b_legend = True
    # run
    if not os.path.exists(args.path + '/initial.xml'):
        sys.exit(f"Error @ pcdl_plot_timeseries : path does not look like a physicell output directory. {args.path}/initial.xml is missing.")
    mcdsts = pcdl.pyMCDSts(
        output_path = args.path,
        #custom_type,
        load = True,
        microenv = False if args.microenv.lower().startswith('f') else True,
        graph = False,
        settingxml = args.settingxml,
        verbose = False if args.verbose.lower().startswith('f') else True,
    )
    s_pathfile = mcdsts.plot_timeseries(
        focus_cat = None if (args.focus_cat.lower() == 'none') else args.focus_cat,
        focus_num = None if (args.focus_num.lower() == 'none') else args.focus_num,
        aggregate_num = o_aggregate_num,
        frame = args.frame,
        z_slice = None if (args.z_slice.lower() == 'none') else float(args.z_slice),
        logy = True if args.logy.lower().startswith('t') else False,
        ylim = None if (args.ylim.lower() == 'none') else args.ylim.split(),
        secondary_y = ls_secondary_y,
        subplots = True if args.subplots.lower().startswith('t') else False,
        sharex = True if args.sharex.lower().startswith('t') else False,
        sharey = True if args.sharey.lower().startswith('t') else False,
        linestyle = args.linestyle,
        linewidth = None if (args.linewidth.lower() == 'none') else int(args.linewidth),
        cmap = None if (args.cmap.lower() == 'none') else args.cmap,
        #color = None if (args.color.lower() == 'none') else args.color,
        color = args.color if type(args.color) is list else None,
        grid = False if args.grid.lower().startswith('f') else True,
        legend = b_legend,
        yunit = None if (args.yunit.lower() == 'none') else args.yunit,
        title = None if (args.title.lower() == 'none') else args.title,
        ax = None,
        figsizepx = [int(i) for i in args.figsizepx.split()],
        ext = args.ext,
        figbgcolor = None if (args.figbgcolor.lower() == 'none') else args.figbgcolor,
    )
    # going home
    return s_pathfile
