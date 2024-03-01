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
def get_anndata():
    # argv
    parser = argparse.ArgumentParser(
        prog = 'pcdl_get_anndata',
        description = 'function to transform mcds time steps into one or many anndata objects for downstream analysis.',
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
        default = 'true',
        help = 'should the microenvironment be extracted? setting microenv to False will use less memory and speed up processing, similar to the original pyMCDS_cells.py script. default is True.'
    )
    # TimeSeries graph False
    # TimeSeries settingxml
    parser.add_argument(
        '--settingxml',
        default = 'PhysiCell_settings.xml',
        help = 'from which settings.xml should the cell type ID label mapping be extracted? set to None or False if the xml file is missing! default is PhysiCell_settings.xml.',
    )
    # TimeSeries verbose
    parser.add_argument(
        '-v', '--verbose',
        default = 'true',
        help = 'setting verbose to False for less text output, while processing. default is True.'
    )
    # get_anndata values
    parser.add_argument(
        'values',
        nargs = '?',
        default = 1,
        type = int,
        help = 'minimal number of values a variable has to have in any of the mcds time steps to be outputted. variables that have only 1 state carry no information. None is a state too. default is 1.'
    )
    # get_anndata drop
    parser.add_argument(
        '--drop',
        nargs = '*',
        default = '',
        help = "set of column labels to be dropped for the dataframe. don't worry: essential columns like ID, coordinates and time will never be dropped. Attention: when the keep parameter is given, then the drop parameter has to be an empty string! default is an empty string."
    )
    # get_anndata keep
    parser.add_argument(
        '--keep',
        nargs = '*',
        default = '',
        help = "set of column labels to be kept in the dataframe. set values=1 to be sure that all variables are kept. don't worry: essential columns like ID, coordinates and time will always be kept. default is an empty string."
    )
    # get_anndata scale
    parser.add_argument(
        '--scale',
        default = 'maxabs',
        help = "specify how the data should be scaled. possible values are None, maxabs, minmax, std. None: no scaling. set scale to None if you would like to have raw data or entirely scale, transform, and normalize the data later. maxabs: maximum absolute value distance scaler will linearly map all values into a [-1, 1] interval. if the original data has no negative values, the result will be the same as with the minmax scaler (except with features with only one value). if the feature has only zeros, the value will be set to 0. minmax: minimum maximum distance scaler will map all values linearly into a [0, 1] interval. if the feature has only one value, the value will be set to 0. std: standard deviation scaler will result in sigmas. each feature will be mean centered around 0. ddof delta degree of freedom is set to 1 because it is assumed that the values are samples out of the population and not the entire population. it is incomprehensible to me that the equivalent sklearn method has ddof set to 0. if the feature has only one value, the value will be set to 0. default is maxabs")
    # get_anndata collapse
    parser.add_argument(
        '--collapse',
        default = 'true',
        help = 'should all mcds time steps from the time series be collapsed into one big anndata h5ad file, or a many h5ad, one h5ad for each time step?, default is True.'
    )

    # parse arguments
    args = parser.parse_args()
    print(args)

    # process arguments
    s_pathfile = args.path
    if not s_pathfile.endswith('.xml'):
        s_pathfile = s_pathfile + '/initial.xml'
    if not os.path.exists(s_pathfile):
        sys.exit(f"Error @ pcdl_get_anndata : {args.path} path does not look like a outputnnnnnnnn.xml file or physicell output directory ({args.path}/initial.xml is missing).")
    # run
    if os.path.isfile(args.path):
        mcds = pcdl.TimeStep(
            xmlfile = args.path,
            output_path = '.',
            microenv = False if args.microenv.lower().startswith('f') else True,
            graph = False,
            settingxml = None if ((args.settingxml.lower() == 'none') or (args.settingxml.lower() == 'false')) else args.settingxml,
            verbose = False if args.verbose.lower().startswith('f') else True
        )
        ann_mcds = mcds.get_anndata(
            values = args.values,
            drop = set(args.drop),
            keep = set(args.keep),
            scale = args.scale,
        )
        s_opathfile = args.path.replace('.xml','_cell.h5ad')
        ann_mcds.write(s_opathfile)
        # going home
        return s_opathfile
    else:
        mcdsts = pcdl.TimeSeries(
            output_path = args.path,
            #custom_type,
            load = True,
            microenv = False if args.microenv.lower().startswith('f') else True,
            graph = False,
            settingxml = None if ((args.settingxml.lower() == 'none') or (args.settingxml.lower() == 'false')) else args.settingxml,
            verbose = False if args.verbose.lower().startswith('f') else True,
        )
        ann_mcdsts = mcdsts.get_anndata(
            values = args.values,
            drop = set(args.drop),
            keep = set(args.keep),
            scale = args.scale,  # ERROR
            collapse = False if args.collapse.lower().startswith('f') else True,
        )
        if (type(ann_mcdsts) is list):
            ls_opathfile = [f"{args.path}/{s_xmlfile.replace('.xml','_cell.h5ad')}" for s_xmlfile in mcdsts.get_xmlfile_list()]
            for i, ann_mcds in enumerate(ann_mcdsts):
                ann_mcds.write(ls_opathfile[i])
            # going home
            return ls_opathfile
        else:
            s_opathfile = f'{args.path}/timeseries_cell.h5ad'
            ann_mcdsts.write(s_opathfile)
            # going home
            return s_opathfile


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
        default = 'true',
        help = 'should the microenvironment be extracted? setting microenv to False will use less memory and speed up processing, similar to the original pyMCDS_cells.py script. default is True.'
    )
    # TimeSeries graph False
    # TimeSeries settingxml
    parser.add_argument(
        '--settingxml',
        default = 'PhysiCell_settings.xml',
        help = 'from which settings.xml should the cell type ID label mapping be extracted? set to None or False if the xml file is missing! default is PhysiCell_settings.xml.',
    )
    # TimeSeries verbose
    parser.add_argument(
        '-v', '--verbose',
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
        nargs = '*',
        default = '',
        help = "set of column labels to be dropped for the dataframe. don't worry: essential columns like ID, coordinates and time will never be dropped. Attention: when the keep parameter is given, then the drop parameter has to be an empty string! default is an empty string."
    )
    # get_cell_df keep
    parser.add_argument(
        '--keep',
        nargs = '*',
        default = '',
        help = "set of column labels to be kept in the dataframe. set values=1 to be sure that all variables are kept. don't worry: essential columns like ID, coordinates and time will always be kept. default is an empty string."
    )
    # get_cell_df collapse
    parser.add_argument(
        '--collapse',
        default = 'true',
        help = 'should all mcds time steps from the time series be collapsed into one big csv, or a many csv, one csv for each time step?, default is True.'
    )

    # parse arguments
    args = parser.parse_args()
    print(args)

    # process arguments
    s_pathfile = args.path
    if not s_pathfile.endswith('.xml'):
        s_pathfile = s_pathfile + '/initial.xml'
    if not os.path.exists(s_pathfile):
        sys.exit(f"Error @ pcdl_get_cell_df : {args.path} path does not look like a outputnnnnnnnn.xml file or physicell output directory ({args.path}/initial.xml is missing).")
    # run
    if os.path.isfile(args.path):
        mcds = pcdl.pyMCDS(
            xmlfile = args.path,
            output_path = '.',
            microenv = False if args.microenv.lower().startswith('f') else True,
            graph = False,
            settingxml = None if ((args.settingxml.lower() == 'none') or (args.settingxml.lower() == 'false')) else args.settingxml,
            verbose = False if args.verbose.lower().startswith('f') else True
        )
        df_cell = mcds.get_cell_df(
            values = args.values,
            drop = set(args.drop),
            keep = set(args.keep),
        )
        s_opathfile = args.path.replace('.xml','_cell.csv')
        df_cell.to_csv(s_opathfile)
        # going home
        return s_opathfile
    else:
        mcdsts = pcdl.pyMCDSts(
            output_path = args.path,
            #custom_type,
            load = True,
            microenv = False if args.microenv.lower().startswith('f') else True,
            graph = False,
            settingxml = None if ((args.settingxml.lower() == 'none') or (args.settingxml.lower() == 'false')) else args.settingxml,
            verbose = False if args.verbose.lower().startswith('f') else True,
        )
        ldf_cell = mcdsts.get_cell_df(
            values = args.values,
            drop = set(args.drop),
            keep = set(args.keep),
            collapse = False if args.collapse.lower().startswith('f') else True,
        )
        if type(ldf_cell) is list:
            ls_opathfile = [f"{args.path}/{s_xmlfile.replace('.xml','_cell.csv')}" for s_xmlfile in mcdsts.get_xmlfile_list()]
            for i, df_cell in enumerate(ldf_cell):
                df_cell.to_csv(ls_opathfile[i])
            # going home
            return ls_opathfile
        else:
            s_opathfile = f'{args.path}/timeseries_cell.csv'
            ldf_cell.to_csv(s_opathfile)
            # going home
            return s_opathfile


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
        default = 'true',
        help = 'should the microenvironment be extracted? setting microenv to False will use less memory and speed up processing, similar to the original pyMCDS_cells.py script. default is True.',
    )
    # TimeSeries graph False
    # TimeSeries settingxml
    parser.add_argument(
        '--settingxml',
        default = 'PhysiCell_settings.xml',
        help = 'from which settings.xml should the cell type ID label mapping be extracted? set to None or False if the xml file is missing! default is PhysiCell_settings.xml.',
    )
    # TimeSeries verbose
    parser.add_argument(
        '-v', '--verbose',
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
        nargs = '*',
        default = '',
        help = "set of column labels to be dropped for the dataframe. don't worry: essential columns like ID, coordinates and time will never be dropped. Attention: when the keep parameter is given, then the drop parameter has to be an empty string! default is an empty string.",
    )
    # get_cell_df_features keep
    parser.add_argument(
        '--keep',
        nargs = '*',
        default = '',
        help = "set of column labels to be kept in the dataframe. set values=1 to be sure that all variables are kept. don't worry: essential columns like ID, coordinates and time will always be kept. default is an empty string.",
    )
    # get_cell_df_features allvalues
    parser.add_argument(
        '--allvalues',
        default = 'false',
        help = 'for numeric data, should only the min and max values or all values be returned? default is false.',
    )

    # parse arguments
    args = parser.parse_args()
    print(args)

    # process arguments
    if not os.path.exists(args.path + '/initial.xml'):
        sys.exit(f"Error @ pcdl_get_cell_df_features : {args.path} path does not look like a physicell output directory ({args.path}/initial.xml is missing).")
    mcdsts = pcdl.pyMCDSts(
        output_path = args.path,
        #custom_type,
        load = True,
        microenv = False if args.microenv.lower().startswith('f') else True,
        graph = False,
        settingxml = None if ((args.settingxml.lower() == 'none') or (args.settingxml.lower() == 'false')) else args.settingxml,
        verbose = False if args.verbose.lower().startswith('f') else True,
    )
    # run
    s_values = 'minmax'
    b_allvalues = True if args.allvalues.lower().startswith('t') else False
    if b_allvalues:
        s_values = 'all'
    dl_variable = mcdsts.get_cell_df_features(
        values = args.values,
        drop = set(args.drop),
        keep = set(args.keep),
        allvalues = b_allvalues,
    )
    s_opathfile = f'{args.path}/timeseries_cell_features_{s_values}.json'
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
    # TimeSeries settingxml None
    # TimeSeries verbose
    parser.add_argument(
        '-v', '--verbose',
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
        nargs = '*',
        default = '',
        help = "set of column labels to be dropped for the dataframe. don't worry: essential columns like ID, coordinates and time will never be dropped. Attention: when the keep parameter is given, then the drop parameter has to be an empty string! default is an empty string.",
    )
    # get_conc_df keep
    parser.add_argument(
        '--keep',
        nargs = '*',
        default = '',
        help = "set of column labels to be kept in the dataframe. set values=1 to be sure that all variables are kept. don't worry: essential columns like ID, coordinates and time will always be kept. default is an empty string.",
    )
    # get_conc_df collapse
    parser.add_argument(
        '--collapse',
        default = 'true',
        help = 'should all mcds time steps from the time series be collapsed into one big csv, or a many csv, one for each time step?, default is True.'
    )

    # parse arguments
    args = parser.parse_args()
    print(args)

    # process arguments
    s_pathfile = args.path
    if not s_pathfile.endswith('.xml'):
        s_pathfile = s_pathfile + '/initial.xml'
    if not os.path.exists(s_pathfile):
        sys.exit(f"Error @ pcdl_get_conc_df : {args.path} path does not look like a outputnnnnnnnn.xml file or physicell output directory ({args.path}/initial.xml is missing).")
    # run
    if os.path.isfile(args.path):
        mcds = pcdl.pyMCDS(
            xmlfile = args.path,
            output_path = '.',
            microenv = True,
            graph = False,
            settingxml = None,
            verbose = False if args.verbose.lower().startswith('f') else True
        )
        df_conc = mcds.get_conc_df(
            values = args.values,
            drop = set(args.drop),
            keep = set(args.keep),
        )
        s_opathfile = args.path.replace('.xml','_conc.csv')
        df_conc.to_csv(s_opathfile, index=False)
        # going home
        return s_opathfile
    else:
        mcdsts = pcdl.pyMCDSts(
            output_path = args.path,
            #custom_type,
            load = True,
            microenv = True,
            graph = False,
            settingxml = None,
            verbose = False if args.verbose.lower().startswith('f') else True,
        )
        ldf_conc = mcdsts.get_conc_df(
            values = args.values,
            drop = set(args.drop),
            keep = set(args.keep),
            collapse = False if args.collapse.lower().startswith('f') else True,
        )
        if type(ldf_conc) is list:
            ls_opathfile = [f"{args.path}/{s_xmlfile.replace('.xml','_conc.csv')}" for s_xmlfile in mcdsts.get_xmlfile_list()]
            for i, df_conc in enumerate(ldf_conc):
                df_conc.to_csv(ls_opathfile[i], index=False)
            # going home
            return ls_opathfile
        else:
            s_opathfile = f'{args.path}/timeseries_conc.csv'
            ldf_conc.to_csv(s_opathfile, index=False)
            # going home
            return s_opathfile


def get_conc_df_features():
    # argv
    parser = argparse.ArgumentParser(
        prog = 'pcdl_get_conc_df_features',
        description = 'function to detect informative substrate concentration variables in a time series. this function detects even variables which have less than the minimal state count in each time step, but different values from time step to time step. the output is a json file with an entry of all non-coordinate column names that at least in one of the time steps or in between time steps, reach the given minimal value count. key is the column name, mapped is a list of all values (bool, str, and, if allvalues is True, int and float) or a list with minimum and maximum values (int, float).',
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
    # TimeSeries microenv True
    # TimeSeries graph False
    # TimeSeries settingxml None
    # TimeSeries verbose
    parser.add_argument(
        '-v', '--verbose',
        default = 'true',
        help = 'setting verbose to False for less text output, while processing. default is True.',
    )
    # get_conc_df_features values
    parser.add_argument(
        'values',
        nargs = '?',
        default = 1,
        type = int,
        help = 'minimal number of values a variable has to have in any of the mcds time steps to be outputted. variables that have only 1 state carry no information. None is a state too. default is 1.',
    )
    # get_conc_df_features drop
    parser.add_argument(
        '--drop',
        nargs = '*',
        default = '',
        help = "set of column labels to be dropped for the dataframe. don't worry: essential columns like ID, coordinates and time will never be dropped. Attention: when the keep parameter is given, then the drop parameter has to be an empty string! default is an empty string.",
    )
    # get_conc_df_features keep
    parser.add_argument(
        '--keep',
        nargs = '*',
        default = '',
        help = "set of column labels to be kept in the dataframe. set values=1 to be sure that all variables are kept. don't worry: essential columns like ID, coordinates and time will always be kept. default is an empty string.",
    )
    # get_conc_df_features allvalues
    parser.add_argument(
        '--allvalues',
        default = 'false',
        help = 'for numeric data, should only the min and max values or all values be returned? default is false.',
    )

    # parse arguments
    args = parser.parse_args()
    print(args)

    # process arguments
    if not os.path.exists(args.path + '/initial.xml'):
        sys.exit(f"Error @ pcdl_get_conc_df_features : {args.path} path does not look like a physicell output directory ({args.path}/initial.xml is missing).")
    mcdsts = pcdl.pyMCDSts(
        output_path = args.path,
        #custom_type,
        load = True,
        microenv = True,
        graph = False,
        settingxml = None,
        verbose = False if args.verbose.lower().startswith('f') else True,
    )
    # run
    s_values = 'minmax'
    b_allvalues = True if args.allvalues.lower().startswith('t') else False
    if b_allvalues:
        s_values = 'all'
    dl_variable = mcdsts.get_conc_df_features(
        values = args.values,
        drop = set(args.drop),
        keep = set(args.keep),
        allvalues = b_allvalues,
    )
    s_opathfile = f'{args.path}/timeseries_conc_features_{s_values}.json'
    json.dump(dl_variable, open(s_opathfile, 'w'), sort_keys=True)
    # going home
    return s_opathfile


def get_unit_se():
    # argv
    parser = argparse.ArgumentParser(
        prog = 'pcdl_get_unit_se',
        description = 'function returns a csv that lists all tracked variables from metadata, cell, and microenvironment and maps them to their unit.',
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
        default = 'true',
        help = 'should the microenvironment be extracted? setting microenv to False will use less memory and speed up processing, similar to the original pyMCDS_cells.py script. default is True.',
    )
    # TimeSeries graph False
    # TimeSeries settingxml
    parser.add_argument(
        '--settingxml',
        default = 'PhysiCell_settings.xml',
        help = 'from which settings.xml should the cell type ID label mapping be extracted? set to None or False if the xml file is missing! default is PhysiCell_settings.xml.',
    )
    # TimeSeries verbose
    parser.add_argument(
        '-v', '--verbose',
        default = 'true',
        help = 'setting verbose to False for less text output, while processing. default is True.',
    )

    # parse arguments
    args = parser.parse_args()
    print(args)

    # process arguments
    s_path = args.path
    s_pathfile = args.path
    if s_pathfile.endswith('.xml'):
        s_path = '/'.join(args.path.replace('\\','/').split('/')[:-1])
    else:
        s_pathfile = s_pathfile + '/initial.xml'
    if not os.path.exists(s_pathfile):
        sys.exit(f"Error @ pcdl_get_conc_df : {args.path} path does not look like a outputnnnnnnnn.xml file or physicell output directory ({args.path}/initial.xml is missing).")
    # run
    mcds = pcdl.pyMCDS(
        xmlfile = s_pathfile,
        output_path = '.',
        microenv = False if args.microenv.lower().startswith('f') else True,
        graph = False,
        settingxml = None if ((args.settingxml.lower() == 'none') or (args.settingxml.lower() == 'false')) else args.settingxml,
        verbose = True if args.verbose.lower().startswith('t') else False
    )
    se_unit = mcds.get_unit_se()
    s_opathfile = f'{s_path}/timeseries_unit.csv'
    se_unit.to_csv(s_opathfile)
    # going home
    return s_opathfile


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
        sys.exit(f"Error @ pcdl_get_version : {args.path} path does not look like a outputnnnnnnnn.xml file or physicell output directory ({args.path}/initial.xml is missing).")
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


def make_graph_gml():
    # argv
    parser = argparse.ArgumentParser(
        prog = 'pcdl_make_graph_gml',
        description = 'function to generate graph files in the gml graph modelling language standard format. gml was the outcome of an initiative that started at the international symposium on graph drawing 1995 in Passau and ended at Graph Drawing 1996 in Berkeley. the networkx python and igraph C and python libraries for graph analysis are gml compatible and can as such read and write this file format.',
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
    # TimeSeries microenv
    parser.add_argument(
        '--microenv',
        default = 'true',
        help = 'should the microenvironment be extracted? setting microenv to False will use less memory and speed up processing, similar to the original pyMCDS_cells.py script. default is True.'
    )
    # TimeSeries graph True
    # TimeSeries settingxml
    parser.add_argument(
        '--settingxml',
        default = 'PhysiCell_settings.xml',
        help = 'from which settings.xml should the cell type ID label mapping be extracted? set to None or False if the xml file is missing! default is PhysiCell_settings.xml.',
    )
    # TimeSeries verbose
    parser.add_argument(
        '-v', '--verbose',
        default = 'true',
        help = 'setting verbose to False for less text output, while processing. default is True.',
    )
    # make_graph_gml graph_type
    parser.add_argument(
        'graph_type',
        nargs = '?',
        default = 'neighbor',
        help = 'to specify which physicell output data should be processed. attached: processes mcds.get_attached_graph_dict dictionary. neighbor: processes mcds.get_neighbor_graph_dict dictionary. default is neighbor.',
    )
    # make_graph_gml edge_attr
    parser.add_argument(
        '--edge_attr',
        default = 'true',
        help = 'specifies if the spatial Euclidean distance is used for edge attribute, to generate a weighted graph. default is True.',
    )
    # make_graph_gml node_attr
    parser.add_argument(
        '--node_attr',
        nargs = '*',
        default = [],
        help = 'listing of mcds.get_cell_df dataframe columns, used for node attributes. default is and empty list.',
    )

    # parse arguments
    args = parser.parse_args()
    print(args)

    # process arguments
    s_pathfile = args.path
    if not s_pathfile.endswith('.xml'):
        s_pathfile = s_pathfile + '/initial.xml'
    if not os.path.exists(s_pathfile):
        sys.exit(f"Error @ pcdl_make_graph_gml : {args.path} path does not look like a outputnnnnnnnn.xml file or physicell output directory ({args.path}/initial.xml is missing).")
    if os.path.isfile(args.path):
        # run
        mcds = pcdl.pyMCDS(
            xmlfile = args.path,
            output_path = '.',
            microenv = False if args.microenv.lower().startswith('f') else True,
            graph = True,
            settingxml = None if ((args.settingxml.lower() == 'none') or (args.settingxml.lower() == 'false')) else args.settingxml,
            verbose = False if args.verbose.lower().startswith('f') else True
        )
        s_opathfile = df_conc = mcds.make_graph_gml(
            graph_type = args.graph_type,
            edge_attr = False if args.edge_attr.lower().startswith('f') else True,
            node_attr = args.node_attr,
        )
        # going home
        return s_opathfile
    else:
        # run
        mcdsts = pcdl.pyMCDSts(
            output_path = args.path,
            #custom_type,
            load = True,
            microenv = False if args.microenv.lower().startswith('f') else True,
            graph = True,
            settingxml = None if ((args.settingxml.lower() == 'none') or (args.settingxml.lower() == 'false')) else args.settingxml,
            verbose = False if args.verbose.lower().startswith('f') else True,
        )
        ls_opathfile = []
        for mcds in mcdsts.get_mcds_list():
            s_opathfile = df_conc = mcds.make_graph_gml(
                graph_type = args.graph_type,
                edge_attr = False if args.edge_attr.lower().startswith('f') else True,
                node_attr = args.node_attr,
            )
            ls_opathfile.append(s_opathfile)
        # going home
        return ls_opathfile


def make_gif():
    # argv
    parser = argparse.ArgumentParser(
        prog = 'pcdl_make_gif',
        description = 'this function generates a gif movie from all image files found in the path directory in the specified interface file format.',
        epilog = 'homepage: https://github.com/elmbeech/physicelldataloader',
    )

    # TimeSeries path
    parser.add_argument(
        'path',
        nargs = '?',
        default = '.',
        help = 'relative or absolute path to where the images are from which the gif will be generated. default is . .',
    )
    # make_movie interface
    parser.add_argument(
        'interface',
        nargs = '?',
        default = 'jpeg',
        help = 'specify the image format from which the gif will be generated. these images have to exist under the given path. they can be generated with the plot_scatter or plot_contour function. default is jpeg.',
    )

    # parse arguments
    args = parser.parse_args()
    print(args)

    # process arguments
    # run
    s_opathfile = pcdl.make_gif(
        path = args.path,
        interface = args.interface,
    )
    # going home
    return s_opathfile


def make_movie():
    # argv
    parser = argparse.ArgumentParser(
        prog = 'pcdl_make_movie',
        description = 'this function generates an mp4 movie file from all image files found in the path directory in the specified interface file format.',
        epilog = 'homepage: https://github.com/elmbeech/physicelldataloader',
    )

    # TimeSeries path
    parser.add_argument(
        'path',
        nargs = '?',
        default = '.',
        help = 'relative or absolute path to where the images are from which the mp4 movie will be generated. default is . .',
    )
    # make_movie interface
    parser.add_argument(
        'interface',
        nargs = '?',
        default = 'jpeg',
        help = 'specify the image format from which the mp4 movie will be generated. these images have to exist under the given path. they can be generated with the plot_scatter or plot_contour function. default is jpeg.',
    )
    # make_movie framerate
    parser.add_argument(
        '--framerate',
        default = 12,
        type = int,
        help = 'specifies how many images per second will be used. humans are capable of processing 12 images per second and seeing them individually. higher rates are seen as motion. default is 12.',
    )

    # parse arguments
    args = parser.parse_args()
    print(args)

    # process arguments
    # run
    s_opathfile = pcdl.make_movie(
        path = args.path,
        interface = args.interface,
        framerate = args.framerate,
    )
    # going home
    return s_opathfile


def plot_contour():
    # argv
    parser = argparse.ArgumentParser(
        prog = 'pcdl_plot_contour',
        description = 'function generates matplotlib contour (or contourf) plots, inclusive color bar, for the substrate specified, under the returned path.',
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
    # TimeSeries settingxml None
    # TimeSeries verbose
    parser.add_argument(
        '-v', '--verbose',
        default = 'true',
        help = 'setting verbose to False for less text output, while processing. default is True.',
    )
    # plot_contour focus
    parser.add_argument(
        'focus',
        nargs = '?',
        help = 'column name within conc dataframe.',
    )
    # plot_contour z_slice
    parser.add_argument(
        '--z_slice',
        default = 0.0,
        type = float,
        help = 'z-axis position to slice a 2D xy-plain out of the 3D mesh. if z_slice position numeric but not an exact mesh center coordinate, then z_slice will be adjusted to the nearest mesh center value, the smaller one, if the coordinate lies on a saddle point. default is 0.0.',
    )
    # plot_contour extrema
    parser.add_argument(
        '--extrema',
        nargs = '+',
        default = ['none'],
        help = 'listing of two floats. None takes min and max from data. default is None.',
    )
    # plot_contour alpha
    parser.add_argument(
        '--alpha',
        default = 1.0,
        type = float,
        help = 'alpha channel transparency value between 1 (not transparent at all) and 0 (totally transparent). default is 1.0.',
    )
    # plot_contour fill
    parser.add_argument(
        '--fill',
        default = 'true',
        help = 'True generates a matplotlib contourf plot. False generates a matplotlib contour plot. default is True.',
    )
    # plot_contour cmap
    parser.add_argument(
        '--cmap',
        default = 'viridis',
        help = 'matplotlib colormap string from https://matplotlib.org/stable/tutorials/colors/colormaps.html . default is viridis.',
    )
    # plot_contour grid
    parser.add_argument(
        '--grid',
        default = 'true',
        help = 'plot axis grid lines. default is True.',
    )
    # plot_contour xlim
    parser.add_argument(
        '--xlim',
        nargs = '+',
        default = ['none'],
        help = 'two floats. x axis min and max value. None takes min and max from mesh x axis range. default is None.',
    )
    # plot_contour ylim
    parser.add_argument(
        '--ylim',
        nargs = '+',
        default = ['none'],
        help = 'two floats. y axis min and max value. None takes min and max from mesh y axis range. default is None.',
    )
    # plot_contour xyequal
    parser.add_argument(
        '--xyequal',
        default = 'true',
        help = 'to specify equal axis spacing for x and y axis. default is true.',
    )
    # plot_contour figsizepx
    parser.add_argument(
        '--figsizepx',
        nargs = '+',
        default = ['none'],
        help = 'size of the figure in pixels (integer), x y. the given x and y will be rounded to the nearest even number, to be able to generate movies from the images. None tries to take the values from the initial.svg file. fall back setting is 640 480. default is None.',
    )
    # plot_contour ext
    parser.add_argument(
        '--ext',
        default = 'jpeg',
        help = 'output image format. possible formats are jpeg, png, and tiff. default is jpeg.',
    )
    # plot_contour figbgcolor
    parser.add_argument(
        '--figbgcolor',
        default = 'none',
        help = 'figure background color. None is transparent (png) or white (jpeg, tiff). default is None.',
    )

    # parse arguments
    args = parser.parse_args()
    print(args)

    # process arguments
    args.path = args.path.replace('\\','/')
    s_pathfile = args.path
    if s_pathfile.endswith('.xml'):
        args.path = '/'.join(s_pathfile.split('/')[:-1])
    else:
        s_pathfile = s_pathfile + '/initial.xml'
    if not os.path.exists(s_pathfile):
        sys.exit(f"Error @ pcdl_plot_contour : {args.path} path does not look like a outputnnnnnnnn.xml file or physicell output directory ({args.path}/initial.xml is missing).")
    if (args.focus == None):
        sys.exit(f"Error @ pcdl_plot_contour : input for positional argument focus is missung! this has to be a column name from the conc dataframe.")
    # run
    mcdsts = pcdl.pyMCDSts(
        output_path = args.path,
        #custom_type,
        load = False,
        microenv = True,
        graph = False,
        settingxml = None,
        verbose = False if args.verbose.lower().startswith('f') else True,
    )
    if s_pathfile.endswith('/initial.xml'):
        ls_xmlfile = mcdsts.get_xmlfile_list()
    else:
        ls_xmlfile = [s_pathfile.split('/')[-1]]
    mcdsts.read_mcds(xmlfile_list=ls_xmlfile)
    s_opath = mcdsts.plot_contour(
        focus = args.focus,
        z_slice = args.z_slice,
        extrema = None if (args.extrema[0].lower() == 'none') else args.extrema,
        alpha = args.alpha,
        fill = False if args.fill.lower().startswith('f') else True,
        cmap = args.cmap,
        grid = False if args.grid.lower().startswith('f') else True,
        xlim = None if (args.xlim[0].lower() == 'none') else args.xlim,
        ylim = None if (args.ylim[0].lower() == 'none') else args.ylim,
        xyequal = False if args.xyequal.lower().startswith('f') else True,
        figsizepx = None if (args.figsizepx[0].lower() == 'none') else [int(n) for n in args.figsizepx],
        ext = args.ext,
        figbgcolor = None if (args.figbgcolor.lower() == 'none') else args.figbgcolor,
    )
    # going home
    return s_opath


def plot_scatter():
    # argv
    parser = argparse.ArgumentParser(
        prog = 'pcdl_plot_scatter',
        description = 'function generates pandas scatter plots, under the returned path.',
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
    # TimeSeries microenv
    parser.add_argument(
        '--microenv',
        default = 'true',
        help = 'should the microenvironment be extracted? setting microenv to False will use less memory and speed up processing, similar to the original pyMCDS_cells.py script. default is True.',
    )
    # TimeSeries graph False
    # TimeSeries settingxml
    parser.add_argument(
        '--settingxml',
        default = 'PhysiCell_settings.xml',
        help = 'from which settings.xml should the cell type ID label mapping be extracted? set to None or False if the xml file is missing! default is PhysiCell_settings.xml.',
    )
    # TimeSeries verbose
    parser.add_argument(
        '-v', '--verbose',
        default = 'true',
        help = 'setting verbose to False for less text output, while processing. default is True.',
    )
    # plot_scatter focus
    parser.add_argument(
        'focus',
        nargs = '?',
        default = 'cell_type',
        help = 'column name within conc dataframe. default is cell_type.',
    )
    # plot_scatter z_slice
    parser.add_argument(
        '--z_slice',
        default = 0.0,
        type = float,
        help = 'z-axis position to slice a 2D xy-plain out of the 3D mesh. if z_slice position numeric but not an exact mesh center coordinate, then z_slice will be adjusted to the nearest mesh center value, the smaller one, if the coordinate lies on a saddle point. default is 0.0.',
    )
    # plot_scatter z_axis
    parser.add_argument(
        '--z_axis',
        nargs = '+',
        default = ['none'],
        help = 'for a categorical focus: list of labels; for a numeric focus: list of two floats; None, depending on the focus column variable dtype, extracts labels or min and max values from data. default is None',
    )
    # plot_contour alpha
    parser.add_argument(
        '--alpha',
        default = 1.0,
        type = float,
        help = 'alpha channel transparency value between 1 (not transparent at all) and 0 (totally transparent). default is 1.0.',
    )
    # plot_scatter cmap
    # nop partly
    parser.add_argument(
        '--cmap',
        default = 'viridis',
        help = 'matplotlib colormap string from https://matplotlib.org/stable/tutorials/colors/colormaps.html . default is viridis.',
    )
    # plot_scatter grid
    parser.add_argument(
        '--grid',
        default = 'true',
        help = 'plot axis grid lines. default is True.',
    )
    # plot_scatter legend_loc
    parser.add_argument(
        '--legend_loc',
        default = 'lower left',
        help = "the location of the categorical legend, if applicable. possible strings are: best, 'upper right', 'upper center', 'upper left', 'center left', 'lower left', 'lower center', 'lower right', 'center right', center. default is 'lower left'",
    )
    # plot_scatter xlim
    parser.add_argument(
        '--xlim',
        nargs = '+',
        default = ['none'],
        help = 'two floats. x axis min and max value. None takes min and max from mesh x axis range. default is None.',
    )
    # plot_scatter ylim
    parser.add_argument(
        '--ylim',
        nargs = '+',
        default = ['none'],
        help = 'two floats. y axis min and max value. None takes min and max from mesh y axis range. default is None.',
    )
    # plot_scatter xyequal
    parser.add_argument(
        '--xyequal',
        default = 'true',
        help = 'to specify equal axis spacing for x and y axis. default is True.',
    )
    # plot_scatter s
    parser.add_argument(
        '--s',
        default = 'none',
        help = "scatter plot dot size in pixel. typographic points are 1/72 inch. the marker size s is specified in points**2. plt.rcParams['lines.markersize']**2 is in my case 36. None tries to take the value from the initial.svg file. fall back setting is 36. default is None.",
    )
    # plot_scatter figsizepx
    parser.add_argument(
        '--figsizepx',
        nargs = '+',
        default = ['none'],
        help = 'size of the figure in pixels (integer), x y. the given x and y will be rounded to the nearest even number, to be able to generate movies from the images. None tries to take the values from the initial.svg file. fall back setting is 640 480. default is None.',
    )
    # plot_scatter ext
    parser.add_argument(
        '--ext',
        default = 'jpeg',
        help = 'output image format. possible formats are jpeg, png, and tiff. default is jpeg.',
    )
    # plot_scatter figbgcolor
    parser.add_argument(
        '--figbgcolor',
        default = 'none',
        help = 'figure background color. None is transparent (png) or white (jpeg, tiff). default is None.',
    )

    # parse arguments
    args = parser.parse_args()
    print(args)

    # process arguments
    args.path = args.path.replace('\\','/')
    s_pathfile = args.path
    if s_pathfile.endswith('.xml'):
        args.path = '/'.join(s_pathfile.split('/')[:-1])
    else:
        s_pathfile = s_pathfile + '/initial.xml'
    if not os.path.exists(s_pathfile):
        sys.exit(f"Error @ plot_scatter : {args.path} path does not look like a outputnnnnnnnn.xml file or physicell output directory ({args.path}/initial.xml is missing).")
    # run
    mcdsts = pcdl.pyMCDSts(
        output_path = args.path,
        #custom_type,
        load = False,
        microenv = False if args.microenv.lower().startswith('f') else True,
        graph = False,
        settingxml = None if ((args.settingxml.lower() == 'none') or (args.settingxml.lower() == 'false')) else args.settingxml,
        verbose = False if args.verbose.lower().startswith('f') else True,
    )
    if s_pathfile.endswith('/initial.xml'):
        ls_xmlfile = mcdsts.get_xmlfile_list()
    else:
        ls_xmlfile = [s_pathfile.split('/')[-1]]
    mcdsts.read_mcds(xmlfile_list=ls_xmlfile)
    s_opath = mcdsts.plot_scatter(
        focus = args.focus,
        z_slice = args.z_slice,
        z_axis = None if (args.z_axis[0].lower() == 'none') else args.z_axis,
        alpha = args.alpha,
        cmap = args.cmap,
        grid = False if args.grid.lower().startswith('f') else True,
        legend_loc = args.legend_loc,
        xlim = None if (args.xlim[0].lower() == 'none') else args.xlim,
        ylim = None if (args.ylim[0].lower() == 'none') else args.ylim,
        xyequal = False if args.xyequal.lower().startswith('f') else True,
        s = None if (args.s.lower() == 'none') else int(args.s),
        figsizepx = None if (args.figsizepx[0].lower() == 'none') else [int(i) for i in args.figsizepx],
        ext = args.ext,
        figbgcolor = None if (args.figbgcolor.lower() == 'none') else args.figbgcolor,
    )
    # going home
    return s_opath


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
        default = 'true',
        help = 'should the microenvironment be extracted? setting microenv to False will use less memory and speed up processing, similar to the original pyMCDS_cells.py script. default is True.',
    )
    # TimeSeries graph
    # nop
    # TimeSeries settingxml
    parser.add_argument(
        '--settingxml',
        default = 'PhysiCell_settings.xml',
        help = 'from which settings.xml should the cell type ID label mapping be extracted? set to None or False if the xml file is missing! default is PhysiCell_settings.xml.',
    )
    # TimeSeries verbose
    parser.add_argument(
        '-v', '--verbose',
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
        default = 'cell',
        help = 'to specifies the data dataframe. cell: dataframe will be retrieved through the mcds.get_cell_df function. conc: dataframe will be retrieved through the mcds.get_conc_df function. default is cell.',
    )
    # plot_timeseries z_slice
    parser.add_argument(
        '--z_slice',
        default = 'none',
        help = 'z-axis position to slice a 2D xy-plain out of the 3D mesh. if z_slice position numeric but not an exact mesh center coordinate, then z_slice will be adjusted to the nearest mesh center value, the smaller one, if the coordinate lies on a saddle point. if set to None, the whole domain is taken. default is None.',
    )
    # plot_timeseries logy
    parser.add_argument(
        '--logy',
        default = 'false',
        help = 'if True, then y axis is natural log scaled. default is False.',
    )
    # plot_timeseries ylim
    parser.add_argument(
        '--ylim',
        nargs = '+',
        default = ['none'],
        help = 'two floats. y axis min and max value. default is None, which automatically detects min and max value. default is None.',
    )
    # plot_timeseries secondary_y
    parser.add_argument(
        '--secondary_y',
        nargs = '+',
        default = ['false'],
        help = 'whether to plot on the secondary y-axis. if a listing of string, which columns to plot on the secondary y-axis. default is False.',
    )
    # plot_timeseries subplots
    # nop partly
    parser.add_argument(
        '--subplots',
        default = 'false',
        help = 'whether to split the plot into subplots, one per column. default is False.',
    )
    # plot_timeseries sharex
    parser.add_argument(
        '--sharex',
        default = 'false',
        help = 'in case subplots is True, share x-axis by setting some x-axis labels to invisible. default is False.',
    )
    # plot_timeseries sharey
    parser.add_argument(
        '--sharey',
        default = 'false',
        help = 'in case subplots is True, share y-axis range and possibly setting some y-axis labels to invisible. default is False.',
    )
    # plot_timeseries linestyle
    parser.add_argument(
        '--linestyle',
        default = '-',
        help = 'matplotlib line style {-, --, .-, :} string. default is - .',
    )
    # plot_timeseries linewidth
    parser.add_argument(
        '--linewidth',
        default = 'none',
        help = 'line width in points, integer. default is None.',
    )
    # plot_timeseries cmap
    parser.add_argument(
        '--cmap',
        default = 'none',
        help = 'matplotlib colormap string from https://matplotlib.org/stable/tutorials/colors/colormaps.html . default is None.',
    )
    # plot_timeseries color
    parser.add_argument(
        '--color',
        nargs = '+',
        default = ['none'],
        help = 'color string or listing of color string referred to by name, RGB or RGBA code. default is None.',
    )
    # plot_timeseries grid
    parser.add_argument(
        '--grid',
        default = 'true',
        help = 'plot axis grid lines. default is True.',
    )
    # plot_timeseries legend
    parser.add_argument(
        '--legend',
        default = 'true',
        help = 'if True or reverse, place legend on axis subplots. default is True.',
    )
    # plot_timeseries yunit
    parser.add_argument(
        '--yunit',
        default = 'none',
        help = 'string to specify y-axis unit. None will not print a unit on the y-axis. default is None.',
    )
    # plot_timeseries title
    # nop partly
    parser.add_argument(
        '--title',
        default = 'none',
        help = 'title to use for the plot. None will print no title. default is None.',
    )
    # plot_timeseries ax
    # nop
    # plot_timeseries figsizepx
    parser.add_argument(
        '--figsizepx',
        nargs = '+',
        default = ['640', '480'],
        help = 'size of the figure in pixels (integer), x y. the given x and y will be rounded to the nearest even number, to be able to generate movies from the images. default is 640 480.',
    )
    # plot_timeseries ext
    # nop partly
    parser.add_argument(
        '--ext',
        default = 'jpeg',
        help = 'output image format. possible formats are jpeg, png, and tiff. default is jpeg.',
    )
    # plot_timeseries figbgcolor
    parser.add_argument(
        '--figbgcolor',
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
    if (args.secondary_y[0].lower() == 'false'): ls_secondary_y = False
    elif (args.secondary_y[0].lower() == 'true'): ls_secondary_y = True
    else: ls_secondary_y = args.secondary_y
    # legend
    if (args.legend.lower() == 'reverse'): b_legend = 'reverse'
    elif args.legend.lower().startswith('f'): b_legend = False
    else: b_legend = True
    # run
    if not os.path.exists(args.path + '/initial.xml'):
        sys.exit(f"Error @ pcdl_plot_timeseries : path does not look like a physicell output directory ({args.path}/initial.xml is missing).")
    mcdsts = pcdl.pyMCDSts(
        output_path = args.path,
        #custom_type,
        load = True,
        microenv = False if args.microenv.lower().startswith('f') else True,
        graph = False,
        settingxml = None if ((args.settingxml.lower() == 'none') or (args.settingxml.lower() == 'false')) else args.settingxml,
        verbose = False if args.verbose.lower().startswith('f') else True,
    )
    s_pathfile = mcdsts.plot_timeseries(
        focus_cat = None if (args.focus_cat.lower() == 'none') else args.focus_cat,
        focus_num = None if (args.focus_num.lower() == 'none') else args.focus_num,
        aggregate_num = o_aggregate_num,
        frame = args.frame,
        z_slice = None if (args.z_slice.lower() == 'none') else float(args.z_slice),
        logy = True if args.logy.lower().startswith('t') else False,
        ylim = None if (args.ylim[0].lower() == 'none') else args.ylim,
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
        figsizepx = [int(i) for i in args.figsizepx],
        ext = args.ext,
        figbgcolor = None if (args.figbgcolor.lower() == 'none') else args.figbgcolor,
    )
    # going home
    return s_pathfile
