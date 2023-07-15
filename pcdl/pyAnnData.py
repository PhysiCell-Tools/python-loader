####
# title: pyAnnData.py
#
# language: python3
# date: 2023-06-24
# license: BSD-3-Clause
# author: Elmar Bucher
#
# description:
#     pyAnnData.py spices up the pyMCDS (TimeStep) and pyMCDSts (TimeSeries)
#     classes with a function to transform mcds time steps and time series
#     into anndata objects.
#     anndata is the de facto python3 single cell data standard.
#     In other words, these functions enable us to analyze PhysiCell output
#     the same way that bioinformatician analyze their data retrieved from
#     single cell wet lab experiments.
#
# + https://www.biorxiv.org/content/10.1101/2021.12.16.473007v1
# + https://anndata.readthedocs.io/en/latest/
# + https://scverse.org/
####


import anndata as ad
import numpy as np
import pandas as pd
from pcdl.pyMCDS import pyMCDS, es_coor_cell
from pcdl.pyMCDSts import pyMCDSts


def scaler(df_x, scale='maxabs'):
    """
    input:
        df_x: pandas dataframe
              one feature per column, one sample per row.

        scale: string; default 'maxabs'
            None: no scaling. set scale to None if you would like to have raw data
                or entirely scale, transform, and normalize the data later.

            maxabs: maximum absolute value distance scaler will linearly map
                all values into a [-1, 1] interval. if the original data
                has no negative values, the result will be the same as with
                the minmax scaler (except with features with only one state).
                if the feature has only zeros, the value will be set to 0.

            minmax: minimum maximum distance scaler will map all values
                linearly into a [0, 1] interval.
                if the feature has only one state, the value will be set to 0.

            std: standard deviation scaler will result in sigmas.
                each feature will be mean centered around 0.
                ddof delta degree of freedom is set to 1 because it is assumed
                that the values are samples out of the population
                and not the entire population. it is incomprehensible to me
                that the equivalent sklearn method has ddof set to 0.
                if the feature has only one state, the value will be set to 0.

    output:
        df_x: pandas dataframe
            scaled df_x dataframe.

    description:
        inspired by scikit-learn's preprocessing scaling method, this function
        offers a re-implementation of the linear re-scaling methods maxabs,
        minmax, and scale.

        the robust scaler methods (quantile based) found there are missing.
        since we deal with simulated data, we don't expect heavy outliers,
        and if they exist, then they are of interest.
        the power and quantile based transformation methods and unit circle
        based normalizer methods found there are missing too.
        if you need to apply any such methods, you can do so to an anndata object
        like this:

        from sklearn import preprocessing
        adata.obsm["X_scaled"] = preprocessing.scale(adata.X)

    + https://scikit-learn.org/stable/auto_examples/preprocessing/plot_all_scaling.html
    + https://scikit-learn.org/stable/modules/classes.html#module-sklearn.preprocessing
    + https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.maxabs_scale.html
    + https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.minmax_scale.html
    + https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.scale.html
    """
    if scale is None:
        pass
    # -1,1
    elif scale == 'maxabs':
        a_x = df_x.values
        a_maxabs = a_x / abs(a_x).max(axis=0)
        a_maxabs[np.isnan(a_maxabs)] = 0  # fix if entier column is 0
        df_x = pd.DataFrame(a_maxabs, columns=df_x.columns, index=df_x.index)
    # 0,1
    elif scale == 'minmax':
        a_x = df_x.values
        a_minmax = (a_x - a_x.min(axis=0)) / (a_x.max(axis=0) - a_x.min(axis=0))
        a_minmax[np.isnan(a_minmax)] = 0  # fix if entier column has same value
        df_x = pd.DataFrame(a_minmax, columns=df_x.columns, index=df_x.index)
    # sigma
    elif scale == 'std':
        a_x = df_x.values
        a_std = (a_x - a_x.mean(axis=0)) / a_x.std(axis=0, ddof=1)
        a_std[np.isnan(a_std)] = 0  # fix if entier column has same value
        df_x = pd.DataFrame(a_std, columns=df_x.columns, index=df_x.index)
    else:
        raise ValueError(f"Error @ scaler : unknown scale algorithm {scale} detected. known are [None, 'maxabs', 'minmax', 'std'].")
    return df_x


def _anndextract(df_cell, scale='maxabs'):
    """
    input:
        df_cell:  pandas dataframe
            data frame retrieved with the mcds.get_cell_df function.

        scale: string; default 'maxabs'
            specify how the data should be scaled.
            possible values are None, maxabs, minmax, std.
            for more input, check out: help(pcdl.scaler).

    output:
        df_count, df_obs, df_spatial pandas dataframes
        ready to be backed into an anndata object.

    description:
        this function takes a pcdl df_cell pandas dataframe and re-formats
        it into a set of three dataframes (df_count, df_obs, and df_spatial),
        which downstream might be transformed into an anndata object.
    """
    # extract discrete cell data
    df_cell.drop({
        'voxel_i', 'voxel_j', 'voxel_k',
        'mesh_center_m', 'mesh_center_n', 'mesh_center_p'
    }, axis=1, inplace=True)  # maybe obs?
    df_cell.index = df_cell.index.astype(str)

    # dectect variable types
    des_type = {'float': set(), 'int': set(), 'bool': set(), 'str': set()}
    for _, se_cell in df_cell.items():
        if str(se_cell.dtype).startswith('float'):
            des_type['float'].add(se_cell.name)
        elif str(se_cell.dtype).startswith('int'):
            des_type['int'].add(se_cell.name)
        elif str(se_cell.dtype).startswith('bool'):
            des_type['bool'].add(se_cell.name)
        elif str(se_cell.dtype).startswith('object'):
            des_type['str'].add(se_cell.name)
        else:
            print(f'Error @ TimeSeries.get_anndata : column {se_cell.name} detected with unknown dtype {str(se_cell.dtype)}.')

    # build anndata object
    df_spatial = df_cell.loc[:,['position_x', 'position_y','position_z','time']].copy()
    df_obs = df_cell.loc[:,sorted(des_type['str'])].copy()
    es_num = des_type['float'].union(des_type['int'].union(des_type['bool']))
    es_num = es_num.difference(set(df_spatial.columns))
    df_count = df_cell.loc[:,sorted(es_num)].copy()
    for s_col in des_type['bool']:
        df_count[s_col] = df_count[s_col].astype(int)
    df_count = scaler(df_count, scale=scale)

    # output
    return(df_count, df_obs, df_spatial)


# class definition
class TimeStep(pyMCDS):
    """
    input:
        xmlfile: string
            name of the xml file with or without path.
            in the with path case, output_path has to be set to the default!

        output_path: string; default '.'
            relative or absolute path to the directory where
            the PhysiCell output files are stored.

        custom_type: dictionary; default is {}
            variable to specify custom_data variable types
            besides float (int, bool, str) like this: {var: dtype, ...}.
            downstream float and int will be handled as numeric,
            bool as Boolean, and str as categorical data.

        microenv: boole; default True
            should the microenvironment be extracted?
            setting microenv to False will use less memory and speed up
            processing, similar to the original pyMCDS_cells.py script.

        graph: boole; default True
            should the graphs be extracted?
            setting graph to False will use less memory and speed up processing.

        settingxml: string; default PhysiCell_settings.xml
            from which settings.xml should the substrate and cell type
            ID label mapping be extracted?
            set to None or False if the xml file is missing!

        verbose: boole; default True
            setting verbose to False for less text output while processing.

    output:
        mcds: TimeStep class instance
            all fetched content is stored at mcds.data.

    description:
        TimeStep.__init__ will call pyMCDS.__init__ that generates a mcds
        class instance, a dictionary of dictionaries data structure that
        contains all output from a single PhysiCell model time step.
        furthermore, the mcds object offers functions to access the stored data.
        the code assumes that all related output files are stored
        in the same directory. data is loaded by reading the xml file for
        a particular time step and the therein referenced files.
    """
    def __init__(self, xmlfile, output_path='.', custom_type={}, microenv=True, graph=True, settingxml='PhysiCell_settings.xml', verbose=True):
        pyMCDS.__init__(self, xmlfile=xmlfile, output_path=output_path, custom_type=custom_type, microenv=microenv, graph=graph, settingxml=settingxml, verbose=verbose)

    def get_anndata(self, states=1, drop=set(), keep=set(), scale='maxabs'):
        """
        input:
            states: integer; default is 1
                minimal number of states a variable has to have to be outputted.
                variables that have only 1 state carry no information.
                None is a state too.

            drop: set of strings; default is an empty set
                set of column labels to be dropped for the dataframe.
                don't worry: essential columns like ID, coordinates
                and time will never be dropped.
                Attention: when the keep parameter is given, then
                the drop parameter has to be an empty set!

            keep: set of strings; default is an empty set
                set of column labels to be kept in the dataframe.
                set states=1 to be sure that all variables are kept.
                don't worry: essential columns like ID, coordinates
                and time will always be kept.

            scale: string; default 'maxabs'
                specify how the data should be scaled.
                possible values are None, maxabs, minmax, std.
                for more input, check out: help(pcdl.scaler)

        output:
            anmcds: anndata object
                for this one time step.

        description:
            function to transform a mcds time step into an anndata object
            for downstream analysis.
        """
        # processing
        print(f'processing: 1/1 {self.get_time()}[min] mcds into anndata obj.')
        df_cell = self.get_cell_df(states=states, drop=drop, keep=keep)
        df_count, df_obs, df_spatial = _anndextract(df_cell=df_cell, scale=scale)
        anmcds = ad.AnnData(X=df_count, obs=df_obs, obsm={"spatial": df_spatial.values})

        # output
        return anmcds


class TimeSeries(pyMCDSts):
    """
    input:
        output_path: string, default '.'
            relative or absolute path to the directory where
            the PhysiCell output files are stored.

        custom_type: dictionary; default is {}
            variable to specify custom_data variable types
            besides float (int, bool, str) like this: {var: dtype, ...}.
            downstream float and int will be handled as numeric,
            bool as Boolean, and str as categorical data.

        load: boole; default True
            should the whole time series data, all time steps, straight at
            object initialization be read and stored to mcdsts.l_mcds?

        microenv: boole; default True
            should the microenvironment be extracted?
            setting microenv to False will use less memory and speed up
            processing, similar to the original pyMCDS_cells.py script.

        graph: boole; default True
            should the graphs be extracted?
            setting graph to False will use less memory and speed up processing.

        settingxml: string; default PhysiCell_settings.xml
            from which settings.xml should the substrate and cell type
            ID label mapping be extracted?
            set to None or False if the xml file is missing!

        verbose: boole; default True
            setting verbose to False for less text output while processing.

    output:
        mcdsts: pyMCDSts class instance
            this instance offers functions to process all stored time steps
            from a simulation.

    description:
        TimeSeries.__init__ will call pyMCDSts.__init__ that generates a mcdsts
        class instance. this instance offers functions to process all time steps
        in the output_path directory.
    """
    def __init__(self, output_path='.', custom_type={}, load=True, microenv=True, graph=True, settingxml='PhysiCell_settings.xml', verbose=True):
        pyMCDSts.__init__(self, output_path=output_path, custom_type=custom_type, load=load, microenv=microenv, graph=graph, settingxml=settingxml, verbose=verbose)

    def get_anndata(self, states=1, drop=set(), keep=set(), scale='maxabs', collapse=True, keep_mcds=True):
        """
        input:
            states: integer; default is 1
                minimal number of states a variable has to have to be outputted.
                variables that have only 1 state carry no information.
                None is a state too.

            drop: set of strings; default is an empty set
                set of column labels to be dropped for the dataframe.
                don't worry: essential columns like ID, coordinates
                and time will never be dropped.
                Attention: when the keep parameter is given, then
                the drop parameter has to be an empty set!

            keep: set of strings; default is an empty set
                set of column labels to be kept in the dataframe.
                don't worry: essential columns like ID, coordinates
                and time will always be kept.

            scale: string; default 'maxabs'
                specify how the data should be scaled.
                possible values are None, maxabs, minmax, std.
                for more input, check out: help(pcdl.scaler)

            collapse: boole; default True
                should all mcds time steps from the time series be collapsed
                into one single anndata object, or a list of anndata objects
                for each time step?

            keep_mcds: boole; default True
                should the loaded original mcds be kept in memory
                after transformation?

        output:
            anmcds or l_anmcds: anndata object or list of anndata objects.
                what is returned depends on the collapse setting.

        description:
            function to transform mcds time steps into one or many
            anndata objects for downstream analysis.
        """
        d_ann = {}
        df_anncount = None
        df_annobs = None
        df_annspatial = None

        # variable triage
        if (states < 2):
            ls_column = self.l_mcds[0].get_cell_df()
        else:
            ls_column = sorted(es_coor_cell.difference({'ID'}))
            ls_column.extend(self.get_cell_df_columns_min_states(states=states, drop=drop, keep=keep))

        # processing
        i_mcds = len(self.l_mcds)
        for i in range(i_mcds):
            # fetch mcds
            if keep_mcds:
                mcds = self.l_mcds[i]
            else:
                mcds = self.l_mcds.pop(0)
            # extract time and dataframes
            i_time = int(round(mcds.get_time()))
            print(f'processing: {i+1}/{i_mcds} {i_time}[min] mcds into anndata obj.')
            df_cell = mcds.get_cell_df()
            df_cell = df_cell.loc[:,ls_column]
            df_count, df_obs, df_spatial = _anndextract(df_cell=df_cell, scale=scale)
            # pack collapsed
            if collapse:
                # count
                df_count.reset_index(inplace=True)
                df_count.index = df_count.ID + f'id_{i_time}min'
                df_count.index.name = 'ID_time'
                df_count.drop('ID', axis=1, inplace=True)
                if df_anncount is None:
                    df_anncount = df_count
                else:
                    df_anncount = pd.concat([df_anncount, df_count], axis=0)
                # obs
                df_obs.reset_index(inplace=True)
                df_obs.index = df_obs.ID + f'id_{i_time}min'
                df_obs.index.name = 'ID_time'
                if df_annobs is None:
                    df_annobs = df_obs
                else:
                    df_annobs = pd.concat([df_annobs, df_obs], axis=0)
                # spatial
                df_spatial.reset_index(inplace=True)
                df_spatial.index = df_spatial.ID + f'id_{i_time}min'
                df_spatial.index.name = 'ID_time'
                df_spatial.drop('ID', axis=1, inplace=True)
                if df_annspatial is None:
                    df_annspatial = df_spatial
                else:
                    df_annspatial = pd.concat([df_annspatial, df_spatial], axis=0)
            # pack not collapsed
            else:
                annmcds = ad.AnnData(X=df_count, obs=df_obs, obsm={"spatial": df_spatial.values})
                d_ann.update({i_time : annmcds})

        # output
        if collapse:
            annmcdsts = ad.AnnData(X=df_anncount, obs=df_annobs, obsm={"spatial": df_annspatial.values})
        else:
            annmcdsts = d_ann
        return annmcdsts

