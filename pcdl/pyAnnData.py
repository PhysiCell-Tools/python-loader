###
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
from scipy import sparse
import warnings


def scaler(df_x, scale='maxabs'):
    """
    input:
        df_x: pandas dataframe
              one attribute per column, one sample per row.

        scale: string; default 'maxabs'
            None: no scaling. set scale to None if you would like to have
                raw data or scale, transform, and normalize the data later.

            maxabs: maximum absolute value distance scaler will linearly map
                all values into a [-1, 1] interval. if the original data
                has no negative values, the result will be the same as with
                the minmax scaler (except with attributes with only one value).
                if the attribute has only zeros, the value will be set to 0.

            minmax: minimum maximum distance scaler will map all values
                linearly into a [0, 1] interval.
                if the attribute has only one value, the value will be set to 0.

            std: standard deviation scaler will result in sigmas.
                each attribute will be mean centered around 0.
                ddof delta degree of freedom is set to 1 because it is assumed
                that the values are samples out of the population
                and not the entire population. it is incomprehensible to me
                that the equivalent sklearn method has ddof set to 0.
                if the attribute has only one value, the value will be set to 0.

    output:
        df_x: pandas dataframe
            scaled df_x dataframe.

    description:
        inspired by scikit-learn's preprocessing scaling method, this function
        offers a re-implementation of the linear re-scaling methods maxabs,
        minmax, and scale.

        the robust scaler methods (quantile based) found in scikit-learn are
        missing. since we deal with simulated data, we don't expect heavy
        outliers, and if they exist, then they are of interest.
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
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        a_maxabs = a_x / abs(a_x).max(axis=0)
        warnings.simplefilter('default')
        a_maxabs[np.isnan(a_maxabs)] = 0  # fix if entier column is 0
        df_x = pd.DataFrame(a_maxabs, columns=df_x.columns, index=df_x.index)
    # 0,1
    elif scale == 'minmax':
        a_x = df_x.values
        warnings.simplefilter("ignore")
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        a_minmax = (a_x - a_x.min(axis=0)) / (a_x.max(axis=0) - a_x.min(axis=0))
        warnings.simplefilter('default')
        a_minmax[np.isnan(a_minmax)] = 0  # fix if entier column has same value
        df_x = pd.DataFrame(a_minmax, columns=df_x.columns, index=df_x.index)
    # sigma
    elif scale == 'std':
        a_x = df_x.values
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        a_std = (a_x - a_x.mean(axis=0)) / a_x.std(axis=0, ddof=1)
        warnings.simplefilter('default')
        a_std[np.isnan(a_std)] = 0  # fix if entier column has same value
        df_x = pd.DataFrame(a_std, columns=df_x.columns, index=df_x.index)
    else:
        raise ValueError(f"Error @ scaler : unknown scale algorithm {scale} detected. known are [None, 'maxabs', 'minmax', 'std'].")

    return df_x


def _anndextract(df_cell, scale='maxabs', graph_attached={}, graph_neighbor={}, graph_spring={}, graph_method='PhysiCell'):
    """
    input:
        df_cell:  pandas dataframe
            data frame retrieved with the mcds.get_cell_df function.

        scale: string; default maxabs
            specify how the data should be scaled.
            possible values are None, maxabs, minmax, std.
            for more input, check out: help(pcdl.scaler).

        graph_attached: dict; default {}
            attached graph dictionary, retrieved with
            with the mcds.get_attched_graph() function.

        graph_neighbor: dict; default {}
            neighbor graph dictionary, retrieved
            with the mcds.get_neighbor_graph() function.

        graph_spring: dict; default {}
            spring_attached graph dictionary, retrieved
            with the mcds.get_spring_graph_dict() function.

        graph_method: string; default PhysiCell
            method how the graphs were generated.

    output:
        df_count, df_obs, d_obsm, d_obsp, d_uns dataframes and dictionaries,
            ready to be backed into an anndata object.

    description:
        this function takes a pcdl df_cell pandas dataframe and re-formats
        it into a set of two dataframes (df_count, df_obs),
        two dictionary of numpy array (d_obsm, d_obsp),
        and one dictionary of string (d_uns),
        which downstream might be transformed into an anndata object.
    """
    # transform index to string
    df_coor = df_cell.loc[:,['position_x','position_y','position_z']].copy()
    df_cell.index = df_cell.index.astype(str)

    # build obs anndata object (annotation of observations)
    df_obs = df_cell.loc[:,['mesh_center_p','time']].copy()
    df_obs.columns = ['z_layer', 'time']

    # buil obsm anndata object spatial (multi-dimensional annotation of observations)
    if (len(set(df_cell.position_z)) == 1):
        df_obsm = df_cell.loc[:,['position_x','position_y']].copy()
    else:
        df_obsm = df_cell.loc[:,['position_x','position_y','position_z']].copy()
    d_obsm = {"spatial": df_obsm.values}

    # build obsp and uns anndata object graph (pairwise annotation of obeservation) and (unstructured data)
    ####
    # acknowledgement:
    #   this code is inspired from the tysserand add_to_AnnData impelmentation
    #   from Alexis Coullomb form the Pancaldi Lab.
    #   https://github.com/VeraPancaldiLab/tysserand/blob/main/tysserand/tysserand.py#L1546
    ####
    # extract cell_id to index mapping (i always loved perl)
    di_ididx = df_cell.reset_index().loc[:,'ID'].reset_index().astype(int).set_index('ID').squeeze().to_dict()
    # transform cell id graph dict to index matrix and pack for anndata
    d_obsp = {}  # pairwise annotation of obeservation
    d_uns = {}  # unstructured data
    for s_graph, dei_graph in [('neighbor', graph_neighbor), ('attached', graph_attached), ('spring', graph_spring)]:
        lli_edge = []
        lr_distance = []
        for i_src, ei_dst in dei_graph.items():
            for i_dst in ei_dst:
                # extract edge
                lli_edge.append([di_ididx[i_src], di_ididx[i_dst]])
                r_distance = ((df_coor.loc[i_src,:].values -  df_coor.loc[i_dst,:].values)**2).sum()**(1/2)
                lr_distance.append(r_distance)
        # if there is a graph
        if (len(lli_edge) > 0):
            # handle edge data
            ai_edge = np.array(lli_edge, dtype=np.uint)
            # handle connection data
            ai_conectivity = np.ones(ai_edge.shape[0], dtype=np.uint16)
            ai_conectivity_sparse = sparse.csr_matrix(
                (ai_conectivity, (ai_edge[:,0], ai_edge[:,1])),
                shape = (df_cell.shape[0], df_cell.shape[0]),
                dtype = np.uint
            )
            # handle distance data
            ar_distance  = np.array(lr_distance, dtype=np.float64)
            ar_distance_sparse = sparse.csr_matrix(
                (ar_distance, (ai_edge[:,0], ai_edge[:,1])),
                shape = (df_cell.shape[0], df_cell.shape[0]),
                dtype = np.float64
            )
            # pack obsp
            d_obsp.update({
                f'physicell_{s_graph}_conectivities': ai_conectivity_sparse,
                f'physicell_{s_graph}_distances': ar_distance_sparse,
            })
            # pack uns
            d_uns.update({
                s_graph : {
                    'connectivities_key': f'physicell_{s_graph}_conectivities',
                    'distances_key': f'physicell_{s_graph}_distances',
                    'params': {
                        'metric': 'euclidean',
                        'method': graph_method,
                    }
                }
            })

    # extract discrete cell data
    es_drop = set(df_cell.columns).intersection({
        'voxel_i', 'voxel_j', 'voxel_k',
        'mesh_center_m', 'mesh_center_n', 'mesh_center_p',
        'position_x', 'position_y','position_z',
        'time', 'runtime', 'xmlfile',
    })
    df_cell.drop(es_drop, axis=1, inplace=True)  # maybe obs?

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

    # build on obs and X anndata object
    df_cat = df_cell.loc[:,sorted(des_type['str'])].copy()
    df_obs = pd.merge(df_obs, df_cat, left_index=True, right_index=True)
    es_num = des_type['float'].union(des_type['int'].union(des_type['bool']))
    df_count = df_cell.loc[:,sorted(es_num)].copy()
    for s_col in des_type['bool']:
        df_count[s_col] = df_count[s_col].astype(int)
    df_count = scaler(df_count, scale=scale)

    # return
    return(df_count, df_obs, d_obsm, d_obsp, d_uns)


# class definition
class TimeStep(pyMCDS):
    def __init__(self, xmlfile, output_path='.', custom_data_type={}, microenv=True, graph=True, physiboss=True, settingxml='PhysiCell_settings.xml', verbose=True):
        """
        input:
            xmlfile: string
                name of the xml file with or without path.
                in the with path case, output_path has to be set to the default!

            output_path: string; default '.'
                relative or absolute path to the directory where
                the PhysiCell output files are stored.

            custom_data_type: dictionary; default is {}
                variable to specify custom_data variable types
                besides float (int, bool, str) like this: {var: dtype, ...}.
                downstream float and int will be handled as numeric,
                bool as Boolean, and str as categorical data.

            microenv: boole; default True
                should the microenvironment data be loaded?
                setting microenv to False will use less memory and speed up
                processing, similar to the original pyMCDS_cells.py script.

            graph: boole; default True
                should neighbor garph, attached graph, and spring attached graph
                be loaded? setting graph to False will use less memory and
                speed up processing.

            physiboss: boole; default True
                should physiboss state data be loaded, if found?
                setting physiboss to False will use less memory and speed up processing.

            settingxml: string; default PhysiCell_settings.xml
                the settings.xml that is loaded, from which the cell type ID
                label mapping, is extracted, if this information is not found
                in the output xml file.
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
        pyMCDS.__init__(
            self,
            xmlfile = xmlfile,
            output_path = output_path,
            custom_data_type = custom_data_type,
            microenv = microenv,
            graph = graph,
            physiboss = physiboss,
            settingxml = settingxml,
            verbose = verbose
        )


    def get_anndata(self, values=1, drop=set(), keep=set(), scale='maxabs'):
        """
        input:
            values: integer; default is 1
                minimal number of values a variable has to have to be outputted.
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
                set values=1 to be sure that all variables are kept.
                don't worry: essential columns like ID, coordinates
                and time will always be kept.

            scale: string; default 'maxabs'
                specify how the data should be scaled.
                possible values are None, maxabs, minmax, std.
                for more input, check out: help(pcdl.scaler)

        output:
            annmcds: anndata object
                for this one time step.

        description:
            function to transform a mcds time step into an anndata object
            for downstream analysis.
        """
        # processing
        if self.verbose:
            print(f'processing: 1/1 {round(self.get_time(),9)}[min] mcds into anndata obj.')
        df_cell = self.get_cell_df(values=values, drop=drop, keep=keep)
        df_count, df_obs, d_obsm, d_obsp, d_uns = _anndextract(
            df_cell = df_cell,
            scale = scale,
            graph_attached = self.get_attached_graph_dict(),
            graph_neighbor = self.get_neighbor_graph_dict(),
            graph_method = self.get_physicell_version(),
        )
        annmcds = ad.AnnData(
            X = df_count,
            obs = df_obs,
            obsm = d_obsm,
            obsp = d_obsp,
            uns = d_uns
        )
        # output
        return annmcds


class TimeSeries(pyMCDSts):
    def __init__(self, output_path='.', custom_data_type={}, load=True, microenv=True, graph=True, physiboss=True, settingxml='PhysiCell_settings.xml', verbose=True):
        """
        input:
            output_path: string, default '.'
                relative or absolute path to the directory where
                the PhysiCell output files are stored.

            custom_data_type: dictionary; default is {}
                variable to specify custom_data variable types
                besides float (int, bool, str) like this: {var: dtype, ...}.
                downstream float and int will be handled as numeric,
                bool as Boolean, and str as categorical data.

            load: boole; default True
                should the whole time series data, all time steps, straight at
                object initialization be read and stored to mcdsts.l_mcds?

            microenv: boole; default True
                should the microenvironment data be loaded?
                setting microenv to False will use less memory and speed up
                processing, similar to the original pyMCDS_cells.py script.

            graph: boole; default True
                should neighbor garph, attached graph, and spring attached graph
                be loaded? setting graph to False will use less memory and
                speed up processing.

            physiboss: boole; default True
                should physiboss state data be loaded, if found?
                setting physiboss to False will use less memory and speed up processing.

            settingxml: string; default PhysiCell_settings.xml
                the settings.xml that is loaded, from which the cell type ID
                label mapping, is extracted, if this information is not found
                in the output xml file.
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
        pyMCDSts.__init__(
            self,
            output_path = output_path,
            custom_data_type = custom_data_type,
            load = load,
            microenv = microenv,
            graph = graph,
            physiboss = physiboss,
            settingxml = settingxml,
            verbose = verbose
        )
        self.l_annmcds = None


    def get_anndata(self, values=1, drop=set(), keep=set(), scale='maxabs', collapse=True, keep_mcds=True):
        """
        input:
            values: integer; default is 1
                minimal number of values a variable has to have to be outputted.
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
            annmcds or self.l_annmcds: anndata object or list of anndata objects.
                what is returned depends on the collapse setting.

        description:
            function to transform mcds time steps into one or many
            anndata objects for downstream analysis.
        """
        # initialize vaiable
        l_annmcds = []
        df_anncount = None
        df_annobs = None
        ar_annobsm = None

        # variable triage
        if (values < 2):
            ls_column = list(self.l_mcds[0].get_cell_df(drop=drop, keep=keep).columns)
        else:
            ls_column = sorted(es_coor_cell.difference({'ID'}))
            ls_column.extend(sorted(self.get_cell_attribute(values=values, drop=drop, keep=keep, allvalues=False).keys()))

        # collapse warning
        if collapse and self.verbose:
            print('Warning @ mcdsts.get_anndata : only df_cell data, but not graph data, can be collapsed.')

        # processing
        lann_mcds = []
        i_mcds = len(self.l_mcds)
        for i in range(i_mcds):
            # fetch mcds
            if keep_mcds:
                mcds = self.l_mcds[i]
            else:
                mcds = self.l_mcds.pop(0)
            # extract physicell version
            s_physicellv = mcds.get_physicell_version(),
            # extract time and dataframes
            r_time = round(mcds.get_time(),9)
            if self.verbose:
                print(f'processing: {i+1}/{i_mcds} {r_time}[min] mcds into anndata obj.')
            df_cell = mcds.get_cell_df()
            df_cell = df_cell.loc[:,ls_column]

            # pack collapsed
            if collapse:
                # extract
                df_count, df_obs, d_obsm, d_obsp, d_uns = _anndextract(
                    df_cell=df_cell,
                    scale = scale,
                    #graph_attached = {},
                    #graph_neighbor = {},
                    #graph_spring = {},
                    #graph_method = s_physicellv,
                )
                # count
                df_count.reset_index(inplace=True)
                df_count.index = df_count.ID + f'id_{r_time}min'
                df_count.index.name = 'id_time'
                df_count.drop('ID', axis=1, inplace=True)
                if df_anncount is None:
                    df_anncount = df_count
                else:
                    df_anncount = pd.concat([df_anncount, df_count], axis=0)
                # obs
                df_obs.reset_index(inplace=True)
                df_obs.index = df_obs.ID + f'id_{r_time}min'
                df_obs.index.name = 'id_time'
                if df_annobs is None:
                    df_annobs = df_obs
                else:
                    df_annobs = pd.concat([df_annobs, df_obs], axis=0)
                # obsm (spatial)
                if ar_annobsm is None:
                    ar_annobsm = d_obsm['spatial']
                else:
                    ar_annobsm = np.vstack([ar_annobsm, d_obsm['spatial']])
                # obsp: nop (graph)
                # uns: nop (graph)

            # pack not collapsed
            else:
                # extract
                df_count, df_obs, d_obsm, d_obsp, d_uns = _anndextract(
                    df_cell=df_cell,
                    scale = scale,
                    graph_attached = mcds.get_attached_graph_dict(),
                    graph_neighbor = mcds.get_neighbor_graph_dict(),
                    graph_spring = mcds.get_spring_graph_dict(),
                    graph_method = s_physicellv,
                )
                # annmcds
                ann_mcds = ad.AnnData(
                    X = df_count,
                    obs = df_obs,
                    obsm = d_obsm,
                    obsp = d_obsp,
                    uns = d_uns,
                )
                lann_mcds.append(ann_mcds)

        # output
        if collapse:
            ann_mcdsts = ad.AnnData(
                X = df_anncount,
                obs = df_annobs,
                obsm = {'spatial': ar_annobsm},
                #obsp = d_obsp,
                #uns = d_uns
            )
            return ann_mcdsts
        else:
            self.l_annmcds = lann_mcds
            return self.l_annmcds


    def get_annmcds_list(self):
        """
        input:
            self: TimeSeries class instance.

        output:
            self.l_annmcds: list of chronologically ordered anndata mcds objects.
                watch out, this is a pointer to the
                self.l_annmcds list of anndata mcds objects, not a copy of self.l_annmcds!

        description:
            function returns a binding to the self.l_annmcds list of anndata mcds objects.
        """
        return self.l_annmcds

