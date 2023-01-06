import anndata as ad
import pandas as pd
#import pcDataLoader as pc


# function
def annD(xml_file, output_path='.', microenv=True, categorical={}):  # numerical, ordinal
    '''
    '''
    # test case
    xml_file = 'output00003696.xml'
    output_path = 'src_github/pcDataLoader/pcDataLoader/data_snapshot/'
    microenv=True

    # definition
    es_spatial = {
        'position_x',
        'position_y',
        'position_z',
    }
    es_voxel = {
        'voxel_i',
        'voxel_j',
        'voxel_k',
        'voxel_position_m',
        'voxel_position_n',
        'voxel_position_o',
    }
    es_obs = {
        'cell_type',  # p26 type
        'current_phase',
        'cycle_model',
        'voxel_i',
        'voxel_j',
        'voxel_k',
        'voxel_position_m',
        'voxel_position_n',
        'voxel_position_o'
    }

    es_motility  = {
        'migration_bias',  # could this have biological information?
        'migration_bias_direction',  # will not have biological information
        #'migration_speed',  # could this have biological information?
    }

    # get data
    #mcds = pc.pyMCDS(xml_file=xml_file, output_path=output_path, microenv=microenv)

    # extract discrete cell data
    df_cell = mcds.get_cell_df()

    # extract continuous microenv data
    if microenv:
        df_cell.reset_index(inplace=True)
        df_conc = mcds.get_concentrations_df()
        df_extract = pd.merge(df_cell, df_conc, on=[sorted(es_voxel)], how='left')
        df_extract.set_index('ID', inplace=True)
        df_cell.set_index('ID', inplace=True)
    else:
        df_count = df_cell

    # build anndata object
    df_count.index = df_count.index.astype(str)
    df_spatial = df_count.loc[:,[sorted(es_spatial)]]
    df_obs = df_count.loc[:,[sorted(es_obs)]]
    es_count = set(df_count.columns).difference(es_spatial.union(es_obs))
    df_count = df_extract.loc[:, sorted(es_count)]
    adata = ad.AnnData(X=df_count, obs=df_obs, obsm={"spatial": df_spatial.values})

    # output
    return adata
