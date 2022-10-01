import anndata as ad
import pandas as pd
import pcDataLoader as pc


# function
def AnnData(xml_file, output_path='.', microenv=True):
    '''
    '''
    xml_file = 'output00003696.xml' 
    output_path = 'src_github/pcDataLoader/pcDataLoader/data_snapshot/'
    microenv=True

    # get data
    mcds = pc.pyMCDS(xml_file=xml_file, output_path=output_path, microenv=microenv)

    # extract discrete cell data
    df_cell = mcds.get_cell_df()

    # extract contiuous microenv data
    if microenv:
        df_cell.reset_index(inplace=True)
        df_conc = mcds.get_concentrations_df()
        df_count = pd.merge(df_cell, df_conc, on=['voxel_i','voxel_j','voxel_k'], how='left')
        df_count.set_index('ID', inplace=True)
        df_cell.set_index('ID', inplace=True)
    else:
        df_count = df_cell

    # build anndata object
    df_count.index = df_count.index.astype(str)
    df_spatial = df_count.loc[:,['position_x', 'position_y','position_z']]
    df_obs = df_count.loc[:,[
        'cell_type','current_phase','cycle_model',
        'voxel_i','voxel_j','voxel_k', 
        'voxel_position_m','voxel_position_n','voxel_position_o'
    ]]
    es_count = set(df_count.columns).difference({
        'cell_type','current_phase','cycle_model',
        'voxel_i', 'voxel_j', 'voxel_k', 
        'voxel_position_m', 'voxel_position_n', 'voxel_position_o',
        'position_x', 'position_y','position_z'
    })
    df_count = df_count.loc[:, sorted(es_count)]
    adata = ad.AnnData(X=df_count, obs=df_obs, obsm={"spatial": df_spatial.values})

    # output
    return adata
