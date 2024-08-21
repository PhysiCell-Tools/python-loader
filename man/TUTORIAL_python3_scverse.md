# Working with PhysiCell Data in Python

+ author: Elmar Bucher
+ data: July 2023

#### description:
+ this notebook should give you an idea about how to work with pcdl in a Python3 REPL (read eval print loop) to analyse PhysiCell data.

#### installation instruction:
1. copy this Jupyter notebook into your PhysiCell folder.
2. open the notebook in Jupyter or JupterLab.

# A) Preparation


```python
# python library installation
!pip3 install -U pcdl
!pip3 install -U scanpy[leiden]  # single cell analysis inclusive leiden graph clustering algorithm.
!pip3 install -U squidpy # single cell spatial analysis.
```


```python
# compile and run the PhysiCell 2D interaction-sample project
# note: this will easily take 10[min].

# uncomment commands below to run!
#!make reset
#!make data-clean
#!make list-projects
#!make interaction-sample
#!make
#!if [ ! -d output ]; then mkdir -p output; fi;
#!./interaction_demo
#!mv output output2d
```


```python
# compile and run the first 24 hours of the PhysiCell 3D cancer-immune-sample project
# note: this will easily take 70[min].

# uncomment commands below to run!
#!make reset
#!make data-clean
#!make list-projects
#!make cancer-immune-sample
#!make
#
#!if [ ! -d output ]; then mkdir -p output; fi;
#import xml.etree.ElementTree as ET
#x_tree = ET.parse('config/PhysiCell_settings.xml')
#x_root = x_tree.getroot()
#x_element = x_root.find('.//overall/max_time')
#x_element.text =  '1440'  # = 24[h] * 60[min/h]
#x_tree.write(
#    'config/PhysiCell_settings.xml',
#    xml_declaration='<?xml version="1.0" encoding="UTF-8"?>'
#)
#
#!./cancer_immune_3D
#!mv output output3d
```

# B) Interactive PhysiCell output data analysis


```python
# library
import anndata as ad  # from the scverse
import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import pcdl
import scanpy as sc  # from the scverse
import squidpy as sq  # from the scverse

# const
s_path_2d = 'output2d/'
s_pathfile_2d = f'{s_path_2d}output00000009.xml'   # 24[h] = 1440[min]
s_pathfile_3d = 'output3d/output00000024.xml'  # 24[h] = 1440[min]

# pcdl version
pcdl.__version__
```

## 1 pcdl Time Series - chronological list of mcds


```python
###################
# LOAD TimeSeries #
###################

#mcdsts = pcdl.TimeSeries?
mcdsts = pcdl.TimeSeries(s_path_2d)

# note to the Warning:
# these custom_data attributes are all numeric (not categorical),
# the default float variable type is ok (if needed, other possible types are: int, bool, str).
```

### 1.1 time series: plot_scatter, plot_contour, make_gif, make_movie


```python
# make images and movies

#help(mcdsts.plot_scatter)
#mcdsts.make_movie?
#mcdsts.make_gif?

#s_pathfile = mcdsts.plot_scatter()
#print('path file:', s_pathfile)
#mcdsts.make_movie(s_pathfile)
#mcdsts.make_gif(s_pathfile)

mcdsts.make_gif(mcdsts.plot_scatter())  # cmap='turbo'
```


```python
# plot_scatter and plot_contour difference

#mcdsts.get_mcds_list()[0].get_substrate_list?
#mcdsts.get_mcds_list()[0].get_conc_df().columns?

#mcdsts.get_mcds_list()[0].get_celltype_list?
#mcdsts.get_mcds_list()[0].get_cell_df().columns?

for s_subs in mcdsts.get_mcds_list()[0].get_substrate_list():
    mcdsts.make_gif(mcdsts.plot_scatter(s_subs, cmap='turbo'))
    mcdsts.make_gif(mcdsts.plot_contour(s_subs, cmap='turbo'))
    break
```

### 1.2 time series: get_cell_attribute, get_conc_attribute


```python
#mcdsts.get_cell_attribute?

dl_list = mcdsts.get_cell_attribute(allvalues=True, values=1)  # values=2
print('total attributes df_cell:', len(dl_list))
json.dump(dl_list, open(f'{s_path_2d}cell_attribute.json', 'w'))
```


```python
#mcdsts.get_conc_attribute?

dl_list = mcdsts.get_conc_attribute(allvalues=True, values=1)  # values=2
print('total attributes df_conc:', len(dl_list))
json.dump(dl_list, open(f'{s_path_2d}conc_attribute.json', 'w'))
```

### 1.3 time series: get_anndata


```python
#mcdsts.get_anndata?
annts = mcdsts.get_anndata()
l_ann = mcdsts.get_anndata(collapse=False)
```


```python
print('annts:\n', annts)
print()
print('l_ann:\n', l_ann[0:3])
```

## 2 pcdl Time Step - one mcds


```python
#mcdsts.get_mcds_list?
#mcdsts.get_mcds_list()
mcds = mcdsts.get_mcds_list()[9]  # 24[h] = 1440[min]  # watch out: this mcds is a pyMCDS object and has no get_anndata function!
```


```python
#mcdsts.get_annmcds_list?
#mcdsts.get_annmcds_list()
mcds = mcdsts.get_annmcds_list()[9]  # 24[h] = 1440[min]  # this mcds is a TimeStep object and has a get_anndata function!
```


```python
#################
# LOAD TimeStep #
#################

# pcdl.TimeStep?
mcds = pcdl.TimeStep(f'{s_path_2d}output00000009.xml')    # 24[h] = 1440[min]  # this mcds is a TimeStep object has a get_anndata functionality!
```

### 2.1 time step: cell_df, conc_df, get_substarte_names, plot_scatter, plot_contour


```python
#df = mcds.get_cell_df?
df_cell = mcds.get_cell_df(values=2)
df_cell.info()
sorted(df_cell.columns)

# here I break with the rule that pcdl is simply an interface.
#mcds.plot_scatter?
mcds.plot_scatter()
mcds.plot_scatter('toxin')
print()
```


```python
df_conc = mcds.get_conc_df(values=2)
df_conc.info()
mcds.get_substrate_list()

# here I break with the rule that pcdl is simply an interface.
# this implementation based on matplotlib contour and contourf.
#mcds.plot_contour?
mcds.plot_contour('toxin', vmin=0.0, vmax=0.16)
print()
```

### 2.2 time step: get_unit_se


```python
#mcds.get_unit_se?
mcds.get_unit_se()
```

### 2.3 time step: get_anndata


```python
# loads only attribute that have not the same value in all cells.
# max absolute scales the attributes into a range between -1 and 1.
ann = mcds.get_anndata(values=2, scale='maxabs')
```


```python
# https://scverse.org/
# https://anndata.readthedocs.io/en/latest/
ann?
```


```python
ann
#ann.obs.head()
#ann.obsm
#ann.obsm['spatial']
#nn.var
#ann.var.keys
#ann.X
```


```python
# save and load anndata objects
ann.write(f'{s_path_2d}interaction_16200min.h5ad')
```


```python
# save and load anndata objects
ann = ad.read(f'{s_path_2d}interaction_16200min.h5ad')
```


```python
# attributes
print('x_axis: genes: numerical attributes:\n', ann.var_names)  # list the numerical attributes we have at hand (alternative way: ann.var.index).
print('y_axis: cells: categorical attributes:\n', ann.obs_keys())  # list the categories attributes we have at hand (alternative way: ann.obs.columns).
```

## 3 Data Analysis with Pandas


```python
df_cell.plot?
```

### 3.1 time step: pandas plot categorical data


```python
# bar plot
df = df_cell.loc[:,['cell_type']]   # ['cell_type','current_phase']
df['count'] = '1'
df_bar = df.groupby(['cell_type']).count()  # ['cell_type','current_phase']
df_bar
df_bar.plot(kind='bar')
df_bar.plot(kind='barh')
```


```python
# pie plot
#df_bar.plot(kind='pie')
#df_bar.plot(kind='pie', subplots=True)
#df_bar.plot(kind='pie', y='count')
df_bar.plot(kind='pie', y='count', legend=False, ylabel='', title='cell_type fraction')
```

### 3.2 time step: pandas plot numerical data


```python
# count histogram
df_cell.loc[:,'toxin'].plot(kind='hist', bins=32, title='cell count toxin cell surrounding')  # series
#df_cell.loc[:,['toxin']].plot(kind='hist', bins=32, title='cell surrounding toxin')  # datafarme
#df_cell.loc[:,['toxin','debris','resource','pro-inflammatory']].plot(kind='hist', bins=32, title='cell surrounding substrate')  # stacked=True, alpha=0.5
```


```python
# probability histogram
se_toxin = df_cell.loc[:,'toxin']
a_ones = np.ones_like(se_toxin.values)
a_weight = a_ones / se_toxin.shape[0]
se_toxin.plot(kind='hist', bins=32, title='cell fraction toxin in cell surrounding', weights=a_weight)
```


```python
# kde kernel density estimation
df_cell.loc[:,['toxin','debris','resource', 'pro-inflammatory']].plot(kind='kde', title='cell surounding substrate')
```


```python
# box plot
df_cell.loc[:,mcds.get_substrate_list()].plot(kind='box', title='cell surrounding substrate')
```

### 3.3 time step: pandas and pcdl plot numerical, spatial data


```python
# scatter plot - plot_scatter() equivalent
#df_cell.plot(kind='scatter', x='position_x', y='position_y')

## map cell_type and colors ##
#sorted(df_cell.loc[:,'cell_type'].unique())
ds_color = {
    'CD8+_T_cell' : 'magenta',
    'macrophage' : 'orange',
    'neutrophil' : 'yellow',
    'bacteria' : 'black',
    'blood_vessel' : 'red',
    'differentiated' : 'green',
    'stem' : 'lime',
}

df_cell['cell_type_color'] = None
#print('df_cell before loop:\n', df_cell.loc[:,['cell_type', 'cell_type_color']])  # show dataframe
for s_celltype in df_cell.loc[:,'cell_type'].unique():
    df_cell.loc[df_cell.loc[:,'cell_type'] == s_celltype, 'cell_type_color'] = ds_color[s_celltype]
#print('df_cell after loop:\n', df_cell.loc[:,['cell_type', 'cell_type_color']])  # show dataframe
df_cell.plot(kind='scatter', x='position_x', y='position_y', c=df_cell.loc[:,'cell_type_color'], s=9)

## add color legend ##
lo_patch = []
for s_label, s_color in sorted(ds_color.items()):
    o_patch = mpatches.Patch(color=s_color, label=s_label)
    lo_patch.append(o_patch)

ax = plt.gca()
ax.legend(
    handles = lo_patch,
    loc = 'lower left',
    fontsize = 'small'
)
```


```python
# hexbin - plot_contour() equivalent
#df_conc.plot(kind='hexbin', x='mesh_center_m', y='mesh_center_n')
df_conc.plot(kind='hexbin', x='mesh_center_m', y='mesh_center_n', C='toxin', gridsize=25, cmap='viridis', title='spatial toxin level')   # capital C!
```

### 3.4 time series: pandas plot ordered data


```python
df_series = None
for timestep in mcdsts.get_mcds_list():
    df_cell = timestep.get_cell_df()
    df_celltype = df_cell.loc[:,['cell_type','time']].copy()
    # get count per cell type
    s_time = str(list(df_celltype.loc[:,'time'])[0])
    df_celltype.columns = ['cell_type', s_time]
    df_count = df_celltype.groupby('cell_type').count()  # pandas dataframe
    # store result
    if (df_series is None):
        df_series = df_count
    else:
        df_series = pd.merge(df_series, df_count, left_index=True, right_index=True, how='outer')
    #break

df_series
```


```python
# line plot
df_series.T.plot(kind='line', xlabel='time [min]', ylabel='count [cell]', logy=False)  # logy=True
```


```python
# area plot
df_series.T.plot(kind='area', xlabel='time [min]', ylabel='count [cell]')
```

## 4 Matplotlib Embedding of Pandas and pcdl Plots


```python
# pandas to matplotlib
fig, ax = plt.subplots(figsize=(9,6))
fig.suptitle('cell_type and toxin')
ax.axis('equal')
mcds.plot_contour('toxin', vmin=0, vmax=0.15, cmap='Blues', ax=ax)
mcds.plot_scatter(ax=ax)
plt.tight_layout()
fig.savefig(f'{s_path_2d}celltype_toxin_fusion.png', facecolor='white')
#plt.close()
```


```python
# pandas to matplotlib
fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(15,6))
fig.suptitle('cell_type and toxin')
ax = ax.ravel()
ax[0].axis('equal')
ax[1].axis('equal')
# scatter
mcds.plot_scatter(ax=ax[0])
# contour
mcds.plot_contour('toxin', vmin=0, vmax=0.15, cmap='viridis', ax=ax[1])
# finalize
plt.tight_layout()
fig.savefig(f'{s_path_2d}celltype_toxin_separate.png', facecolor='white')
#plt.close()
```

## 5 Data Analysis in Scverse - the single cell universe

### 5.1 Scanpy - basic single cell data analyis

+ pp: preprocessing
+ tl: tools
+ pl: plotting

#### 5.1.1 Scanpy:  pca, t-sne, umap, leiden


```python
# principal component analysis
sc.tl.pca(ann)  # process anndata object with the pca tool.
```


```python
#sc.pl.pca?
#sc.pl.pca(ann)  # plot pca result.
#sc.pl.pca(ann, color=list(ann.var_names)+list(ann.obs_keys()))  # gotta catch 'em all! # ncols=3 projection='3d' color_map='turbo' palette='turbo'
#sc.pl.pca(ann, color=['cell_type','current_phase'])  # plot the pca results colored by some attributes.
sc.pl.pca(ann, color='cell_type')  # plot the pca results colored by one attribute.
#sc.pl.pca_variance_ratio(ann)  # plot how much of the variation each principal component captures.
#sc.pl.pca_loadings(ann)  # principal component attribute loading (derived from sc.pl.ranking).
#sc.pl.pca_overview(ann, color=['cell_type','current_phase'])
#sc.pl.pca(ann, save='interaction_16200min_pca.png')  # plot is saved to figures directory.
```


```python
# generate and cluster neighborhood graph
sc.pp.neighbors(ann, n_neighbors=15)  # compute the neighborhood graph with the neighbors preprocess step.
sc.tl.leiden(ann, resolution=0.01)  # cluster the neighborhood graph with the leiden tool.
```


```python
sc.pl.pca(ann, color=['leiden','cell_type'])  # plot the pca results colored by leiden clusters and cell_type.
```


```python
# umap dimensional reduction embedding
sc.tl.umap(ann)  # process anndata object with the umap tool.
sc.pl.umap(ann, color=['current_phase','cell_type','leiden'])  # plot the umap result colored by some attributes.
```


```python
# t-sne dimensional reduction embedding
sc.tl.tsne(ann)  # process anndata object with the tsne tool.
sc.pl.tsne(ann, color=['current_phase','cell_type','leiden'])  # plot the tsne result colored by some attributes.
```


```python
# force-directed graph and other dimensional reduction embedding.
# force-directed is an alternative to tSNE that often preserves the topology of the data better.

# igraph layouts:
# ‘fa’ (ForceAtlas2) is default.
# ‘fr’ (Fruchterman Reingold).
## ‘grid_fr’ (Grid Fruchterman Reingold, faster than ‘fr’).
# ‘kk’ (Kamadi Kawai’, slower than ‘fr’).
# ‘lgl’ (Large Graph, very fast).
# ‘drl’ (Distributed Recursive Layout, pretty fast).
# ‘rt’ (Reingold Tilford tree layout).
# ‘rt_circular’ (Reingold Tilford circular layout).

#sc.tl.draw_graph?
sc.tl.draw_graph(ann, layout='fa')
sc.pl.draw_graph(ann, color=['current_phase', 'cell_type', 'leiden'])
```


```python
# pl.pca, pl.umap, pl.tsne, and pl.draw_graph is based on pl.embedding
sc.pl.embedding(ann, basis='X_pca', color=['current_phase', 'cell_type', 'leiden'])
```


```python
#sc.tl.embedding_density?
#sc.tl.embedding_density(ann, basis='umap')
#sc.pl.embedding_density(ann, basis='umap')

sc.tl.embedding_density(ann, basis='umap', groupby='current_phase')
sc.tl.embedding_density(ann, basis='umap', groupby='cell_type')
sc.tl.embedding_density(ann, basis='umap', groupby='leiden')
sc.pl.embedding_density(ann, basis='umap', key='umap_density_current_phase')
sc.pl.embedding_density(ann, basis='umap', key='umap_density_cell_type')
sc.pl.embedding_density(ann, basis='umap', key='umap_density_leiden')
```


```python
# leiden cluster celltype mapping
d_leiden = {
#    0: 'differentiated',
#    1: 'stem',
#    2: 'neutrophil',
#    3: 'CD8+_T_cell',
#    4: 'macrophage',
#    5: 'bacteria_a',
#    6: 'blood_vessel',
#    7: 'bacteria_b',
}
ls_label = [s_label for _, s_label in sorted(d_leiden.items())]
ann.rename_categories('leiden', ls_label)
```


```python
sc.metrics.confusion_matrix("cell_type", "leiden", ann.obs)  # pandas dataframe
import seaborn as sns
ax = sns.heatmap(sc.metrics.confusion_matrix("cell_type", "leiden", ann.obs), cmap='viridis')
ax.set_title('celltype leiden cluster confusion matrix')
```


```python
#help(sc.tl.dendrogram)
#help(sc.pl.dendrogram)
# https://flensted-mobiles.com/
# https://github.com/elmbeech/physicelldataloader/tree/master/man/img/dendrogram_mobile_rabbits.png
sc.tl.dendrogram(ann, groupby='cell_type')
sc.pl.dendrogram(ann, groupby='cell_type')
sc.tl.dendrogram(ann, groupby='leiden')
sc.pl.dendrogram(ann, groupby='leiden')
```


```python
#sc.pl.matrixplot?
#sc.pl.matrixplot(ann, var_names=ann.var_names, groupby='cell_type')
sc.pl.matrixplot(ann, var_names=ann.var_names, groupby='cell_type', cmap='RdBu_r', dendrogram=True)
#sc.pl.matrixplot(ann, var_names=ann.var_names, groupby='leiden', cmap='RdBu_r', dendrogram=True)
```


```python
#sc.pl.dotplot?
sc.pl.dotplot(ann, ann.var_names, groupby='cell_type', cmap='RdBu_r', dendrogram=True)  # cmap='magma'
#sc.pl.dotplot(ann, ann.var_names, groupby='leiden', cmap='RdBu_r', dendrogram=True)  # cmap='magma'
```


```python
#sc.pl.stacked_violin?
sc.pl.stacked_violin(ann, ann.var_names, groupby='cell_type', cmap='inferno', dendrogram=True)
#sc.pl.stacked_violin(ann, ann.var_names, groupby='leiden', cmap='inferno', dendrogram=True)
```


```python
sc.pl.heatmap(ann, var_names=ann.var_names, groupby='cell_type', cmap='RdBu_r', show_gene_labels=True,  dendrogram=True)
#sc.pl.heatmap(ann, var_names=ann.var_names, groupby='leiden', cmap='RdBu_r', show_gene_labels=True,  dendrogram=True)
```


```python
sc.pl.tracksplot(ann, var_names=ann.var_names, groupby='cell_type', dendrogram=True)
#sc.pl.tracksplot(ann, var_names=ann.var_names, groupby='leiden', dendrogram=True)
```


```python
sc.pl.clustermap(ann, cmap='RdBu_r', obs_keys='cell_type', yticklabels=False, xticklabels=True)
#sc.pl.clustermap(ann, cmap='RdBu_r', obs_keys='leiden')
```

#### 5.1.2 scanpy spatial - the precursor for squidpy


```python
# spatial
#sc.pl.spatial?  # e.g. H&E image could be put as background.
sc.pl.spatial(ann, spot_size=10)
sc.pl.spatial(ann, spot_size=10, color=['cell_type','leiden','toxin','cell_density_micron3'], ncols=2)
```


```python
# autocorrelation in 2D
# https://en.wikipedia.org/wiki/Geary's_C
# https://en.wikipedia.org/wiki/Moran%27s_I
gc = sc.metrics.gearys_c(ann)
mi = sc.metrics.morans_i(ann)
len(ann.var_names)
pd.DataFrame([gc,mi], columns=ann.var_names, index=['gc','mi'])
```

#### 5.1.3 scanpy gene expression

the basic physicell attribute data is quite different from log transformed gene expression data.
these functions might become interesting for physiboss output analysis.

genes
+ sc.pl.highest_expr_genes
+ sc.pl.filter_genes_dispersion
+ sc.pp.highly_variable_genes
+ sc.pl.highly_variable_genes

gene groups
+ sc.tl.rank_genes_groups
+ sc.tl.filter_rank_genes_groups
+ sc.pl.rank_genes_groups
+ sc.pl.rank_genes_groups_violin
+ sc.pl.rank_genes_groups_stacked_violin
+ sc.pl.rank_genes_groups_heatmap
+ sc.pl.rank_genes_groups_dotplot
+ sc.pl.rank_genes_groups_matrixplot
+ sc.pl.rank_genes_groups_tracksplot

gene score
+ sc.tl.score_genes
+ sc.tl.score_genes_cell_cycle
+ sc.pl.ranking

marker genes
+ sc.tl.marker_gene_overlap (overlap score between data-deriven marker genes and provided markers.)

simulate dynamic gene expression data
+ sc.tl.sim  (sample from a stochastic differential equation model, for example built from literature-curated boolean gene regulatory networks.)


#### 5.1.4 scanpy cell differentiation

in physicell, we do not have this problem because we know the lineage trace.
nevertheless, it might be interesting to try to apply these functions on adequate physicell output.

diffusion map - for denoising the graph
+ sc.tl.diffmap(ann, n_comps=15)
+ sc.pl.diffmap

paga - partition based graph abstraction
+ sc.tl.paga(ann, groups='leiden')
+ sc.pl.paga
+ sc.pl.paga_compare
+ sc.pl.paga_path

dpt - diffusion pseudo time
+ sc.tl.dpt
+ sc.pl.dpt_groups_pseudotime
+ sc.pl.dpt_timeseries

#### 5.1.5 scanpy other functions

in physicell, we do not have this problem because it is synthetic data.
nevertheless, it might be interesting to try to apply these functions on adequate physicell output.

+ sc.pp.combat (batch correction.)
+ sc.tl.ingest (integrates embeddings and annotations of an adata with a reference dataset adata_ref.)

#### 5.1.6 scanpy logging


```python
# version logging
sc.logging.print_header()
print()
sc.logging.print_versions()
```

#### 5.1.6 scanpy - time series data analysis


```python
# time series - not collapsed and collapsed

# not collapsed:
# a list of time step anndata objects can be processed in the same way as a single time step.
#l_ann = mcdsts.get_anndata(collapse=False)

# collapsed:
# one anndata object for the whole time can be processed in the same way as a single time step, using time as a attribute.
# this was not so obvious to me!
annts = mcdsts.get_anndata()
sc.pp.neighbors(annts, n_neighbors=15)  # compute the neighborhood graph with the neighbors preprocess step.
sc.tl.leiden(annts, resolution=0.01)  # cluster the neighborhood graph with the leiden tool.
sc.tl.umap(annts)  # process anndata object with the umap tool.
sc.pl.umap(annts, color=['time','current_phase','cell_type','leiden'], ncols=2)  # plot the umap result colored by some attributes.
```

### 5.2 Squidpy - spatial single cell data analysis

+ gr: graph
+ im: image
+ tl: tools
+ pl: plotting


```python
# load anndata object
ann = mcds.get_anndata(values=2, scale='maxabs')
print(ann)
```

#### 5.2.1 squidpy - very basic


```python
#sq.pl.spatial_scatter?
sq.pl.spatial_scatter(ann, shape=None, color="cell_type", size=42) #color=ds_color  # H&E image could be put under ann.uns as background,
print(ann)
```

#### 5.2.2 squidpy graph


```python
#sq.gr.spatial_neighbors?

sq.gr.spatial_neighbors(ann)  # generate a graph from spatial coordinates.
print(ann)
sq.pl.spatial_scatter(ann, connectivity_key="spatial_connectivities", shape=None, size=32) #color="cell_type"
```


```python
#sq.gr.nhood_enrichment?
#sq.pl.nhood_enrichment?

sq.gr.nhood_enrichment(ann, cluster_key='cell_type')  # compute neighborhood enrichment by permutation test.
print(ann)
sq.pl.nhood_enrichment(ann, cluster_key='cell_type')  # plot neighborhood enrichment.
```


```python
#sq.gr.co_occurrence?
#sq.pl.co_occurrence?

sq.gr.co_occurrence(ann, cluster_key='cell_type')  # compute co-occurrence probability of clusters.
print(ann)
sq.pl.co_occurrence(
    ann,
    cluster_key='cell_type',
    clusters=['bacteria','differentiated','blood_vessel','stem'],
)  # plot co-occurrence probability ratio for each cluster.
```


```python
#sq.gr.centrality_scores?
#help(sq.pl.centrality_scores)

sq.gr.centrality_scores(ann, cluster_key='cell_type')  # compute centrality scores per cluster or cell type.
print(ann)
sq.pl.centrality_scores(ann, cluster_key='cell_type', figsize=(12,4))  # plot centrality scores.
```


```python
#sq.gr.interaction_matrix?
#sq.pl.interaction_matrix?

sq.gr.interaction_matrix(ann, cluster_key='cell_type')  # compute interaction matrix for clusters.
print(ann)
sq.pl.interaction_matrix(ann, cluster_key='cell_type')  # plot cluster interaction matrix.
```


```python
#sq.gr.ripley?
#sq.pl.ripley?

sq.gr.ripley(ann, cluster_key='cell_type')  # calculate various Ripley's statistics for point processes.
print(ann)
sq.pl.ripley(ann, cluster_key='cell_type')  # plot Ripley's statistics for each cluster.
```


```python
#sq.gr.spatial_autocorr?

sq.gr.spatial_autocorr(ann, mode='moran')  # calculate global autocorrelation atatistic Moran’s (default).
sq.gr.spatial_autocorr(ann, mode='geary')  # calculate global autocorrelation statistic Geary's C.
print(ann)
ann.uns["moranI"].head(10)
#ann.uns["gearyC"].head(10)
```


```python
# sq.tl.var_by_distance?
# sq.pl.var_by_distance?

sq.tl.var_by_distance(ann, groups=['bacteria','blood_vessel','differentiated','stem'], cluster_key='cell_type')  # build a design matrix consisting of distance measurements to selected anchor point(s) for each observation.
print(ann)
ann.obsm["design_matrix"]
sq.pl.var_by_distance(
    ann,
    var='damage',
    anchor_key='blood_vessel',
) # plot a variable using a smooth regression line with increasing distance to an anchor point.
```

#### 5.2.5 squidpy other fuctions

functions hardly applicable on physicell output data.

**sepal** to identify spatially variable genes
works only on square or hexagonal lattice data, and because of that, not on PhysiCell output data!
+ sq.gr.sepal(ann, max_neighs[, genes, n_iter, ...]) # max_neighs=4 max_neighs=6

**ligrec** cellphonedb related ligand-receptor analysis.
+ sq.gr.ligrec(ann, cluster_key='cell_type', use_raw=False)  # perform the permutation test as described in doi:10.1038/s41596-020-0292-x.
+ sq.pl.ligrec(ann, cluster_key='cell_type')  # plot the result of a receptor-ligand permutation test.

**extract** a temporary anndata.AnnData object for plotting.
+ sq.pl.extract(adata[, obsm_key, prefix])


#### 5.2.2 squidpy imaging

for H&E and microscopy images

+ sq.im.ImageContainer
+ sq.im.process(img[, layer, library_id, method, ...]) # process an image by applying a transformation.
+ sq.im.segment(img[, layer, library_id, method, ...])  # segment an image.
+ sq.im.calculate_image_attribute(adata, img[, ...])  # calculate image attributes for all observations in adata.
+ sq.pl.spatial_segment(adata[, color, groups, ...])  # plot spatial omics data with segmentation masks on top.
+ https://napari.org/stable/  # a fast, interactive viewer for multi-dimensional images in Python

#### 5.2.4: squidpy - 3D


```python
# load data
mcds3d = pcdl.TimeStep(s_pathfile_3d)  # 24[h] = 1440[min]
ann3d = mcds3d.get_anndata()
```


```python
#ann3d
#ann3d.obs.z_layer.unique()
sc.pl.embedding(ann3d, basis="spatial", projection="3d", color=["cell_type","pressure","oxygen","z_layer"]) # ncols=2
```


```python
sq.pl.spatial_scatter(ann3d[ann3d.obs.z_layer == 0], shape=None, color=["cell_type","pressure","oxygen","z_layer"], size=128, ncols=2)
#sc.pl.spatial(ann3d[ann3d.obs.z_layer == 0], shape=None, color=["cell_type","pressure","oxygen","z_layer"], size=128)
#ann3d.uns.keys()
```

### 5.3 the other scverse libraries

https://scverse.org/

data objects
+ spatialdata  # spatial data object based on anndata.
+ mudata  # multimodal data object based on anndata.

analysis frame works
+ muon  # multi-omics analysis,
+ scvi-tools  # single cell machine learning: https://doi.org/10.1038/s41587-021-01206-w
+ scirpy  # single cell immune receptor sequence analysis: https://doi.org/10.1093/bioinformatics/btaa611
+ scverse ecosystem: https://scverse.org/packages/#ecosystem
