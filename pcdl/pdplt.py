####
# title: biotransistor.pdplt.py
#
# language: python3
# date: 2019-06-29
# license: GPL>=v3
# author: Elmar Bucher
#
# description:
#    library with the missing pandas plot features.
#    https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.plot.html
####


# library
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import matplotlib.patches as mpatches
import numpy as np
import random


# pandas to matplotlib
#fig, ax = plt.subplots()
#ax = ax.ravel()
#ax.axis('equal')
#df.plot(ax=ax)
#plt.tight_layout()
#fig.savefig(s_filename, facecolor='white')
#plt.close()


# plot stuff
def df_label_to_color(df_abc, s_focus, es_label=None, s_nolabel='gray', s_cmap='viridis', b_shuffle=False):
    '''
    input:
        df_abc: dataframe to which the color column will be added.
        s_focus: column name with sample labels for which a color column will be generated.
        es_label: set of labels to color. if None, es_label will be extracted for the s_focus column.
        s_nolabel: color for labels not defined in es_label.
        s_cmap:  matplotlib color map label.
            https://matplotlib.org/stable/tutorials/colors/colormaps.html
        b_shuffle: should colors be given by alphabetical order,
            or should the label color mapping order be random.

    output:
        df_abc: dataframe updated with color column.
        ds_color: lable to hex color string mapping dictionary

    description:
        function adds for the selected label column
        a color column to the df_abc dataframe.
    '''
    if (es_label is None):
        es_label = set(df_abc.loc[:,s_focus])
    if b_shuffle:
       ls_label = list(es_label)
       random.shuffle(ls_label)
    else:
       ls_label = sorted(es_label)
    a_color = plt.get_cmap(s_cmap)(np.linspace(0, 1, len(ls_label)))
    do_color = dict(zip(ls_label, a_color))
    df_abc[f'{s_focus}_color'] = s_nolabel
    ds_color = {}
    for s_category, o_color in do_color.items():
        s_color = colors.to_hex(o_color)
        ds_color.update({s_category : s_color})
        df_abc.loc[(df_abc.loc[:,s_focus] == s_category), f'{s_focus}_color'] = s_color
    # output
    return(ds_color)

def ax_colorlegend(ax, ds_color, s_loc='lower left', s_fontsize='small'):
    '''
    input:
        ax: matplotlib axis object to which a color legend will be added.
        ds_color: lables to color strings mapping dictionary
        s_loc: the location of the legend.
            possible strings are: best,
            upper right, upper center, upper left, center left,
            lower left, lower center, lower right, center right,
            center.
        s_fontsize: font size used for the legend. known are:
            xx-small, x-small, small, medium, large, x-large, xx-large.

    output:
        ax: matplotlib axis object updated with color legend.

    description:
        function to add color legend to a figure.
    '''
    lo_patch = []
    for s_label, s_color in sorted(ds_color.items()):
        o_patch = mpatches.Patch(color=s_color, label=s_label)
        lo_patch.append(o_patch)
    ax.legend(
        handles = lo_patch,
        loc = s_loc,
        fontsize = s_fontsize
    )

def ax_colorbar(ax, r_vmin, r_vmax, s_cmap='viridis', s_text=None, o_fontsize='medium', b_axis_erase=False):
    '''
    input:
        ax: matplotlib axis object to which a colorbar will be added.
        r_vmin: colorbar min value.
        r_vmax: colorbar max value.
        s_cmap: matplotlib color map label.
            https://matplotlib.org/stable/tutorials/colors/colormaps.html
        s_text: to label the colorbar axis.
        o_fontsize: font size used for the legend. known are:
            xx-small, x-small, small, medium, large, x-large, xx-large.
        b_axis_erase: should the axis ruler be erased?

    output:
        ax: matplotlib axis object updated with colorbar.

    description"
        function to add colorbar to a figure.
    '''
    if b_axis_erase:
        ax.axis('off')
    if not (s_text is None):
        ax.text(0.5,0.5, s_text, fontsize=o_fontsize)
    plt.colorbar(
        cm.ScalarMappable(
            norm=colors.Normalize(vmin=r_vmin, vmax=r_vmax, clip=False),
            cmap=s_cmap,
        ),
        ax=ax,
    )

