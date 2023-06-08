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
#df.plot(ax=ax)
#ax.axis('equal')
#plt.tight_layout()
#fig.savefig(s_filename, facecolor='white')
#plt.close()


# plot stuff
def df_label_to_color(df_abc, s_label, s_cmap='viridis', b_shuffle=False):
    '''
    input:
        df_abc: dataframe to which the color column will be added.
        s_label: column name with sample labels for which a color column will be generated.
        s_cmap:  matplotlib color map label.
            https://matplotlib.org/stable/tutorials/colors/colormaps.html
        b_shuffle: should colors be given by alphabetical order,
            or should the label color mapping order be random.

    output:
        df_abc: dataframe updated with color column.

    description:
        function adds for the selected label column
        a color column to the df_abc dataframe.
    '''
    ls_label = sorted(df_abc.loc[:,s_label].unique())
    if b_shuffle:
       random.shuffle(ls_label)
    a_color = plt.get_cmap(s_cmap)(np.linspace(0, 1, len(ls_label)))
    d_color = dict(zip(ls_label, a_color))
    df_abc[f'{s_label}_color'] = [colors.to_hex(d_color[s_scene]) for s_scene in df_abc.loc[:,s_label]]

def ax_colorlegend(ax, df_abc, s_label, s_color, r_x_figure2legend_space=0.01, s_fontsize='small'):
    '''
    input:
        ax: matplotlib axis object to which a color legend will be added.
        df_abc: dataframe
        s_label: column name with sample labels.
        s_color: column name with hex colors or as web color labels.
        r_x_figure2legend_space: space between plot and legend.
            -1 is left plot border, 0 is right plot border.
        s_fontsize: font size used for the legend. known are:
            xx-small, x-small, small, medium, large, x-large, xx-large.

    output:
        ax: matplotlib axis object updated with color legend.

    description:
        function to add color legend to a figure.
    '''
    d_color = df_abc.loc[:,[s_label,s_color]].drop_duplicates().set_index(s_label).loc[:,s_color].to_dict()
    lo_patch = []
    for s_label, s_color in sorted(d_color.items()):
        o_patch = mpatches.Patch(color=s_color, label=s_label)
        lo_patch.append(o_patch)
    ax.legend(
        handles = lo_patch,
        bbox_to_anchor = (1+r_x_figure2legend_space, 0, 0, 0),
        loc = 'lower left',
        borderaxespad = 0.00,
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

