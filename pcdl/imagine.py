###
# title: biotransistor.imagine.py
#
# language Python3
# license: GPLv3
# author: bue, jenny eng
# date: 2019-01-31
#
# run:
#    form pcdl import imagine
#
# description:
#    my image analysis library
#
# note:
#    load basins tiff into numpy array
#
#    import matplotlib.pyplot as plt
#    from skimage import io
#
#    a_tiff = ski.io.imread('segmentation_label_file.tiff')
#    io.imshow(a_tiff)
#    plt.imshow(a_tiff)
####

# library
import numpy as np
import pandas as pd
import sys

# function
def slide_up(a):
    """
    input:
      a: numpy array

    output:
      a: input numpy array shifted one row up.
        top row get deleted,
        bottom row of zeros is inserted.

    description:
      inspired by np.roll function, though elements that roll
      beyond the last position are not re-introduced at the first.
    """
    a = np.delete(np.insert(a, a.shape[0], 0, axis=0), 0, axis=0)
    return(a)


def slide_down(a):
    """
    input:
      a: numpy array

    output:
      a: input numpy array shifted one row down.
        top row of zeros is inserted.
        bottom row get deleted,

    description:
      inspired by np.roll function, though elements that roll
      beyond the last position are not re-introduced at the first.
    """
    a = np.delete(np.insert(a, 0, 0, axis=0), -1, axis=0)
    return(a)


def slide_left(a):
    """
    input:
      a: numpy array

    output:
      a: input numpy array shifted one column left.
        left most column gets deleted,
        right most a column of zeros is inserted.

    description:
      inspired by np.roll function, though elements that roll
      beyond the last position are not re-introduced at the first.
    """
    a = np.delete(np.insert(a, a.shape[1], 0, axis=1), 0, axis=1)
    return(a)


def slide_right(a):
    """
    input:
      a: numpy array

    output:
      a: input numpy array shifted one column right.
        left most a column of zeros is inserted.
        right most column gets deleted,

    description:
      inspired by np.roll function, though elements that roll
      beyond the last position are not re-introduced at the first.
    """
    a = np.delete(np.insert(a, 0, 0, axis=1), -1, axis=1)
    return(a)


def slide_upleft(a):
    """
    input:
      a: numpy array

    output:
      a: input numpy array shifted one row up and one column left.

    description:
      inspired by np.roll function.
    """
    a = slide_left(slide_up(a))
    return(a)


def slide_upright(a):
    """
    input:
      a: numpy array

    output:
      a: input numpy array shifted one row up and one column right.

    description:
      inspired by np.roll function.
    """
    a = slide_right(slide_up(a))
    return(a)


def slide_downleft(a):
    """
    input:
      a: numpy array

    output:
      a: input numpy array shifted one row down  and one column left.

    description:
      inspired by np.roll function.
    """
    a = slide_left(slide_down(a))
    return(a)


def slide_downright(a):
    """
    input:
      a: numpy array

    output:
      a: input numpy array shifted one row down and one column right.

    description:
      inspired by np.roll function.
    """
    a = slide_right(slide_down(a))
    return(a)


def border(ai_segment):
    """
    input:
      ai_segment: numpy array representing a cells or nuclei basin file.
        it is assumed that basin borders are represented by 0 values,
        and basins are represented with any values different from 0.
        ai_segment = skimage.io.imread("cells_basins.tif")

    output:
      ai_border: numpy array containing only the cell or nuclei basin border.
        border value will be 1, non border value will be 0.

    description:
      algorithm to extract the basin borders form basin numpy arrays.
    """
    ab_border_up = (ai_segment - slide_up(ai_segment)) != 0
    ab_border_down = (ai_segment - slide_down(ai_segment)) != 0
    ab_border_left = (ai_segment - slide_left(ai_segment)) != 0
    ab_border_right = (ai_segment - slide_right(ai_segment)) != 0
    ab_border_upleft = (ai_segment - slide_upleft(ai_segment)) != 0
    ab_border_upright = (ai_segment - slide_upright(ai_segment)) != 0
    ab_border_downleft = (ai_segment - slide_downleft(ai_segment)) != 0
    ab_border_downright = (ai_segment - slide_downright(ai_segment)) != 0
    ab_border = ab_border_up | ab_border_down | ab_border_left | ab_border_right | ab_border_upleft | ab_border_upright | ab_border_downleft | ab_border_downright
    ai_border = ab_border * 1
    return(ai_border)


def grow(ai_segment, i_step=1, b_verbose=True):
    """
    input:
      ai_segment: numpy array representing a cells basin file.
        it is assumed that basin borders are represented by 0 values,
        and basins are represented with any values different from 0.
        ai_segment = skimage.io.imread("cells_basins.tif")

      i_step: integer which specifies how many pixels the basin
        to each direction should grow.
        function can handle shrinking. enter negative steps like -1.

      b_verbose: boolean which specifies if, while processing,
          text should be outputted.

    output:
      ai_grown: numpy array with the grown basins

    description:
      algorithm to grow the basins in a given basin numpy array.
      growing happens counterclockwise, starting at noon.
    """
    ai_tree = ai_segment.copy()  # initialize output
    if (i_step > -1):
        # growing
        for i in range(i_step):
            # next grow cycle
            if b_verbose:
                print(f'grow {i+1}[px] ring ...')
            ai_treering = ai_tree.copy()
            for o_slide in [slide_up, slide_upleft, slide_left, slide_downleft, slide_down, slide_downright, slide_right, slide_upright]:
                ai_evolve = o_slide(ai_tree)
                ai_treering[(ai_evolve != ai_tree) & (ai_treering == 0)] = ai_evolve[(ai_evolve != ai_tree) & (ai_treering == 0)]
                #print(ai_treering)
            # update output
            ai_tree = ai_treering
    else:
        # shrinking
        ai_border = border(ai_segment)
        ai_membrane = grow(ai_border, i_step=abs(i_step) - 1)
        ai_tree[ai_membrane != 0] = 0

    # output
    return(ai_tree)


def grow_seed(ai_segment, i_step=1, b_verbose=True):
    """
    input:
      ai_segment: numpy array representing a cells center basin file.
        it is assumed that basin borders are represented by 0 values,
        and basins are represented with any values different from 0.
        ai_segment = skimage.io.imread("cells_basins.tif")

      i_step: integer which specifies how many pixels the seed pixel
          should grow.

      b_verbose: boolean which specifies if, while processing,
          text should be outputted.

    output:
      ai_grown: numpy array with the grown basins

    description:
      algorithm to grow the basins in a given cell center seed numpy array.
      growing happens counterclockwise, starting at noon.
    """
    # initialize output
    ai_tree = ai_segment.copy()

    # growing
    for i in range(i_step):
        # next grow cycle
        if b_verbose:
            print(f'grow {i+1}[px] ring ...')
        ai_treering = ai_tree.copy()

        # calculate pythagoras
        b_circle = ((i + 1) * 2**(1/2)) < (i_step + 1)

        # up
        ai_evolve = slide_up(ai_tree)
        ai_treering[(ai_evolve != ai_tree) & (ai_treering == 0)] = ai_evolve[(ai_evolve != ai_tree) & (ai_treering == 0)]
        # upleft
        if b_circle:
            ai_evolve = slide_upleft(ai_tree)
            ai_treering[(ai_evolve != ai_tree) & (ai_treering == 0)] = ai_evolve[(ai_evolve != ai_tree) & (ai_treering == 0)]
        # left
        ai_evolve = slide_left(ai_tree)
        ai_treering[(ai_evolve != ai_tree) & (ai_treering == 0)] = ai_evolve[(ai_evolve != ai_tree) & (ai_treering == 0)]
        # downleft
        if b_circle:
            ai_evolve = slide_downleft(ai_tree)
            ai_treering[(ai_evolve != ai_tree) & (ai_treering == 0)] = ai_evolve[(ai_evolve != ai_tree) & (ai_treering == 0)]
        # down
        ai_evolve = slide_down(ai_tree)
        ai_treering[(ai_evolve != ai_tree) & (ai_treering == 0)] = ai_evolve[(ai_evolve != ai_tree) & (ai_treering == 0)]
        # downright
        if b_circle:
            ai_evolve = slide_downright(ai_tree)
            ai_treering[(ai_evolve != ai_tree) & (ai_treering == 0)] = ai_evolve[(ai_evolve != ai_tree) & (ai_treering == 0)]
        # right
        ai_evolve = slide_right(ai_tree)
        ai_treering[(ai_evolve != ai_tree) & (ai_treering == 0)] = ai_evolve[(ai_evolve != ai_tree) & (ai_treering == 0)]
        # upright
        if b_circle:
            ai_evolve = slide_upright(ai_tree)
            ai_treering[(ai_evolve != ai_tree) & (ai_treering == 0)] = ai_evolve[(ai_evolve != ai_tree) & (ai_treering == 0)]

        # update output
        #print(ai_treering)
        ai_tree = ai_treering
    # output
    return(ai_tree)


def collision(ai_segment, i_step_size=1):
    """
    input:
      ai_segment: numpy array representing a cells basin file.
        it is assumed that basin borders are represented by 0 values,
        and basins are represented with any values different from 0.
        ai_segment = skimage.io.imread("cells_basins.tif")

    i_step_size: integer that specifies the distance from a basin
        where collisions with other basins are detected.
        increasing the step size behind > 1 will result in faster processing
        but less certain results. step size < 1 make no sense.
        default step size is 1.

    output:
        eti_collision: a set of tuples representing colliding basins.

    description:
        algorithm to detect which basin collide a given step size away.
    """
    eti_collision = set()
    for o_slide in [slide_up, slide_down, slide_left, slide_right, slide_upleft, slide_upright, slide_downleft, slide_downright]:
        ai_walk = ai_segment.copy()
        for _ in range(i_step_size):
            ai_walk = o_slide(ai_walk)
        ai_alice = ai_walk[(ai_segment != 0) & (ai_walk != 0)]
        ai_bob = ai_segment[(ai_segment != 0) & (ai_walk != 0)]
        eti_collision = eti_collision.union(set(
            zip(
                ai_alice[(ai_alice != ai_bob)],
                ai_bob[(ai_bob != ai_alice)]
            )
        ))
    # return
    return(eti_collision)


def touching_cells(ai_segment, i_border_width=0, i_step_size=1):
    """
    input:
      ai_segment: numpy array representing a cells basin file.
        it is assumed that basin borders are represented by 0 values,
        and basins are represented with any values different from 0.
        ai_segment = skimage.io.imread("cells_basins.tif")

      i_border_width: maximal acceptable border with in pixels.
        this is half of the range how far two the adjacent cell maximal
        can be apart and still are regarded as touching each other.

      i_step_size: step size by which the border width is sampled for
        touching cells.
        increase the step size behind > 1 will result in faster processing
        but less certain results. step size < 1 make no sense.
        default step size is 1.

    output:
      dei_touch: a dictionary that for each basin states
        which other basins are touching.

    description:
      algorithm to extract the touching basins from a cell basin numpy array.
      algorithm inspired by C=64 computer games with sprit collision.
    """

    # detect neighbors
    eti_collision = set()
    ai_evolve = ai_segment.copy()
    for _ in range(-1, i_border_width, i_step_size):
        # detect cell border collision
        eti_collision = eti_collision.union(
            collision(ai_segment=ai_evolve, i_step_size=i_step_size)
        )
        # grow basin
        ai_evolve = grow(ai_segment=ai_evolve, i_step=i_step_size)

    # transform set of tuple of alice and bob collision to dictionary of sets
    dei_touch = {}
    ei_alice = set(np.ndarray.flatten(ai_segment))
    ei_alice.remove(0)
    for i_alice in ei_alice:
        dei_touch.update({i_alice : set()})
    for i_alice, i_bob in eti_collision:
        ei_bob = dei_touch[i_alice]
        ei_bob.add(i_bob)
        dei_touch.update({i_alice : ei_bob})

    # output
    return(dei_touch)


def detouch_to_df(deo_touch, ls_column=["cell_center","cell_touch"]):
    """
    input:
        deo_touch: touching_cells generated dictionary
        ls_column: future dictionary_key dictionary_value column name

    output:
        df_touch: dataframe which contains the same information
          as the input deo_touch dictionary.

    description:
        transforms dei_touch dictionary into a two column dataframe.
    """
    lo_key_total= []
    lo_value_total = []
    for o_key, eo_value in deo_touch.items():
        try:
            lo_value = sorted(eo_value, key=int)
        except ValueError:
            lo_value = sorted(eo_value)
        # extract form dictionary
        if (len(lo_value) == 0):
            lo_key_total.append(o_key)
            lo_value_total.append(0)
        else:
            lo_key_total.extend([o_key] * len(lo_value))
            lo_value_total.extend(lo_value)
    # generate datafarme
    df_touch = pd.DataFrame([lo_key_total,lo_value_total], index=ls_column).T
    return(df_touch)


# bue: 202021016: maybe refracture this into segment_px and membrane_px function
# the segement px function can then be used too on simple cell and  nucles data
def membrane_px(ai_segment, dai_value, i_step=1, b_approximation=False):
    '''
    input:
      ai_segment: numpy array representing a cells basin tiff file.
        it is assumed that basin borders are represented by 0 values,
        and basins are represented with any values different from 0.
        ai_segement = skimage.io.imread("cells_basins.tif")

      dai_value: dictionary of numpy array representing a
        protein expression value tiff file.
        the dictionary key should be the protein label.

      i_step: number of pixel the cell border, which is 2 pixel,
        in both direction is grown, to cover the membrane segment.
        default is 1.

      b_approximation: if set True, calculation works with
        non-overlapping membranes and will as such be much faster.

    output:
      df_membrane_px: dictionary of pandas datafarame whit cell_id,
        image absolute xy coordinate, extended cell relative coordinate,
        and membrane pixel related expression values for all proteins.

    description:
        function extracts protein membrane expression values,
        for each protein submitted to the function.
    '''
    # get a fix sensor order
    ls_sensor = sorted(dai_value.keys())

    # get cell border
    ab_border = border(ai_segment).astype(bool)
    ai_segment_border = ai_segment.copy()
    ai_segment_border[~ab_border] = 0

    # aproximative algorithm
    if (b_approximation) or (i_step <= 0):

        # membrane segment
        ai_membrane = grow(ai_segment_border, i_step=i_step)
        ab_membrane = ai_membrane.astype(bool)
        ai_cell = ai_membrane[ab_membrane]

        # coordinates
        tai_coor_absolute = np.where(ab_membrane)
        ai_ycoor_absolute = tai_coor_absolute[0]
        ai_xcoor_absolute = tai_coor_absolute[1]

        ai_membrane_coor = np.stack([
            ai_cell,
            ai_ycoor_absolute,
            ai_xcoor_absolute,
        ], axis=1)

        # values
        lai_membrane_value = []
        for s_sensor in ls_sensor:
            ai_value = dai_value[s_sensor]
            lai_membrane_value.append(
                ai_value[ab_membrane]
            )
        ai_membrane_value = np.stack(
            lai_membrane_value,
            axis=1,
        )

        # concatenate coor and value
        ai_membrane_px = np.concatenate([
            ai_membrane_coor,
            ai_membrane_value,
        ], axis=1)

        # pack output
        df_membrane_px = pd.DataFrame(
            ai_membrane_px,
            columns = ['cell', 'y_absolute', 'x_absolute'] + ls_sensor,
        )#.astype({})
        df_membrane_px.index.name = f'approximation_membrane_{2*i_step + 1}px'
        print(df_membrane_px.info())

    # exact algorithm
    else:
        # empty result object
        ai_membrane_px = None

        # for each cell
        ei_cell = set(ai_segment.flatten())
        ei_cell.discard(0) # kick cell0 which is background
        i_total = len(ei_cell)
        #for i, i_cell in enumerate(sorted(ei_cell)[0:16]):
        for i, i_cell in enumerate(sorted(ei_cell)):
            #print(f'processing cell{i_cell}: {i} / {i_total}')

            # membrane segment
            tai_coor = np.where(ai_segment_border == i_cell)
            ai_ycoor = tai_coor[0]
            ai_xcoor = tai_coor[1]
            i_ymin = np.array([ai_ycoor.min() - i_step]).clip(min=0)[0]
            i_xmin = np.array([ai_xcoor.min() - i_step]).clip(min=0)[0]
            i_ymax = np.array([ai_ycoor.max() + i_step + 1]).clip(max=ai_segment_border.shape[0])[0]
            i_xmax = np.array([ai_xcoor.max() + i_step + 1]).clip(max=ai_segment_border.shape[1])[0]
            ai_border = (ai_segment_border == i_cell)[i_ymin:i_ymax,i_xmin:i_xmax].astype(int)
            ab_membrane = grow(ai_border, i_step=i_step).astype(bool)

            # coordiantes
            tai_coor_relative = np.where(ab_membrane)
            ai_ycoor_relative = tai_coor_relative[0]
            ai_xcoor_relative = tai_coor_relative[1]
            ai_ycoor_absolute = ai_ycoor_relative + i_ymin
            ai_xcoor_absolute = ai_xcoor_relative + i_xmin
            ai_cell = np.array([i_cell] * len(ai_ycoor_relative))
            ai_membrane_coor = np.stack([
                ai_cell,
                ai_ycoor_relative,
                ai_xcoor_relative,
                ai_ycoor_absolute,
                ai_xcoor_absolute,
            ], axis=1)

            # values
            lai_membrane_value = []
            for s_sensor in ls_sensor:
                ai_value = dai_value[s_sensor]
                lai_membrane_value.append(
                    ai_value[i_ymin:i_ymax,i_xmin:i_xmax][ab_membrane]
                )
            ai_membrane_value = np.stack(
                lai_membrane_value,
                axis=1,
            )

            # concatenate coor and value
            ai_membrane_coorvalue = np.concatenate([
                ai_membrane_coor,
                ai_membrane_value,
            ], axis=1)

            # update ouput
            if (ai_membrane_px is None):
                ai_membrane_px = ai_membrane_coorvalue
            else:
                ai_membrane_px = np.concatenate([
                    ai_membrane_px,
                    ai_membrane_coorvalue,
                ], axis=0)

            # sanity check
            print(f'cell {i_cell} / {i_total}: min {ai_ycoor_relative.min()},{ai_xcoor_relative.min()} x max {ai_ycoor_relative.max()},{ai_xcoor_relative.max()}')

        # pack output
        df_membrane_px = pd.DataFrame(
            ai_membrane_px,
            columns = ['cell', 'y_relative', 'x_relative', 'y_absolute', 'x_absolute'] + ls_sensor,
        )#.astype({})
        df_membrane_px.index.name = f'membrane_{2*i_step + 2}px'

    # output
    return(df_membrane_px)


# bue 20201016: i should add all statistical moments
# bue 20201024: this should become dfpx_stats and be applicable to cell, nucleus, cytoplasm, and membrane
def dfmembranepx_stats(df_membrane_px, di_threshold_raw=None):
    '''
    input:
      df_membrane_px: membrane_px output which is a
        pandas datafarame whit cell_id,
        image absolute and cell relative xy coordinate and
        membrane pixel related expression values for all proteins.

      di_threshold_raw: dictionary of raw intensity threshold values
        for each protein. if set to none, cell membrane positive fraction
        will not be calculated and be missing from the output.
        default is None.

    output:
      ddf_out: dictionary of pandas datafarames. one dataframe per protein.
        each dataframe stores mean, min, max and whole range of quantile values
        measured for each cell membrane.

    description:
        function calculates statistical key numbers for membrane protein
        expression values, for each protein submitted to the function.
    '''
    # handle input
    es_sensor = set(df_membrane_px.columns)
    es_sensor = es_sensor.difference(
        {'cell', 'y_relative', 'x_relative', 'y_absolute', 'x_absolute'}
    )

    # empty result object
    ddlr_membrane = {}
    for s_sensor in es_sensor:
        ddlr_membrane.update({s_sensor: {}})

    # for each cell
    ei_cell = set(df_membrane_px.cell.unique())
    i_total = len(ei_cell) - 1

    #for i, i_cell in enumerate(sorted(ei_cell)[0:16]):
    for i, i_cell in enumerate(ei_cell):
        print(f'processing cell{i_cell}: {i} / {i_total}')

        # for each sensor
        df_membrane_cell_px = df_membrane_px.loc[df_membrane_px.cell == i_cell,:]
        for s_sensor in es_sensor:
            ai_membrane = df_membrane_cell_px.loc[:,s_sensor].values

            # quantile calculation
            lr_membrane = [
                np.mean(ai_membrane),
                np.std(ai_membrane, ddof=0),
                np.min(ai_membrane),
                np.quantile(ai_membrane, 0.01),
                np.quantile(ai_membrane, 0.05),
                np.quantile(ai_membrane, 0.1),
                np.quantile(ai_membrane, 0.2),
                np.quantile(ai_membrane, 0.3),
                np.quantile(ai_membrane, 0.4),
                np.quantile(ai_membrane, 0.5),
                np.quantile(ai_membrane, 0.6),
                np.quantile(ai_membrane, 0.7),
                np.quantile(ai_membrane, 0.8),
                np.quantile(ai_membrane, 0.9),
                np.quantile(ai_membrane, 0.95),
                np.quantile(ai_membrane, 0.99),
                np.max(ai_membrane),
                ai_membrane.shape[0],
            ]

            # membrane fraction positive calculation
            if not (di_threshold_raw is None):
                ab_membrane = ai_membrane > di_threshold_raw[s_sensor]
                i_sumpos = ab_membrane.sum()
                r_fractpos = i_sumpos / ab_membrane.shape[0]
                lr_membrane.extend([i_sumpos, r_fractpos])

            # update result object
            ddlr_membrane[s_sensor].update({i_cell: lr_membrane})

    # pack output
    ls_index = [
        'mean','std',
        'min','q0_01','q0_05',
        'q0_1','q0_2','q0_3','q0_4','q0_5',
        'q0_6','q0_7','q0_8','q0_9','q0_95',
        'q0_99','max',
        'membrane_size_px'
    ]
    if not (di_threshold_raw is None):
        ls_index.extend(['px_pos','fract_pos'])
    ddf_out = {}
    for s_sensor, dlr_membrane in ddlr_membrane.items():
        df_membrane = pd.DataFrame(
            dlr_membrane,
            index=ls_index,
        ).T
        ddf_out.update({s_sensor: df_membrane})
        df_membrane.index.name = f'cell'

    # output
    return(ddf_out)


def imgfuse(laaai_in):
    """
    input:
        laaai_in: list of 3 channel (RGB) images

    output:
       aaai_out: fused 3 channel image

    description:
       code to fuse many RGB images into one.
    """
    # check shape
    ti_shape = None
    for aaai_in in laaai_in:
        if (ti_shape is None):
            ti_shape = aaai_in.shape
        else:
           if (aaai_in.shape != ti_shape):
               sys.exit(f"Error: input images have not the same shape. {aaai_in.shape} != {aaai_in}.")

    # fuse images
    llli_channel = []
    for i_channel in range(ti_shape[0]):
        lli_matrix = []
        for i_y in range(ti_shape[1]):
            li_row = []
            for i_x in range(ti_shape[2]):
                #print(f"{i_channel} {i_y} {i_x}")
                li_px = []
                for aaai_in in laaai_in:
                    i_in = aaai_in[i_channel,i_y,i_x]
                    if (i_in != 0):
                        li_px.append(i_in)
                if (len(li_px) != 0):
                    i_out = np.mean(li_px)
                else:
                    i_out = 0
                li_row.append(int(i_out))
            lli_matrix.append(li_row)
        llli_channel.append(lli_matrix)

    # output
    aaai_out = np.array(llli_channel)
    return(aaai_out)


# main code
#if __name__ == "__main__":

