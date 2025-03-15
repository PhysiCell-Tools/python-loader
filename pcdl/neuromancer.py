#########
# title: ometiff2neuro.py
#
# language: python3
# date: 2022-02-16
# license: MIT
# author: jason lu, viviana kwong, elmar bucher
#
# installation:
#     pip install neuroglancer   # the basics
#     pip install ipython   # for coding
#     pip install matplotlib   # for expression value color maps
#     pip install numpy   # for basic code manipulation
#     pip install scikit-image   # signal thresh and image manipulation.
#     pip install bioio bioio-ome-tiff   # for loading the image and extracting ometiff metadata
#
# run:
#     ipython3 -i ometiff2neuro.py <path/filename.ome.tiff>
#
# description:
#   script to render multi channel multi slice ome.tiff files into the
#   neuroglancer software.
#
#   the script here makes use of the neuroglancer python library and is based
#   on the neuroglancer/python/examples/example.py code.
#   with this script, it is possible to load three-dimensional single
#   time step ome-tiff files straight into the neuroglancer software.
#   for each channel, mesh generation, rendering, and expression intensity
#   coloring is done on the fly.
#   channels can be viewed together or alone by toggling them on or off in
#   the neuroglancer user interface.
#
# references:
#   + https://www.openmicroscopy.org/ome-files/
#   + https://github.com/google/neuroglancer
#   + https://github.com/google/neuroglancer/tree/master/python
#########


# library
import argparse
from bioio import BioImage
import matplotlib as mpl
import neuroglancer
import neuroglancer.cli
import numpy as np
from skimage import exposure, filters, util
import sys


# functions
def ometiff2neuro(
        o_state,
        s_pathfile_tiff,
        s_intensity_cm = 'gray',  # None
        b_intensity_norm = False,  # True # False
        o_thresh = None,  # filters.threshold_li, # None
        e_render = None,  # {0, 'CD4', 'CD8A', 'GZMB', 28}, # None
    ):
    '''
    input:
        o_state: neuroglancer viewer state object.

        s_pathfile_tiff: file name and path to input tiff file.

        s_intensity_cm: matlab color map object, used to display expression intensity values.
            default is cm.gray.
            if None, no intensity layers will be generated.

        b_intensity_norm: boolean.
            default is True.
            if True, expression intensity values will be cut by the 2.5 and 97.5 percentiles, to remove outliers,
                and intensity will be stretched over the whole range.

        o_thresh: skimage.filter thresh method object to threshold non-threshed data.
            default is filters.threshold_li.
            this is dataset dependent information.
            set None if the tiff contains no raw but already threshed or segmentation mask data!


        e_render: set of channel marker labels and/or integers to specify which markers should be rendered into neurogalncer.
            default is None.
            if None, all channels will be rendered into neuroglancer (if enough RAM is available).

    output:
        local url where the rendered data can be viewed.

    description:
        function to render multi channel multi slice (ome)tiff files into the neuroglancer software.
    '''
    # off we go #
    print(f'\nprocessing: {s_pathfile_tiff}')

    # load tiff image as numpy array
    o_img =  BioImage(s_pathfile_tiff)

    # extract dimension order
    di_coor = {}
    [di_coor.update({coor: i}) for i, coor in enumerate(list(o_img.dims.order.lower()))]
    i_c = di_coor['c']
    i_x = di_coor['x']
    i_y = di_coor['y']
    i_z = di_coor['z']
    if di_coor['x'] > i_c:
        i_x -= 1
    if di_coor['y'] > i_c:
        i_y -= 1
    if di_coor['z'] > i_c:
        i_z -= 1

    # extract channel count
    i_channel = o_img.shape[di_coor['c']]

    # handle set with labels from channels that will be rendered
    if e_render is None:
        e_render = set(range(i_channel))

    # extract data for time step #
    if set(di_coor.keys()) == {'t','c','z','y','x'}:
        i_t = di_coor['t']
        if di_coor['x'] > i_t:
            i_x -= 1
        if di_coor['y'] > i_t:
            i_y -= 1
        if di_coor['z'] > i_t:
            i_z -= 1

        i_timestep = o_img.shape[di_coor['t']] - 1
        if i_timestep != 0:
            sys.exit(f'Error @ ometiff2neuro : {i_timestep} time steps detected. cannot handle ome tiff with more or less than 1 time step.')

        if i_t == 0:
            a_img =  o_img.data[i_timestep,:,:,:,:]
        elif i_t == 1:
            a_img =  o_img.data[:,i_timestep,:,:,:]
        elif i_t == 2:
            a_img =  o_img.data[:,:,i_timestep,:,:]
        elif i_t == 3:
            a_img =  o_img.data[:,:,:,i_timestep,:]
        elif i_t == 4:
            a_img =  o_img.data[:,:,:,:,i_timestep]
        else:
            sys.exit(f'Error @ ometiff2neuro : the extract time setp source code is broken. please, fix the script.')

    elif set(di_coor.keys()) == {'c','z','y','x'}:
        a_img = o_img.data

    else:
        sys.exit(f'Error : cannot handle ome tiff file with dimension {o_img.dims.order.lower()}')

    # generate neuroglancer layers #
    for i_n in range(i_channel):
        s_label = o_img.channel_names[i_n]
        print(f'check channel {i_n}/{i_channel}: {s_label}')
        if (i_n in e_render) or (s_label in e_render):
            print(f'rendering channel {i_n}/{i_channel}: {s_label}')

            # extract data and xyz coordinate columns
            if i_c == 0:
                a_channel = a_img[i_n,:,:,:]
            elif i_c == 1:
                a_channel = a_img[:,i_n,:,:]
            elif i_c == 2:
                a_channel = a_img[:,:,i_n,:]
            elif i_c == 3:
                a_channel = a_img[:,:,:,i_n]
            else:
                sys.exit(f'Error @ ometiff2neuro : the extract channel source code is broken. please, fix the script.')

            # 3D rendering #
            # thresh data
            a_thresh = a_channel.copy()
            if not (o_thresh is None):
                r_thresh = o_thresh(a_thresh)
                a_thresh[a_thresh < r_thresh] = 0
            ab_thresh = a_thresh > 0

            # shape
            a_shape = np.zeros(ab_thresh.shape, dtype=np.uint32)
            a_shape[ab_thresh] = i_n + 1

            # generate neuroglancer object
            ls_name = [None, None, None]
            ls_name[i_x] = 'x'
            ls_name[i_y] = 'y'
            ls_name[i_z] = 'z'
            li_scale = [None, None, None]
            li_scale[i_x] = o_img.physical_pixel_sizes.X
            li_scale[i_y] = o_img.physical_pixel_sizes.Y
            li_scale[i_z] = o_img.physical_pixel_sizes.Z
            o_state.layers.append(
                name = s_label,
                layer = neuroglancer.LocalVolume(
                    data = a_shape,
                    dimensions = neuroglancer.CoordinateSpace(
                        # rgb, x, y, z
                        names = ls_name,
                        scales = li_scale,
                        units = ['nm', 'nm', 'nm'],
                    ),
                ),
                visible = False,  # special thanks to Jeremy Maitin-Shepard.
            )

            # expression intensity rendering #
            if not (s_intensity_cm is None):

                # normalize expression values by clip by two sigma and scale over the whole uint range
                if (b_intensity_norm):
                    i_min_clip = int(np.percentile(a_channel, 2.5))
                    i_max_clip = int(np.percentile(a_channel, 97.5))
                    a_clipped = np.clip(a_channel, a_min=i_min_clip, a_max=i_max_clip)
                    a_channel = exposure.rescale_intensity(a_clipped, in_range='image')  # 16 or 8[bit] normalized
                else:
                    a_channel = exposure.rescale_intensity(a_channel, in_range='image')  # 16 or 8[bit] normalized

                # translate intensity by color map
                a_8bit = util.img_as_ubyte(a_channel)
                a_intensity =  mpl.colormaps[s_intensity_cm](a_8bit, alpha=None, bytes=True)[:,:,:,0:3]

                # generate neuroglancer object
                ls_name = [None, None, None,'c^']
                ls_name[i_x] = 'x'
                ls_name[i_y] = 'y'
                ls_name[i_z] = 'z'
                li_scale = [None, None, None, 3]
                li_scale[i_x] = o_img.physical_pixel_sizes.X
                li_scale[i_y] = o_img.physical_pixel_sizes.Y
                li_scale[i_z] = o_img.physical_pixel_sizes.Z
                o_state.layers.append(
                    name = s_label + '_intensity',
                    layer = neuroglancer.LocalVolume(
                        data = a_intensity,
                        dimensions = neuroglancer.CoordinateSpace(
                            # rgb, x, y, z
                            names = ls_name,
                            scales = li_scale,
                            units = ['nm', 'nm', 'nm', ''],
                        ),
                    ),
                    visible = False,  # special thanks to Jeremy Maitin-Shepard.
                    shader=
"""
void main() {
  emitRGB(
    vec3(toNormalized(getDataValue(0)),
    toNormalized(getDataValue(1)),
    toNormalized(getDataValue(2)))
  );
}
""",
            )


# run the code from the command line
if __name__ == '__main__':
    import neuroglancer
    o_parser = argparse.ArgumentParser(description='Script to render an ome.tiff file into the neuroglancer software.')
    # request path to ometiff and file name as command line argument
    o_parser.add_argument('ometiff', type=str, nargs=1, help='ome.tiff path/filename')
    # start neuroglancer
    neuroglancer.cli.add_server_arguments(o_parser)
    neuroglancer.cli.handle_server_arguments(o_parser.parse_args())
    viewer = neuroglancer.Viewer()
    with viewer.txn() as state:
        # render ometiff
        ometiff2neuro(
            o_state = state,
            s_pathfile_tiff = o_parser.parse_args().ometiff[0],
        )
    # print neuroglancer viewer url
    print(viewer)
