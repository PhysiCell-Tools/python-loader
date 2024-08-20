# mcdsts.make_ome_tiff()


## input:
```
            cell_attribute: strings; default is 'ID', with will result in a segmentation mask.
                column name within cell dataframe.
                the column data type has to be numeric (bool, int, float) and can't be string.

            file: boolean; default True
                if True, an ome.tiff file is output.
                if False, a numpy array with shape tczyx is output.

            collapse: boole; default True
                should all mcds time steps from the time series be collapsed
                into one ome tiff file (numpy array),
                or an ome tiff file (numpy array) for each time step?

```

## output:
```
            a_tczyx_img: numpy array or ome.tiff file.

```

## description:
```
            function to transform chosen mcdsts output into an 1[um] spaced
            tczyx (time, channel, z-axis, y-axis, x-axis) ome tiff file or numpy array,
            one substrate or cell_type per channel.
            a ome tiff file is more or less:
            a numpy array, containing the image information
            and a xml, containing the microscopy metadata information,
            like the channel labels.
            the ome tiff file format can for example be read by the napari
            or fiji (imagej) software.

            https://napari.org/stable/
            https://fiji.sc/
        
```