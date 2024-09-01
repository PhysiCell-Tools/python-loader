# PhysiCell Data Loader Tutorial: pcdl and Python and Napari

[Napari](https://napari.org/stable/) is both, a python library and a GUI software.
Napari is used by wetlab scientist and bioinformatician to analyse fluorescent microscopy data.
Napari can read [ome.tiff](https://www.openmicroscopy.org/ome-files/) files.
https://github.com/AllenCellModeling/napari-aicsimageio


## Install napari

And install the [aicsimageio](https://allencellmodeling.github.io/aicsimageio/) library,
which installs the "ome-types" napari plugin,
that napari can read ome.tiff images inclusive ome metadata.

```bash
pip3 install napari[all]
pip3 install aicsimageio
```


### Save pcdl data constructs from the command line into ome.tiff files

```bash
pcdl_make_ome_tiff('output/')
```


### Save pcdl data constructs from within python into ome.tiff files

```python
import pcdl

mcdsts = pcdl.TimeSeries('output/')
mcdsts.make_ome_tiff()
```


### Open ome.tiff files in napari from withon python

```python
import napari

viewer = napari.Viewer()
viewer.open('output/timeseries_ID.ome.tiff')
```


### Open ome.tiff files in napari from the comandline

```bash
napari output/timeseries_ID.ome.tiff
```

That's it. The rest is analysis within napari!
