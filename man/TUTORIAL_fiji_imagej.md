# PhysiCell Data Loader Tutorial: pcdl and Fiji Imagej

[Fiji](https://fiji.sc/) is an image processing package,
a "batteries-included" distribution of [ImageJ](https://en.wikipedia.org/wiki/ImageJ),
able to read [ome.tiff](https://www.openmicroscopy.org/ome-files/) files.

Fiji Imagej is used by wetlab scientist and bioinformatician to analyse fluorescent microscopy data.

## Install Fiji Imagej ~ the 64[bit] version!

Please follow the installation instruction on the official hompage.
+ https://imagej.net/software/fiji/


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

### Linux: open ome.tiff files in Fiji Imagej from the comandline

```bash
./ImageJ-linux64 path/to/output/timeseries_ID.ome.tiff
```


### Windows:
<!-- Jenny can you write here something inteligent on how to open ome.tiff files? -->


### MacOS X:
<!-- Jenny can you write here something inteligent on how to open ome.tiff files? -->


## Running Fiji Imagej

Please follow the official documentation to learn how to run the software.
+ https://imagej.net/learn/

That's it. The rest is analysis within fiij imagej!
