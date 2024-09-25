# PhysiCell Data Loader Tutorial: pcdl and Fiji Imagej

[Fiji](https://fiji.sc/) is a free and open source image processing library,
a "batteries-included" distribution of [ImageJ](https://en.wikipedia.org/wiki/ImageJ),
able to read [ome.tiff](https://www.openmicroscopy.org/ome-files/) files.

Fiji Imagej is used by wetlab scientists and bioinformatician to analyze fluorescent microscopy data.

## Install Fiji Imagej ~ the 64[bit] version!

Please follow the installation instructions on the official homepage.
+ https://imagej.net/software/fiji/


## Generate ome.tiff files from the command line

```bash
pcdl_make_ome_tiff('output/')
```


## Generate ome.tiff files from within python

```python
import pcdl

mcdsts = pcdl.TimeSeries('output/')
mcdsts.make_ome_tiff()
```

## Linux: open ome.tiff files in Fiji Imagej from the command line

```bash
./ImageJ-linux64 path/to/output/timeseries_ID.ome.tiff
```


## Windows:
<!-- Jenny can you write here something inteligent on how to open ome.tiff files? -->


## MacOS X:
<!-- Jenny can you write here something inteligent on how to open ome.tiff files? -->


## Running Fiji Imagej

Please work through the official documentation to learn how to run the software.
+ https://imagej.net/learn/

That's it. The rest is analysis within fiij imagej!
