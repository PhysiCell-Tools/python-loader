# PhysiCell Data Loader Tutorial: pcdl and Python and the Ome.tiff, Tiff, Png, and Jpeg File Format

In pcdl output from [TimeSteps](https://github.com/elmbeech/physicelldataloader/blob/master/man/TUTORIAL_python3_timestep.md) and [TimeSeries](https://github.com/elmbeech/physicelldataloader/blob/master/man/TUTORIAL_python3_timeseries.md) can be stotred as [tiff](https://www.loc.gov/preservation/digital/formats/fdd/fdd000022.shtml), [png](http://libpng.org/pub/png/), and [jpeg](https://jpeg.org/jpeg/) files.

Tiff, png, and jpeg are raster graphic file formats.

Additionally tiff, png, and jpeg files can as well be loaded back in to python as [numpy](https://numpy.org/) array, for example with the [sci-kit image](https://scikit-image.org/) library (image data only).


### Save pcdl data constructs from the command line into tiff files

```bash
pcdl_plot_contour output/output00000012.xml oxygen
```
```bash
pcdl_plot_scatter output/output00000012.xml
```


### Save pcdl data constructs from within python into tiff files

```python
import pcdl

mcdsts = pcdl.TimeSeries('output/')
mcdsts.get_mcds_list()[12].plot_contour(focus='oxygen', ext='tiff')
mcdsts.get_mcds_list()[12].plot_scatter(ext='tiff')
```


### Load tiff files as numpy array into python

```python
from skimage import io

a_conc = io.imread('output/conc_oxygen_z0.0/output00000012_oxygen.tiff')
a_conc.shape  # (480, 640, 4)
```
```python
from skimage import io

a_cell = io.imread('output/cell_cell_type_z0.0/output00000012_cell_type.tiff')
a_cell.shape  # (480, 640, 4)
```


That's it. The rest is analysis!
