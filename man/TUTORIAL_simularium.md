# PhysiCell Data Loader Tutorial: pcdl and Simularium

The [simularium](https://simularium.allencell.org/) viewer was developed to share and interrogate interactive 3D visualizations of biological simulation trajectories and related plots directly in a web browser.

Simularium files can easily be generated from mcds time series.


### &#x2728; command line

#### command line time series

Generate a simularium trajectory file.

```bash
pcdl_make_simularium output_2d/
```

#### command line man page

```bash
pcdl_make_simularium -h
```


### &#x2728; python

#### python time series

Load a time series and generate a simularium trajectory file.

```python
import pcdl

mcdsts = pcdl.TimeSeries('output_2d/')
mcdsts.make_simularium()
```

#### python docstrings

```python
import pcdl

mcdsts = pcdl.TimeSeries('output_2d/')
help(mcdsts.make_simularium())
```

### &#x2728; simularium viewer

1. Open the simularium viewer web page: https://simularium.allencell.org/
2. On the top right corner, click the "Load models" dropdown menu and choose "Simularium file".
3. In the popup window choose "From your device" "Select File".
4. Select the simularium file and click "Open".
5. Click "Load".
6. Explore the loaded file.


### &#x2728; Further readings'

Please work through the official documentation to learn more about how simularium software.

+ https://simularium.allencell.org/
+ https://github.com/simularium
+ https://simularium.github.io/simulariumio/#
+ https://doi.org/10.1038/s41592-022-01442-1

That's it. The rest is analysis within simularium!
