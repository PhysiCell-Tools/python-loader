# pcDataLoader

## Abstarct

python-loader

[Python3](https://www.python.org/) data loader for [PhysiCell](http://physicell.org/) digital snapshots.

Additional scripts for plotting PhysiCell output can be found in the `physicell/plotting` directory.
You can copy these python scripts into your PhysiCell `output` directory and run them there.
`anim_svg_substrate.py` additionally requires the `scipy` library and `cells3D_fury.py` additionally requires the `fury` libraries to be installed.
Both of them can be installed with pip.

## HowTo Guide

How to install pcDataLoder?
```bash
pip3 install pcDataLoder
```

How to check for the installed pcDataLoder version?
```python
import pcDataLoder
pcDataLoder.__version__
```

How to load the pcDataLoder library?
```python
from pcDataLoder import pyMCDS
from pcDataLoder import pyMCDS_timeseries
from pcDataLoder import read_MultiCellDS_xml
```

How to get reference information about how to use each bokehheat module?
```python
from pcDataLoder import pyMCDS, pyMCDS_timeseries, read_MultiCellDS_xml

help(pyMCDS,)
help(pyMCDS_timeseries,)
help(read_MultiCellDS_xml)
```


## Tutorial
+ http://www.mathcancer.org/blog/python-loader/


## Contributions
+ implementation: Patrick Wall, Elmar Bucher, Randy Heiland, Paul Macklin


## Release Notes:
+ version 2.0.0 (2022-09-16): Physicell-Tools/pcDataLoader pip installable release derived from PhysiCell-Tools/python-loader release 1.1.0 (2022-07-20).
+ version 1.1.0 (2022-05-09): Physicell-Tools/python-loader release compatible with pre-v1.10.x of PhysiCell
+ version 1.0.1 (2020-01-25): Physicell-Tools/python-loader time-series related bug fix
+ version 1.0.0 (2019-09-28): Physicell-Tools/python-loader first public release!
