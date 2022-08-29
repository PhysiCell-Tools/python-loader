# python-loader

## Abstarct

Python data loader for PhysiCell digital snapshots.

Additional scripts for plotting PhysiCell output can be found in the `plotting` directory.
You can copy these into your PhysiCell `output` directory and run them there. They require the `scipy` and `matplotlib` Python modules, so we recommend installing the Anaconda Python 3.x distribution (https://www.anaconda.com/distribution/).


## HowTo Guide

How to install physicellloader?
```bash
pip3 install physicellloader
```

How to check the installed version?
```python
import physicellloader
physicellloader.__version__
```

How to load the bokehheat library?
```python
from physicellloader import pyMCDS
from physicellloader import pyMCDS_timeseries
from physicellloader import read_MultiCellDS_xml
```

## Tutoreal
+ Ref. http://www.mathcancer.org/blog/python-loader/


## Contributions
+ implementation: Randy Heiland


## release notes:
+ version 1.1.1 (2022-08-26): pip installable release derived from PhysiCell-Tools python-loader release 1.1.0 (2022-07-20).
+ version 1.1.0 (2022-05-09): release compatible with pre-v1.10.x of PhysiCell
+ version 1.0.1 (2020-01-25): time-series related bug fix
+ version 1.0.0 (2019-09-28): first public release!


