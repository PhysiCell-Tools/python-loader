# pcDataLoader

## Abstract

pcDataLoader provides a pure, platform independent, python3 based, [pip](https://en.wikipedia.org/wiki/Pip_(package_manager)) installable, interface
to load output generated with the [PhysiCell](https://github.com/MathCancer/PhysiCell) agent based modeling framework
into [python3](https://en.wikipedia.org/wiki/Python_(programming_language)).

pcDataLoader is a fork from the original [PhysiCell-Tools](https://github.com/PhysiCell-Tools) [python-loader](https://github.com/PhysiCell-Tools/python-loader) implementation.

The pcDataLoader python library will maintain two main branches:

+ The **version 2 branch** will strictly be compatible with the original PhysiCell-Tools/python-loader code, although pip installable.
+ The **version 3 branch** might break with old habits, although tries to be as downwards compatible as possible.
  The aim of the v3 branch is to get a very lean and agile physicell output interface, for the ones coming from the python world to physicell.

Note: there can only be one version of pcDataLoader installed in each python environment.


## Header
+ Language: python >= 3.6
+ Library dependencies: matplotlib, numpy, pandas
+ Programmer: Patrick Wall, Elmar Bucher, Randy Heiland, Paul Macklin
+ Date of origin original PhysiCell-Tools python-loader: 2019-09-02
+ Date of origin pcDataLoader fork: 2022-08-30
+ License: [BSD-3-Clause](https://en.wikipedia.org/wiki/BSD_licenses)
+ User manual: this README.md file
+ Source code: [https://github.com/elmbeech/pcDataLoader](https://github.com/elmbeech/pcDataLoader)


## HowTo Guide - branch v2 specific

**How to install the latest pcDataLoder from the v2 branch?**
```bash
pip3 install pcDataLoder<3
```

**How to update to the latest pcDataLoder from the v2 branch?**\
This works, even when you have a v3 branch installation.
```bash
pip3 install -U pcDataLoder<3
```


## HowTo Guide - branch v3 specific

**How to install the latest pcDataLoder from the v3 branch?**
```bash
pip3 install pcDataLoder
```

**How to update to the latest pcDataLoder from the v3 branch?**\
This works, even when you have a v2 branch installation.
```bash
pip3 install -U pcDataLoder
```


## HowTo Guide - branch generic

**How to uninstall pcDataLoder from your python environment?**
```bash
pip3 uninstall pcDataLoder
```

**How to check for the current installed pcDataLoder version?**
```python
import pcDataLoder
pcDataLoder.__version__
```

**How to load the pcDataLoder library?**
```python
from pcDataLoder import pyMCDS
from pcDataLoder import pyMCDS_timeseries
from pcDataLoder import read_MultiCellDS_xml
```

**How to use the addition plotting scripts for plotting PhysiCell output?**\
This plotting scripts can be found in the [pcDataLoader/plotting](https://github.com/elmbeech/pcDataLoader/tree/master/pcDataLoader/plotting) directory.\
You can copy these python scripts into your PhysiCell `output` directory and run them there.\
`anim_svg_substrate.py` additionally requires the `scipy` library and `cells3D_fury.py` additionally requires the `fury` libraries to be installed.
Both of them can be installed with pip.

**Howto run the unit test code?**\
You can always run the unit test suit, to check if you have an integer pcDataLoader installation.
The test code can be found in the [test](https://github.com/elmbeech/pcDataLoader/tree/master/test) directory.
```bash
python3 test_snapshot.py
```
```bash
python3 test_timeseries.py
```


## Reference

This is the technical descriptions of the machinery and how to operate it.
References are maintained in each module`s [docstring](https://en.wikipedia.org/wiki/Docstring).

For each pcDataLoader module, get on the fly reference information with the help command.
```python
from pcDataLoder import pyMCDS, pyMCDS_timeseries, read_MultiCellDS_xml

help(pyMCDS)
help(pyMCDS_timeseries)
help(read_MultiCellDS_xml)
```


## Tutorial
+ http://www.mathcancer.org/blog/python-loader/


## Discussion
To be developed.


## About Documentation
+ [what nobody tells you about documentation](https://www.youtube.com/watch?v=azf6yzuJt54)
Within the pcDataLoader library, I try to stick to the policy lined out by Daniele Procida in his talk at the PyCon 2017 in Portland, Oregon.


## Contributions
+ original PhysiCell-Tools python-loader implementation: Patrick Wall, Randy Heiland, Paul Macklin
+ fork pcDataLoader implementation: Elmar Bucher


## Release Notes:
+ version 2.0.0 (2022-08-30): elmbeech/pcDataLoader pip installable release, derived from and compatible with PhysiCell-Tools/python-loader release 1.1.0 (2022-07-20).
+ version 1.1.0 (2022-05-09): Physicell-Tools/python-loader release compatible with pre-v1.10.x of PhysiCell
+ version 1.0.1 (2020-01-25): Physicell-Tools/python-loader time-series related bug fix
+ version 1.0.0 (2019-09-28): Physicell-Tools/python-loader first public release!
