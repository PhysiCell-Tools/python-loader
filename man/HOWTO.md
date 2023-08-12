# PhysiCell Data Loader How To Man Page

## How to install the latest physicelldataloader?

Full-fledged installation, with all library dependencies installed.
```bash
pip3 install pcdl[all]
```

If necessary, you can tweak your installation to make it more lightweight.
```bash
pip3 install pcdl  # The bare minimum. Installs only the pcdl core library dependencies.
pip3 install pcdl[data]  # Installs pcdl core and test data library dependencies.
pip3 install pcdl[scverse]  # Installs pcdl core and anndata library dependencies.
pip3 install pcdl[all]  # Installs pcdl core, test data, and anndata library dependencies.
```

## How to update to the latest physicelldataloader?

```bash
pip3 install -U pcdl[all]
```


## How to uninstall physicelldataloader from your python3 environment?

Note: For the pcdl library this is a two-step procedure.
First, possibly installed test data and tutorial output will be removed.
Then, the software will be uninstalled.
Keep in mind that pcdl library dependencies (like anndata, matplotlib, numpy, pandas, scipy) will never be uninstalled automatically.

```bash
python3 -c "import pcdl; pcdl.uninstall_data()"
pip3 uninstall pcdl
```


## How to load the physicelldataloader library?

```python3
import pcdl
```


## How to check for the current installed physicelldataloader version?

```python3
import pcdl
pcdl.__version__
```


## How to keep using the legacy library name pcDataLoader?

It is still possible to import pcdl with the legacy library name pcDataLoader, the way it was common before version 3.2.0 (June 2023).

```python3
import pcDataLoader as pc
```

If you need to do so, please update to the latest pcDataLoader package as follows.
The pcDataLoader library will thereafter act as a simple gateway to the latest installed pcdl library.
In future, you can just update the pcdl package to go with the latest version.

```
pip install -U pcDataLoader pcdl[all]
```


## How to run physicelldataloader like in the early days (before autumn 2022)?

1. Copy the latest [pyMCDS.py](https://raw.githubusercontent.com/elmbeech/physicelldataloader/master/pcdl/pyMCDS.py) file into the PhysiCell or PhysiCell/output folder.
2. In the same directory, fire up a python3 shell (core [python](https://docs.python.org/3/tutorial/interpreter.html#interactive-mode), [ipython](https://en.wikipedia.org/wiki/IPython), or a [jupyter](https://en.wikipedia.org/wiki/Project_Jupyter) notebook that is running an ipython kernel).
3. Load the pyMCDS class as follows:

```python3
from pyMCDS import pyMCDS
``

Now you're rolling. \
On the one hand, pyMCDS is very lightweight.
Besides the python3 core library, this code has only matplotlib, numpy, pandas, and scipy library dependencies.
On the other hand, the `pyMCDS` class lacks features present in the `pcdl.TimeStep` and `pcdl.TimeSeries` class.

