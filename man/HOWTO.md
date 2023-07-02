# PhysiCell Data Loader How To Man Page

## How to install the latest physicelldataloader?
```bash
pip3 install pcdl
```

## How to update to the latest physicelldataloader?
```bash
pip3 install -U pcdl
```

## How to uninstall physicelldataloader from your python3 environment?
Note: For the pcdl library this is a two-step procedure.
First, possibly installed test data and tutorial output will be removed.
Then, the software will be uninstalled.
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
