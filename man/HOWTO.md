# PhysiCell Data Loader How To Man Page


## HowTo - branch v2 specific:

### How to install the latest pcdl from the v2 branch?
```bash
pip3 install "pcdl<3"
```

### How to update to the latest physicelldataloader from the v2 branch?
Note: this works, even when you have a v3 branch installation.
```bash
pip3 install -U "pcdl<3"
```


## HowTo - branch v3 specific:

### How to install the latest physicelldataloader from the v3 branch?
```bash
pip3 install pcdl
```

### How to update to the latest physicelldataloader from the v3 branch?
Note: this works, even when you have a v2 branch installation.
```bash
pip3 install -U pcdl
```

## HowTo - branch generic:

### How to uninstall physicelldataloader from your python3 environment?
```bash
pip3 uninstall pcdl
```

### How to load the physicelldataloader library?
```python3
import pcdl
```

### How to check for the current installed physicelldataloader version?
```python3
import pcdl
pcdl.__version__
```
