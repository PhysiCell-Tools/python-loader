# Pcdataloader How To Man Page


## HowTo - branch v2 specific:

### How to install the latest pcdataloader from the v2 branch?
```bash
pip3 install "pcdataloader<3"
```

### How to update to the latest pcdataloader from the v2 branch?
Note: this works, even when you have a v3 branch installation.
```bash
pip3 install -U "pcdataloader<3"
```


## HowTo - branch v3 specific:

### How to install the latest pcdataloader from the v3 branch?
```bash
pip3 install pcdataloader
```

### How to update to the latest pcdataloader from the v3 branch?
Note: this works, even when you have a v2 branch installation.
```bash
pip3 install -U pcdataloader
```

## HowTo - branch generic:

### How to uninstall pcdataloader from your python3 environment?
```bash
pip3 uninstall pcdataloader
```

### How to load the pcdataloader library?
```python3
import pcdataloader as pc
```

### How to check for the current installed pcdataloader version?
```python3
import pcdataloader as pc
pc.__version__
```

