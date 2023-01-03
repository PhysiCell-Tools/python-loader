# pcDataLoader HowTo Man Page


## HowTo - branch v2 specific:

### How to install the latest pcDataLoader from the v2 branch?
```bash
pip3 install "pcDataLoader<3"
```

### How to update to the latest pcDataLoader from the v2 branch?
This works, even when you have a v3 branch installation.
```bash
pip3 install -U "pcDataLoader<3"
```


## HowTo - branch v3 specific:

### How to install the latest pcDataLoader from the v3 branch?
```bash
pip3 install pcDataLoader
```

### How to update to the latest pcDataLoader from the v3 branch?
This works, even when you have a v2 branch installation.
```bash
pip3 install -U pcDataLoader
```

### Howto run the unit test code?
The pcDataLoader version3 sorce code has full unit test coverage.\
Unitest were written with in [pytest](https://docs.pytest.org/).\
All unit test code is in the [test](https://github.com/elmbeech/pcDataLoader/tree/master/test) directory.
```python3
import os
import pathlib
import pcDataLoader as pc
s_path = str(pathlib.Path(pc.__file__).parent.parent.resolve())  # get path to the installed libray
os.system(f'pytest {s_path}')  # run tests
```


## HowTo - branch generic:

### How to uninstall pcDataLoader from your python3 environment?
```bash
pip3 uninstall pcDataLoader
```

### How to load the pcDataLoader library?
```python3
import pcDataLoader as pc
```

### How to check for the current installed pcDataLoader version?
```python3
import pcDataLoader as pc
pc.__version__
```

