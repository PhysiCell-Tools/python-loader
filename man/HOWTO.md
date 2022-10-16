# pcDataLoader HowTo Man Page


## HowTo - branch v2 specific:

**How to install the latest pcDataLoader from the v2 branch?**
```bash
pip3 install "pcDataLoader<3"
```

**How to update to the latest pcDataLoader from the v2 branch?**\
This works, even when you have a v3 branch installation.
```bash
pip3 install -U "pcDataLoader<3"
```


## HowTo - branch v3 specific:

**How to install the latest pcDataLoader from the v3 branch?**
```bash
pip3 install pcDataLoader
```

**How to update to the latest pcDataLoader from the v3 branch?**\
This works, even when you have a v2 branch installation.
```bash
pip3 install -U pcDataLoader
```


## HowTo - branch generic:

**How to uninstall pcDataLoader from your python3 environment?**
```bash
pip3 uninstall pcDataLoader
```

**How to check for the current installed pcDataLoader version?**
```python3
import pcDataLoader as pc
pc.__version__
```

**How to load the pcDataLoader library?**
```python3
import pcDataLoader as pc
```

**How to use the addition plotting scripts for plotting PhysiCell output?**\
This plotting scripts can be found in the [pcDataLoader/plotting](https://github.com/elmbeech/pcDataLoader/tree/master/pcDataLoader/plotting) directory.\
You can copy these python3 scripts into your PhysiCell `output` directory and run them there.\
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
