# pcDataLoader Tutorial Man Page

## Tutorial:
+ http://www.mathcancer.org/blog/python-loader/

```python3
import pcDataLoader as pc

mcds = pc.pyMCDS('data_snapshot/output00003696.xml', microenv=False)
mcds = pc.pyMCDS('data_snapshot/output00003696.xml')
```
```python3
import pcDataLoader as pc

l_mcds = pc.pyMCDS_timeseries('data_snapshot/output00003696.xml', microenv=False)
l_mcds = pc.pyMCDS_timeseries('data_snapshot/output00003696.xml')
```

```python3
import pcDataLoader as pc

mcds.get_substrate_names()
mcds.get_concentrations_df()
mcds.get_concentrations(mcds.get_substrate_names()[0])
mcds.get_concentrations_at(x=0, y=0, z=0)
```

```python3
import pcDataLoader as pc

mcds.get_cell_variables()
mcds.get_cell_df()
mcds.get_cell_df_at(x=0,y=0,z=0)
```
