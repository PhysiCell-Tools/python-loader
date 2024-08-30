# PhysiCell Data Loader Tutorial: pcdl and Python and Pandas


# PhysiCell Data Analyis with Pandas

Pandas is the python3 library you will probably spend a lot of the time when you are analysis you data.
Pandas mimics the computer language R.
Pandas provides us with the Series and Dataframe data types.

For getting started with pandas, getiting an idea what it is all about and capable of,
I recommend to work through this pandas cookbook from Julia Evens:
+ https://jvns.ca/blog/2013/12/22/cooking-with-pandas/

The pandas library as such is very well documented.
Please spend some time to familarize your self with the homepage, the users guides, and the API reference.
+ https://pandas.pydata.org/
+ http://pandas.pydata.org/pandas-docs/stable/user_guide/index.html
+ http://pandas.pydata.org/pandas-docs/stable/reference/index.html


Load libraries.

```python
import pandas as pd
import pcdl
print('pcdl version:', pcdl.__version)
```

Load a mcds time step.

```python
mcds = pcdl.TimeStep('output/output00000000.xml')
```

Extract the cell and conc dataframe.

```python
df_conc = mcds.get_conc_df()
df_conc.info()
```
```python
df_cell = mcds.get_cell_df()
df_cell.info()
```

That's it. The rest is analysis.
