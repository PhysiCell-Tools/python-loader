# PhysiCell Data Loader Tutorial: pcdl and Python and Pandas

In pcdl output form
[TimeStep](https://github.com/elmbeech/physicelldataloader/blob/master/man/TUTORIAL_python3_timestep.md)
and [TimeSeries](https://github.com/elmbeech/physicelldataloader/blob/master/man/TUTORIAL_python3_timeseries.md)
**mcds.get_conc_df** and **mcds.get_cell_df** can store,
and from the [pcdl command line command](https://github.com/elmbeech/physicelldataloader/blob/master/man/TUTORIAL_commandline.md)
**pcdl_get_conc_df**, **pcdl_get_cell_df**, and **pcdl_get_unit_dict** store data into a [csv](https://en.wikipedia.org/wiki/Comma-separated_values) file.

[Pandas](https://pandas.pydata.org/) is the python3 library you will probably spend a lot of the time when you are analysis your data.
Pandas mimics the computer language [R](https://en.wikipedia.org/wiki/R_(programming_language)).
Pandas provides us with the spreadsheet like DataFrame and Series data types.
Csv files can easily be loaded as a DataFrames or Series.

For getting started with pandas, getting an idea of what it is all about and capable of,
I recommend working through the pandas cookbook from Julia Evens:
+ https://jvns.ca/blog/2013/12/22/cooking-with-pandas/

The pandas library as such is very well documented.
Please spend some time to familiarize yourself with the homepage, the users guides, and the API reference.
+ https://pandas.pydata.org/
+ http://pandas.pydata.org/pandas-docs/stable/user_guide/index.html
+ http://pandas.pydata.org/pandas-docs/stable/reference/index.html


## Dump pcdl data construct from the command line into a csv file

```bash
pcdl_get_conc_df output 2
```
```bash
pcdl_get_cell_df output 2
```


## Dump pcdl data construct from the command line into a csv file

```python
import pcdl

mcdsts = pcdl.TimeSeries('output/')
df_conc = mcdsts.get_conc_df(values=2)
df_conc.to_csv('output/timeseries_conc.csv')
```
```python
import pcdl

mcdsts = pcdl.TimeSeries('output/')
df_cell = mcdsts.get_cell_df(values=2)
df_cell.to_csv('output/timeseries_cell.csv')
```


## Load csv files into python

```python
import pandas as pd

df_conc = pd.read_csv('output/timeseries_conc.csv')
df_conc.info()
```

```python
import pandas as pd

df_cell = pd.read_csv('output/timeseries_cell.csv')
df_cell.info()
```

## Plot data from panadas dataframe

It is easy  to plot data from pandas dataframes.
Have a look at [TUTORIAL_python3_matplotlib.md](https://github.com/elmbeech/physicelldataloader/blob/master/man/TUTORIAL_python3_matplotlib.md).

That's it. The rest is analysis!
