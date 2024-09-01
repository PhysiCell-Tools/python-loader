# PhysiCell Data Loader Tutorial: pcdl and Python and the Json File Format

[Json](https://www.json.org/json-en.html) is a light weight data exchange file format.

In pcdl output form the [TimeSeries](https://github.com/elmbeech/physicelldataloader/blob/master/man/TUTORIAL_python3_timeseries.md) **mcdsts.get_conc_attribute()** and **mcdsts.get_cell_attribute()** can be dumped into a json files, and output form the [pcdl command line command](https://github.com/elmbeech/physicelldataloader/blob/master/man/TUTORIAL_commandline.md) **pcdl_get_conc_attribute**  and **pcdl_get_cell_attribute** is dumped as json file.

The [json library](https://docs.python.org/3/library/json.html) is a part of core python.

Json is the ideal datafromat for unstructured data constructs that not can be stored in csv file format, like a dictionary of lists.
Please note, python objecta are not per se json comaptible.
For example:
Complex numbers are a standard data type in python, but complex numbers can not be stored in json.
Python dictonary keys can be of almost any data type, but json object keys have to be strings.


### Dump pcdl data construnct from the command line into a json file

```bash
pcdl_get_conc_attribute output 2
```
```bash
pcdl_get_cell_attribute output 2
```

### Dump pcdl data construnct from within python into a json file

```python
import json
import pcdl

mcdsts = pcdl.TimeSeries('output/')
dl_conc = mcdsts.get_conc_attribute(values=2)

fp = open('output/timeseries_conc_attribute_minmax.json', 'w')
json.dump(dl_conc, fp)
f.close()
```
```python
import json
import pcdl

mcdsts = pcdl.TimeSeries('output/')
dl_cell = mcdsts.get_cell_attribute(values=2)

fp = open('output/timeseries_cell_attribute_minmax.json', 'w')
json.dump(dl_cell, fp)
f.close()
```


### Load json files into python

```python
import json

dl_conc = json.load('output/timeseries_conc_attribute_minmax.json')
```
```python
import json

dl_cell = json.load('output/timeseries_cell_attribute_minmax.json')
```
The python object we retieve from this pcdl conc attribute and cell attribute files is a dictionary of lists.


That's it. The rest is analysis!
