from pyMCDS import pyMCDS

mcds1 = pyMCDS('output00000001.xml', 'timeseries_set')
mcds2 = pyMCDS('output00000008.xml', 'timeseries_set')

print(mcds1.get_time())
print(mcds2.get_menv_species_list())
print(mcds2.get_concentrations('quorum'))