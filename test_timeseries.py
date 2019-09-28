from pyMCDS import *
from pyMCDS_timeseries import *

mcds = pyMCDS_timeseries('timeseries_set')

chem_list = mcds.timeseries[0].get_menv_species_list()
mcds.plot_menv_total(chem_list)
mcds.plot_cell_type_counts()        
