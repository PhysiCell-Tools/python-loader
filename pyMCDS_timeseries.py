import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from pyMCDS import pyMCDS

class pyMCDS_timeseries:
    '''
    This class contains a np.array of pyMCDS objects as well as functions for
    extracting information from that list.

    Parameters
    ----------
    output_path : string
        String containing the path (relative or absolute) to the directory 
        containing the PhysiCell output files

    Attributes
    ----------
    timeseries : array-like (pyMCDS) [n_timesteps,]
        Numpy array of pyMCDS objects sorted by time.
    '''
    def __init__(self, output_path='.'):
        # get a generator of output xml files sorted alphanumerically
        sorted_files = sorted(Path(output_path).glob('output*.xml'))

        # make a big array to keep them all, create a pyMCDS object with each
        # self.timeseries = np.empty(len(sorted_files), dtype=pyMCDS)
        # first we'll do a list, if its slow we'll figure out the array
        ts_list = []
        for f in sorted_files:
            ts_list.append(pyMCDS(f.name, output_path))

        self.timeseries = np.array(ts_list)

    def get_times(self):
        '''
        Helper function to easily return the cumulative time at each step.
        
        Returns
        -------
        time_array : ndaray (np.float) [n_timesteps,]
            Contains the 'current time' at each recorded point.
        '''
        time_array = np.zeros(self.timeseries.shape[0])
        for i in range(len(self.timeseries)):
            time_array[i] = self.timeseries[i].get_time()
        
        return time_array

    def plot_menv_total(self, species):
        '''
        Plot total concentration of chemicals in the system over time.

        Parameters
        ----------
        species : array-like (string) or string
            Name os names of chemical species to plot
        '''
        # set up arrays, quicker than lists
        times = self.get_times()
        total_conc = np.zeros((times.shape[0], len(species)))
        # iterate over the times to fill array
        for time in range(times.shape[0]):
            # iterate over the chemical species to fill array
            for i in range(len(species)):
                conc_arr = self.timeseries[time].get_concentrations(species[i])
                total_conc[time, i] = np.sum(conc_arr)
        # establish plots
        fig, ax = plt.subplots(figsize=(10, 8))
        # add plot elements and show plot
        for i, chem in enumerate(species):
            ax.plot(times, total_conc[:, i], label=chem)
        ax.set_xlabel('Time (min)')
        ax.set_ylabel('Total concentration')
        ax.legend()
        plt.show()
    
    def plot_cell_type_counts(self):
        '''
        Plot the absolute counts of different types of cells over time.
        '''
        times = self.get_times()
        unique_types = np.unique(self.timeseries[0].data['discrete_cells']['cell_type'])
        cell_counts = np.zeros((times.shape[0], unique_types.shape[0]))
        for time_i in range(times.shape[0]):
            df = self.timeseries[time_i].get_cell_df()
            counts = df['cell_type'].value_counts()
            cell_counts[time_i, :] = counts
        fig, ax = plt.subplots(figsize=(10, 8))
        for i in range(cell_counts.shape[1]):
            ax.scatter(times, cell_counts[:, i], label=i, s=10)
        ax.set_xlabel('Time (min)')
        ax.set_ylabel('# of cells in system')
        ax.legend()
        plt.show()


