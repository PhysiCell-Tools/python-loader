try:
    from pcdl.data_timeseries import install_data, uninstall_data
except ModuleNotFoundError:
    print('Warning @ pcdl/__init__.py : you are running a lightweight installation of pcdl.\nTo run pcdl.install_data() or pcdl.uninstall_data() do: pip3 install -U pcdl[data] or pip3 install -U pcdl[all].')
try:
    from pcdl.pyAnnData import TimeStep, TimeSeries, scaler
except ModuleNotFoundError:
    print('Warning @ pcdl/__init__.py : you are running a lightweight installation of pcdl.\nTo run mcds/mcdsts.get_anndata() do: pip3 install -U pcdl[scverse] or pip3 install -U pcdl[all].')
from pcdl.pyMCDS import pyMCDS, graphfile_parser
from pcdl.pyMCDSts import pyMCDSts, make_gif, make_movie
from pcdl.VERSION import __version__
