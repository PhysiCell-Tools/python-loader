try:
    from pcdl.data_timeseries import install_data, uninstall_data
except ModuleNotFoundError:
    print('Warning @ pcdl/__init__.py : you are running a lightweight installation of pcdl.\nTo run pcdl.install_data() or pcdl.uninstall_data() do: pip3 install pcdl[data].')
try:
    from pcdl.pyAnnData import TimeStep, TimeSeries, _anndextract, scaler
except ModuleNotFoundError:
    print('Warning @ pcdl/__init__.py : you are running a lightweight installation of pcdl.\nTo run mcds.get_anndata() or mcdsts.get_anndata() do: pip3 install pcdl[scverse].')
from pcdl.pyMCDS import pyMCDS, graphfile_parser
from pcdl.pyMCDSts import pyMCDSts
from pcdl.VERSION import __version__
