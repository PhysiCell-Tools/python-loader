import importlib.metadata
from .pyMCDS import pyMCDS
from .pyMCDSts import pyMCDSts

# extract version
try: 
    __version__ = importlib.metadata.version("pcdl")
except importlib.metadata.PackageNotFoundError:
    __version__ = "development"
