try:
    from .TruchasConfigBuild import *
except ImportError:
    from .TruchasConfigInstall import *

# TruchasMappedData must be imported first, to ensure Truchas's desired HDF5 is
# loaded prior to h5py. Some systems compile h5py against a serial version of
# HDF5, which will prevent Truchas's parallel version from loading and result in
# unresolved symbol errors.
#
# Another issue may occur if the HDF5 libraries used by Truchas (grid_mapping)
# and h5py have different versions (e.g. 1.10.1 and 1.10.2), but the same major
# library version (e.g., libhdf5.so.101.0 and libhdf5.so.101.1). If this
# happens, HDF5 will print an error message and kill the process. For now, there
# are a few ways around this:
#
# 1. Ensure Truchas and h5py use the exact same version of HDF5.
# 2. Ensure Truchas and h5py's HDF5 dependencies have different major versions
#    of the shared object file.
# 3. Install h5py via pip. Here, h5py is compiled against HDF5 libraries with
#    hashed names, which sidesteps the issue.
from .TruchasMappedData import *
from .TruchasData import *
from .TruchasEnvironment import *
from .TruchasTest import *
from .TruchasStudy import *
from .TruchasDatabase import *
