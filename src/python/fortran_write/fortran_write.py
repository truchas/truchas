#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

import ctypes
import os

import numpy as np

import truchas

# load library
_libfwrite = None
def _init_libfwrite():
    global _libfwrite
    if _libfwrite: return
    _libfwrite = ctypes.CDLL(truchas.TruchasConfig.libfwrite)

    # fortran interfaces
    _libfwrite.fopen.argtypes = [ctypes.c_char_p]
    _libfwrite.fopen.restype = ctypes.c_int

    _libfwrite.fclose.argtypes = [ctypes.c_int]
    _libfwrite.fclose.restype = None

    _libfwrite.fwrite_int.argtypes = [ctypes.c_int, ctypes.c_int]
    _libfwrite.fwrite_int.restype = None

    _libfwrite.fwrite_r8.argtypes = [ctypes.c_int, ctypes.c_double]
    _libfwrite.fwrite_r8.restype = None

    _libfwrite.fwrite_str.argtypes = [ctypes.c_int, ctypes.c_char_p]
    _libfwrite.fwrite_str.restype = None

    _libfwrite.fwrite_i4x1.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.c_int]
    _libfwrite.fwrite_i4x1.restype = None

    _libfwrite.fwrite_r8x1.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_int]
    _libfwrite.fwrite_r8x1.restype = None

    _libfwrite.fwrite_i8x1.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_int8), ctypes.c_int]
    _libfwrite.fwrite_i8x1.restype = None

# python interface
class FortranWrite:
    def __init__(self, filename):
        _init_libfwrite()
        self._unit = _libfwrite.fopen(bytes(filename, "utf-8"))

    def close(self):
        _libfwrite.fclose(self._unit)
        self._unit = None

    def write_i4x0(self, x):
        _libfwrite.fwrite_int(self._unit, x)

    def write_r8x0(self, x):
        _libfwrite.fwrite_r8(self._unit, x)

    def write_str(self, x):
        _libfwrite.fwrite_str(self._unit, bytes(x, "utf-8"))

    def write_i4x1(self, x):
        y = np.ascontiguousarray(x)
        y = y.astype(np.int32, casting="same_kind", copy=False)
        _libfwrite.fwrite_i4x1(self._unit, y.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), y.size)

    def write_r8x1(self, x):
        y = np.ascontiguousarray(x)
        y = y.astype(np.float64, casting="same_kind", copy=False)
        _libfwrite.fwrite_r8x1(self._unit, y.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), y.size)

    def write_i8x2(self, x):
        y = np.ascontiguousarray(x)
        _libfwrite.fwrite_i8x1(self._unit, y.ctypes.data_as(ctypes.POINTER(ctypes.c_int8)), y.size)
