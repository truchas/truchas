#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

# This file provides a Python interface to the Truchas grid mapper.

import ctypes
import collections

import numpy as np

import truchas

# load library
libgridmap = ctypes.CDLL(truchas.TruchasConfig.libgridmap,
                         mode=ctypes.RTLD_LOCAL)


# Fortran interfaces
class _MapData(ctypes.Structure):
    _fields_ = [("ncell", ctypes.c_int),
                ("nnode", ctypes.c_int),
                ("coord", ctypes.POINTER(ctypes.c_double)),
                ("connect", ctypes.POINTER(ctypes.c_int)),
                ("blockid", ctypes.POINTER(ctypes.c_int)),
                ("mapper", ctypes.c_void_p)]

libgridmap.mapper_init.argtypes = [ctypes.c_int, ctypes.c_int,
                                   ctypes.POINTER(ctypes.c_int),
                                   ctypes.POINTER(ctypes.c_int),
                                   ctypes.POINTER(ctypes.c_double),
                                   ctypes.c_char_p,
                                   ctypes.c_double]
libgridmap.mapper_init.restype = _MapData

libgridmap.map_field.argtypes = [ctypes.c_void_p, ctypes.c_int,
                                 ctypes.c_int, ctypes.c_int,
                                 ctypes.c_double,
                                 ctypes.POINTER(ctypes.c_double),
                                 ctypes.POINTER(ctypes.c_double)]
libgridmap.map_field.restype = None


# Python interfaces
MapData = collections.namedtuple("MapData", ("ncell",
                                             "nnode",
                                             "node_coordinates",
                                             "cnode",
                                             "blockid",
                                             "mapper"))

def mapper_init(nnode, ncell, cnode, blockid, coords, filename, scale_factor):
    # Copy data to contiguous arrays and cast as needed
    cnodec = np.ascontiguousarray(cnode)
    cnodec = cnodec.astype(np.int32, casting="same_kind", copy=False)

    blockidc = np.ascontiguousarray(blockid)
    blockidc = blockidc.astype(np.int32, casting="same_kind", copy=False)

    coordsc = np.ascontiguousarray(coords)
    coordsc = coordsc.astype(np.float64, casting="same_kind", copy=False)

    # send back data on the mesh mapping and the new mesh
    cdata = libgridmap.mapper_init(nnode, ncell,
                                   cnodec.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                                   blockidc.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                                   coordsc.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                   bytes(filename, "utf-8"), scale_factor)

    # pull data out of the addresses indicated, and put in numpy arrays
    return MapData(cdata.ncell,
                   cdata.nnode,
                   np.ctypeslib.as_array(cdata.coord, shape=(cdata.nnode,3)),
                   np.ctypeslib.as_array(cdata.connect, shape=(cdata.ncell,8)),
                   np.ctypeslib.as_array(cdata.blockid, shape=(cdata.ncell,)),
                   cdata.mapper)


def map_field(src, mapper, dest_ncell):
    # Copy to contiguous arrays and cast as needed,
    # then send pointer
    srcc = np.ascontiguousarray(src)
    srcc = srcc.astype(np.float64, casting="same_kind", copy=False)
    dest = np.empty(dest_ncell)
    libgridmap.map_field(mapper, 1, len(srcc), dest_ncell, 0.,
                         srcc.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                         dest.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return dest
