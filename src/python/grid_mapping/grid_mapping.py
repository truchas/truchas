#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

# This file provides a Python interface to the Truchas grid mapper.

import ctypes
import collections

import scipy as sp

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
                ("mesh_map", ctypes.c_void_p)]

libgridmap.mesh_map.argtypes = [ctypes.c_int, ctypes.c_int,
                                ctypes.POINTER(ctypes.c_int),
                                ctypes.POINTER(ctypes.c_int),
                                ctypes.POINTER(ctypes.c_double),
                                ctypes.c_char_p,
                                ctypes.c_double]
libgridmap.mesh_map.restype = _MapData

libgridmap.map_cell_field_c.argtypes = [ctypes.c_void_p, ctypes.c_int,
                                        ctypes.c_int, ctypes.c_int,
                                        ctypes.c_double,
                                        ctypes.POINTER(ctypes.c_double),
                                        ctypes.POINTER(ctypes.c_double)]
libgridmap.map_cell_field_c.restype = None


# Python interfaces
MapData = collections.namedtuple("MapData", ("ncell",
                                             "nnode",
                                             "node_coordinates",
                                             "cnode",
                                             "blockid",
                                             "mesh_map"))

def mesh_map(nnode, ncell, cnode, blockid, coords, filename, scale_factor):
    # Copy data to contiguous arrays and cast as needed
    cnodec = sp.ascontiguousarray(cnode)
    cnodec = cnodec.astype(sp.int32, casting="same_kind", copy=False)

    blockidc = sp.ascontiguousarray(blockid)
    blockidc = blockidc.astype(sp.int32, casting="same_kind", copy=False)

    coordsc = sp.ascontiguousarray(coords)
    coordsc = coordsc.astype(sp.float64, casting="same_kind", copy=False)

    # send back data on the mesh mapping and the new mesh
    cdata = libgridmap.mesh_map(nnode, ncell,
                                cnodec.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                                blockidc.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                                coordsc.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                bytes(filename, "utf-8"), scale_factor)

    # pull data out of the addresses indicated, and put in scipy arrays
    return MapData(cdata.ncell,
                   cdata.nnode,
                   sp.ctypeslib.as_array(cdata.coord, shape=(cdata.nnode,3)),
                   sp.ctypeslib.as_array(cdata.connect, shape=(cdata.ncell,8)),
                   sp.ctypeslib.as_array(cdata.blockid, shape=(cdata.ncell,)),
                   cdata.mesh_map)


def map_cell_field(src, mesh_map, dest_ncell):
    # Copy to contiguous arrays and cast as needed,
    # then send pointer
    srcc = sp.ascontiguousarray(src)
    srcc = srcc.astype(sp.float64, casting="same_kind", copy=False)
    dest = sp.empty(dest_ncell)
    libgridmap.map_cell_field_c(mesh_map, 1, len(srcc), dest_ncell, 0.,
                                srcc.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                dest.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
    return dest
