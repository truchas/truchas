#!/usr/bin/env python3

#===============================================================================
#
#  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
#  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
#  in the LICENSE file found in the top-level directory of this distribution.
#
#===============================================================================

import struct
import argparse

import scipy as sp

import truchas

def write_gmv(args, truchas_data):
    filename_base = args.base if args.base else args.simulation
    writetype = "wb" if args.binary else "w"
    gmvtype = "ieeei4r8" if args.binary else "ascii"

    # mesh file
    outfile = filename_base + "-mesh.gmv"
    with open(outfile, writetype) as f:
        fw = FileWriter(f, args.binary)

        fw.write("gmvinput", gmvtype)

        # nodes
        fw.write("nodes   ", truchas_data.nnode)
        for x in truchas_data.node_coordinates().transpose():
            fw.write(x)

        # cells
        fw.write("cells   ", truchas_data.ncell)
        cnode = truchas_data.cell_node_map()
        for cn in cnode:
            if cn[0] == cn[1]: # TET
                fw.write("ptet4   ", 4, cn[1:5])
            elif cn[4] == cn[5]: # PYR
                fw.write("ppyrmd5 ", 5, cn[:5])
            elif cn[4] == cn[7]: # WED
                fw.write("pprism6 ", 6, cn[[0,3,4,1,2,5]])
            else: # HEX
                fw.write("phex8   ", 8, cn)

        fw.write("endgmv  ")

    # series data files
    # prepare non-series data
    nsd_view = truchas_data.non_series_data_view()
    nproc = nsd_view["NUMPROCS"][0]
    cellpart = nsd_view["CELLPART"][:][truchas_data._cellmap] \
        if "CELLPART" in nsd_view.keys() \
        else None
    nodepart = nsd_view["NODEPART"][:][truchas_data._nodemap] \
        if "NODEPART" in nsd_view.keys() \
        else None
    blockids, block_names = \
        block_data(nsd_view["BLOCKID"][:][truchas_data._cellmap])

    # write series
    for series_id in range(1,truchas_data.num_series()+1):
        filesuffix = args.start + args.stride*(series_id-1)
        outfile = filename_base + "-data.gmv.{:04d}".format(filesuffix)
        with open(outfile, writetype) as f:
            fw = FileWriter(f, args.binary)

            fw.write("gmvinput", gmvtype)

            fw.write("nodes   fromfile", "\"" + filename_base + "-mesh.gmv\"")
            fw.write("cells   fromfile", "\"" + filename_base + "-mesh.gmv\"")
            fw.write("cycleno ", truchas_data.cycle(series_id))
            fw.write("probtime", truchas_data.time(series_id))

            series_view = truchas_data._series(series_id)

            # blocks
            fw.write("material", len(block_names), 0)
            for name in block_names:
                fw.write(gmv_str(name))
            fw.write(blockids)

            # flags
            if nodepart is not None or cellpart is not None:
                fw.write("flags   ")

                if cellpart is not None:
                    fw.write("cellpart", nproc, 0)
                    for i in range(nproc): fw.write("PE{:04d}  ".format(i+1))
                    fw.write(cellpart)

                if nodepart is not None:
                    fw.write("nodepart", nproc, 1)
                    for i in range(nproc): fw.write("PE{:04d}  ".format(i+1))
                    fw.write(nodepart)

                fw.write("endflag ")

            # scalars
            fw.write("variable")
            for field in truchas_data.field_names():
                x = truchas_data.field(series_id, field)
                if len(x.shape) == 1:
                    fw.write(gmv_str(field), 0)
                    fw.write(x)
            fw.write("endvars ")

            # vectors
            fw.write("vectors ")
            for field in truchas_data.field_names():
                vector = truchas_data.field(series_id, field)
                if len(vector.shape) == 2 and vector.shape[1] == 3:
                    fw.write(gmv_str(field), 0, 3, 1)
                    fw.write(gmv_str(series_view[field].attrs["FIELDNAME1"]))
                    fw.write(gmv_str(series_view[field].attrs["FIELDNAME2"]))
                    fw.write(gmv_str(series_view[field].attrs["FIELDNAME3"]))
                    for x in vector.transpose():
                        fw.write(x)
            fw.write("endvect ")

            fw.write("endgmv  ")


def block_data(x):
    """Given a block id array, returns a "GMV-corrected" list with ,
    and a list of names."""
    present_blocks = sp.unique(x)

    # Convert to GMV-style blockids. This transforms, for instance,
    # blocks 1, 2, 7, 10 -> 1, 2, 3, 4
    gmv_map = sp.empty(present_blocks[-1]+1, dtype=present_blocks.dtype)
    gmv_map[present_blocks] = sp.arange(1, present_blocks.size+1)
    x_gmv = gmv_map[x]

    return x_gmv, ["blk_" + str(i) for i in present_blocks]


def gmv_str(field):
    """Converts a generic truchas field name to a GMV-compatible string. Strips
    the Truchas 'Z_' prefix if applicable, and pads or trims to return
    8-character string. Accepts strings and byte-like objects."""

    name = field.decode() if type(field) == bytes or type(field) == sp.bytes_ \
        else field
    if name[:2] == "Z_": name = name[2:]
    name = "{:8s}".format(name)[:8]
    return name


class FileWriter:
    """
    This is a class to encapsulate logic for writing to either a binary
    or an ascii file. It is initialized with the file and a bool for whether or
    not a binary file is used. Then the 'write' method is called with any number
    of integers, floats, strings, or Numpy arrays as arguments.

    This function will write a single line, followed by a newline, when an ascii
    format is used. It will also insert spaces between each of the arguments.

    When a binary output file is used, it expects strings with length a multiple
    of 8, corresponding to the GMV format. It does no special handling as when
    ascii is selected.
    """

    def __init__(self, f, is_binary):
        self._f = f
        self._is_binary = is_binary

        self._write = {}
        if is_binary:
            self._write[int] = lambda x: self._f.write(struct.pack('<i', x))
            self._write[float] = lambda x: self._f.write(struct.pack('<d', x))
            self._write[str] = lambda x: self._f.write(x.encode())
            self._write[sp.ndarray] = lambda x: x.tofile(self._f)
            self._write[sp.int32] = self._write[int]
            self._write[sp.float64] = self._write[float]
        else:
            self._write[int] = lambda x: self._f.write(str(x) + " ")
            self._write[str] = lambda x: self._f.write(x + " ")
            self._write[sp.ndarray] = lambda x: x.tofile(self._f, sep=" ")
            self._write[float] = self._write[int]
        self._write[sp.int32] = self._write[int]
        self._write[sp.float64] = self._write[float]

    def write(self, *args):
        for x in args:
            self._write[type(x)](x)
        if not self._is_binary:
            self._f.write("\n")


def parse_argv():
    parser = argparse.ArgumentParser(description="Truchas GMV Writer")
    parser.add_argument("-s", "--simulation", type=str, default="MAIN",
                        metavar="STRING", help="Simulation name.")
    parser.add_argument("-m", "--mesh", type=str, default="DEFAULT",
                        metavar="STRING", help="Mesh name.")
    parser.add_argument("-b", "--binary", action="store_true",
                        help="Write GMV files in a binary format.")
    parser.add_argument("--base", type=str, default=None, metavar="STRING",
                        help=("Use this string instead of the simulation "
                              "name for filenames."))
    parser.add_argument("--stride", type=int, default=1, metavar="NUM",
                        help="GMV data file suffix stride.")
    parser.add_argument("--start", type=int, default=0, metavar="NUM",
                        help="GMV data file initial suffix start value.")
    parser.add_argument("h5file", help="Input H5 file.")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_argv()
    truchas_data = truchas.TruchasData(args.h5file)
    write_gmv(args, truchas_data)
