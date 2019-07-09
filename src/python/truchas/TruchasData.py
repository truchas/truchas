#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

import re
import os

import h5py
import scipy as sp

class TruchasData:
    def __init__(self, filename):
        self.mapped = False
        self.filename = filename
        self.directory = os.path.dirname(filename)
        self._root = h5py.File(filename, 'r')
        self._centroid = None

        # load the cell map, correct for Fortran ordering, and invert
        cm = self._root["Simulations/MAIN/Non-series Data/CELLMAP"][:] - 1
        self._cellmap = sp.empty(cm.size, cm.dtype)
        self._cellmap[cm] = sp.arange(cm.size)

        # repeat for node map
        nm = self._root["Simulations/MAIN/Non-series Data/NODEMAP"][:] - 1
        self._nodemap = sp.empty(nm.size, nm.dtype)
        self._nodemap[nm] = sp.arange(nm.size)

        self.ncell = self._cellmap.size
        self.nnode = self._nodemap.size


    def field(self, series_id, field_name):
        """Return the requested field of the requested series as a ndarray.
        This will read the entire field from disk."""
        length = self._series(series_id)[field_name].shape[0]
        if length == self.ncell:
            field = sp.array(self._series(series_id)[field_name])[self._cellmap]
        elif length == self.nnode:
            field = sp.array(self._series(series_id)[field_name])[self._nodemap]
        else:
            print("Field is not either cell centered or node centered. Skip for now.")
            assert False
        return field


    def field_node(self, series_id, field_name):
        """Return the requested field of the requested series as a ndarray.
        Assumed to be node-centered. This will read the entire field from
        disk."""
        return sp.array(self._series(series_id)[field_name])[self._nodemap]


    def probes(self):
        """Return a list of available probes."""
        return list(self._root["Simulations/MAIN/Probes"].keys())


    def probe_view(self, name):
        """Return a probe H5 view from it's identifier."""
        return self._root["Simulations/MAIN/Probes/" + name]


    def non_series_data_view(self):
        return self._root["Simulations/MAIN/Non-series Data"]


    def probe(self, name, field):
        """Return probe data that matches the given name and field."""
        field_regex = r"P\d+:" + field
        probe = next(probe for probe_name, probe in self._root["Simulations/MAIN/Probes"].items() \
                     if re.match(field_regex, probe_name) and name == probe.attrs["NAME"].decode())
        return probe[:]


    def cell_node_map(self):
        # Load the cnode map, convert to 0-based indexing,
        # reorder for cell and node ordering, then convert
        # back to 1-based indexing.
        cnode = self._root["Meshes/DEFAULT/Element Connectivity"][:] - 1
        nm = self._root["Simulations/MAIN/Non-series Data/NODEMAP"][:] - 1
        cnode = nm[cnode[self._cellmap,:]] + 1
        return cnode


    def time(self, series_id):
        """Return the time at a given series number."""
        return self._series(series_id).attrs["time"]


    def time_step(self, series_id):
        """Return the time at a given series number."""
        return self._series(series_id).attrs["time step"]


    def cycle(self, series_id):
        """Return the cycle id of a given series."""
        return self._series(series_id).attrs["cycle"]

    def series_id(self, cycle):
        """Return the series id of a given cycle. If not found, returns None."""
        nseries = self.num_series()
        for sid in range(1,nseries+1):
            if cycle == self.cycle(sid):
                return sid
        return None

    def centroids(self):
        """Return a list of cell centroids"""
        if self._centroid is None:
            node = self._root["Meshes/DEFAULT/Nodal Coordinates"][:]
            cnode = self._root["Meshes/DEFAULT/Element Connectivity"][:] - 1
            self._centroid = sp.zeros((cnode.shape[0],3))
            for i, cn in enumerate(cnode):
                if cn[0] == cn[1]: # TET
                    self._centroid[i,:] = sp.sum(node[cn[1:5]], axis=0) / 4
                elif cn[4] == cn[5]: # PYR
                    self._centroid[i,:] = sp.sum(node[cn[:5]], axis=0) / 5
                elif cn[4] == cn[7]: # WED
                    self._centroid[i,:] = sp.sum(node[cn[:6]], axis=0) / 6
                else: # HEX
                    self._centroid[i,:] = sp.sum(node[cn], axis=0) / 8
            self._centroid = self._centroid[self._cellmap]
        return self._centroid


    def node_coordinates(self):
        """Return a list of node coordinates."""
        return self._root["Meshes/DEFAULT/Nodal Coordinates"][:][self._nodemap]


    def region(self, *region_blockids):
        """Return a list of bools indicating which cells are part of a given
        set of block ids. Input is an arbitrary number of block id integers."""
        cell_in_region = sp.array([bid in region_blockids
                                   for bid in self.blockid()])
        # for c,cn in zip(cell_in_region,cnode):
        #     if c: node_in_region[cn] = True
        return cell_in_region


    def region_node(self, *region_blockids):
        """Return a list of bools indicating which nodes are part of a given
        set of block ids. Input is an arbitrary number of block id integers."""
        region_cell = self.region(*region_blockids)
        cnode = self.cell_node_map() - 1

        region_node = sp.zeros(self.nnode, dtype=bool)
        for cn, rc in zip(cnode, region_cell):
            if rc:
                region_node[cn] = True

        return region_node


    def blockid(self):
        blockid = self._root["Simulations/MAIN/Non-series Data/BLOCKID"][:][self._cellmap] \
            if "BLOCKID" in self._root["Simulations/MAIN/Non-series Data/"] \
            else []
        return blockid


    def field_names(self, series_id=1):
        """Return a list of field names present in the first series."""
        return list(self._series(series_id).keys())


    def num_species(self):
        return self._root["Simulations/MAIN"].attrs["NUM_SPECIES"]


    def num_series(self):
        """Return the number of series in this file."""
        return len(self._root["Simulations/MAIN/Series Data/"])


    def _series(self, n):
        """Return the requested series as an h5 type."""
        return self._root["Simulations/MAIN/Series Data/Series {:d}".format(n)]
