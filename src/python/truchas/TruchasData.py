#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

import re
import os
import sys

import h5py
import numpy as np
import numpy.linalg as npla

# Try to import through PYTHONPATH, and if not found, look in the truchas
# install directory, set by CMake at configure-time.
try:
    import fortran_write
except ImportError:
    sys.path.append("@TruchasPython_INSTALL_PREFIX@")
    import fortran_write


class TruchasData:
    """A class for reading and interacting with Truchas output data. It provides
    access to fields, and allows users to manipulate them.

    Field manipulation is in-memory, not saved to the H5 file. The Truchas
    output is read-only. Modifying data is intended for modified restarts, and
    is done through the assign_field and assign_value functions. The field
    routine will return a user-modified field, if one has been assigned via one
    of those two methods.

    TruchasData objects may be directly initialized by its initializer.
    TruchasData objects are also provided by :func:`TruchasEnvironment.truchas`
    and :func:`TruchasEnvironment.output`.

    :param filename: Filename of the Truchas h5 output file.
    :type filename: str
    """

    def __init__(self, filename):
        self.mapped = False
        self.filename = filename
        """The filename of the HDF5 Truchas output."""
        self.directory = os.path.dirname(filename)
        self._root = h5py.File(filename, 'r')
        self._centroid = None
        self._volume = None
        self._blockid = None
        self._node_coordinates = None
        self._cnode = None

        # load the cell map, correct for Fortran ordering, and invert
        cm = self._root["Simulations/MAIN/Non-series Data/CELLMAP"][:] - 1
        self._cellmap = np.empty(cm.size, cm.dtype)
        self._cellmap[cm] = np.arange(cm.size)

        # repeat for node map
        nm = self._root["Simulations/MAIN/Non-series Data/NODEMAP"][:] - 1
        self._nodemap = np.empty(nm.size, nm.dtype)
        self._nodemap[nm] = np.arange(nm.size)

        self.ncell = self._cellmap.size
        self.nnode = self._nodemap.size

        # These will not get overwritten if we map the data.
        self._h5_ncell = self.ncell
        self._h5_nnode = self.nnode

        self._modified_fields = {}


    def field(self, series_id, field_name):
        """Return the requested field of the requested series as a ndarray.
        This will read the entire field from disk. If the field has been
        reassigned, returns the user-defined copy.

        :param series_id: The ID of a series in the Truchas output file. Valid
            values are in the range ``[1, TruchasData.num_series()]``.
        :type series_id: int

        :param field_name: The name of a field in the Truchas output file. Valid
            values include ``"Z_TEMP"``, ``"Z_RHO"`` (density),
            ``"Z_ENTHALPY"``, ``"phi1"`` (species 1 concentration), ``"sigma"``
            (stress tensor), ``"epsilon"`` (strain tensor), ``"epstherm"``
            (thermal strain), ``"Rotation"``, ``"Displacement"``, ``"Gap
            Displacement"``, ``"Gap Normal Traction"``, ``"VOF"``, ``"Z_P"``
            (pressure), and ``"Z_VC"`` (velocity). Valid values for the active
            h5 file can be found with the :func:`TruchasData.field_names()`
            function.
        :type field_name: str

        :return: The requested field from the H5 file. If the field has been
            reassigned, returns the user-defined copy.
        :rtype: :class:`numpy.ndarray`
        """
        return self._modified_fields[(series_id, field_name)] \
            if (series_id, field_name) in self._modified_fields \
            else self._field(series_id, field_name)


    def assign_field(self, series_id, field_name, field):
        """Reassign a field with a user-defined copy. Future calls to
        :func:`TruchasData.field` will return this copy, and this copy will
        be used in any restarts.

        Note: VOF is a 2D array. The order of the elements is the same as listed
        in the ``PHYSICS`` namelist ``materials`` parameter and ``MATERIAL``
        namelists ``phases`` parameter.

        :param series_id: The ID of a series in the Truchas output file. Valid
            values are in the range ``[1, TruchasData.num_series()]``.

        :type series_id: int

        :param field_name: The name of a field in the Truchas output file. Valid
            values include ``"Z_TEMP"``, ``"Z_RHO"`` (density),
            ``"Z_ENTHALPY"``, ``"phi1"`` (species 1 concentration), ``"sigma"``
            (stress tensor), ``"epsilon"`` (strain tensor), ``"epstherm"``
            (thermal strain), ``"Rotation"``, ``"Displacement"``, ``"Gap
            Displacement"``, ``"Gap Normal Traction"``, ``"VOF"``, ``"Z_P"``
            (pressure), and ``"Z_VC"`` (velocity). Valid values for the active
            h5 file can be found with the :func:`TruchasData.field_names()`
            function.
        :type field_name: str

        :param field: When assigning to a scalar field (.e.g ``"Z_TEMP"``,
            ``"Z_RHO"``), expected is either a 1D Numpy array with the same
            length as the field being overwritten (``ncell`` or ``nnode``), or a
            number. When assigning to a 2D array (e.g. ``"VOF"`` or
            ``"sigma"``), expected is either a 2D Numpy array with the same
            shape as the field being written, or a 1D Numpy array with the same
            length as the second dimension of the field being overwritten
            (``field.shape[1]``).
        :type field: :class:`numpy.ndarray` or number
        """
        self._modified_fields[(series_id, field_name)] = np.full(original_shape, field)


    def assign_value_block(self, series_id, field_name, blockid, value):
        """Assign a value to a field within an element block.

        :param series_id: The ID of a series in the Truchas output file. Valid
            values are in the range ``[1, TruchasData.num_series()]``.

        :type series_id: int

        :param blockid: The ID of cell block in the Truchas output file.
        :type blockid: int

        :param field_name: The name of a field in the Truchas output file. Valid
            values include ``"Z_TEMP"``, ``"Z_RHO"`` (density),
            ``"Z_ENTHALPY"``, ``"phi1"`` (species 1 concentration), ``"sigma"``
            (stress tensor), ``"epsilon"`` (strain tensor), ``"epstherm"``
            (thermal strain), ``"Rotation"``, ``"Displacement"``, ``"Gap
            Displacement"``, ``"Gap Normal Traction"``, ``"VOF"``, ``"Z_P"``
            (pressure), and ``"Z_VC"`` (velocity). Valid values for the active
            h5 file can be found with the :func:`TruchasData.field_names()`
            function.
        :type field_name: str

        :param field: When assigning to a scalar field (.e.g ``"Z_TEMP"``,
            ``"Z_RHO"``), expected is either a 1D Numpy array with the same
            length as the field being overwritten (``ncell`` or ``nnode``), or a
            number. When assigning to a 2D array (e.g. ``"VOF"`` or
            ``"sigma"``), expected is either a 2D Numpy array with the same
            shape as the field being written, or a 1D Numpy array with the same
            length as the second dimension of the field being overwritten
            (``field.shape[1]``).
        :type field: :class:`numpy.ndarray` or number
        """
        field = self.field(series_id, field_name)
        blockids = self.blockid()
        if ((type(value) == list or type(value) == np.array)
            and len(value) > 1 and len(value) == field.shape[1]):
            for i in range(len(field)):
                if (blockid == blockids[i]):
                    field[i] = value
        else:
            # In this case, the field is either a mesh-wide field
            # or a scalar, and we can do a simple broadcast. It
            # might also be an erroneous input.
            field = np.where(blockids == blockid, value, field)
        self.assign_field(series_id, field_name, field)


    def assign_value_cell(self, series_id, field_name, cellid, value):
        """Assign a value to a field at a given cell. This is equivalent to::

            f = TruchasData.output(sid, "xxx")
            f[cellid] = x
            TruchasData.assign_field(sid, "xxx", f)

        :param series_id: The ID of a series in the Truchas output file. Valid
            values are in the range ``[1, TruchasData.num_series()]``.

        :type series_id: int

        :param blockid: The ID of cell block in the Truchas output file.
        :type blockid: int

        :param field_name: The name of a field in the Truchas output file. Valid
            values include ``"Z_TEMP"``, ``"Z_RHO"`` (density),
            ``"Z_ENTHALPY"``, ``"phi1"`` (species 1 concentration), ``"sigma"``
            (stress tensor), ``"epsilon"`` (strain tensor), ``"epstherm"``
            (thermal strain), ``"Rotation"``, ``"Displacement"``, ``"Gap
            Displacement"``, ``"Gap Normal Traction"``, ``"VOF"``, ``"Z_P"``
            (pressure), and ``"Z_VC"`` (velocity). Valid values for the active
            h5 file can be found with the :func:`TruchasData.field_names()`
            function.
        :type field_name: str

        :param field: When assigning to a scalar field (.e.g ``"Z_TEMP"``,
            ``"Z_RHO"``), expected is a scalar number. When assigning to a 2D
            array (e.g. ``"VOF"`` or ``"sigma"``), expected is a 1D Numpy array
            with the same length as the second dimension of the field being
            overwritten (``field.shape[1]``).

        :type field: :class:`numpy.ndarray` or number

        """
        # WARN: Need to confirm this is assigns to the same cell ID
        #       shown in Paraview, and do a mapping if not.
        field = self.field(series_id, field_name)
        field[cellid] = value
        self.assign_field(series_id, field_name, field)


    def _field(self, series_id, field_name):
        """Return the requested field of the requested series as a ndarray.
        This will read the entire field from disk."""
        length = self._series(series_id)[field_name].shape[0]
        if length == self._h5_ncell:
            field = np.array(self._series(series_id)[field_name])[self._cellmap]
        elif length == self._h5_nnode:
            field = np.array(self._series(series_id)[field_name])[self._nodemap]
        else:
            raise NotImplementedError(("Reading fields which are neither "
                                       "cell-centered nor node-centered is not "
                                       "yet supported"))
        return field


    def probe_data(self, probe_filename):
        """Read data from a probe file.

        :param probe_filename: Name of a probe file. Expected to be
            relative to the directory holding the h5 file.
        :type probe_filename: str

        :return: 2D array containing the text file data.
        :rtype: :class:`numpy.ndarray`
        """
        abs_filename = os.path.join(self.directory, probe_filename)
        return np.loadtxt(abs_filename)


    def non_series_data_view(self):
        return self._root["Simulations/MAIN/Non-series Data"]


    def _cell_node_map(self):
        """Load the cnode map, convert to 0-based indexing,
        then reorder for cell and node ordering."""
        if self._cnode is None:
            self._cnode = self._root["Meshes/DEFAULT/Element Connectivity"][:] - 1
            nm = self._root["Simulations/MAIN/Non-series Data/NODEMAP"][:] - 1
            self._cnode = nm[self._cnode[self._cellmap,:]]
        return self._cnode


    def time(self, series_id):
        """Return the time at a given series number.

        :param series_id: The ID of a series in the Truchas output file. Valid
            values are in the range ``[1, TruchasData.num_series()]``.

        :type series_id: int

        :return: The time at the given series_id
        :rtype: number
        """
        return self._series(series_id).attrs["time"]


    def time_step(self, series_id):
        """Return the time step size at a given series number.

        :param series_id: The ID of a series in the Truchas output file. Valid
            values are in the range ``[1, TruchasData.num_series()]``.

        :type series_id: int

        :return: The time step size at the given series_id
        :rtype: number
        """
        return self._series(series_id).attrs["time step"]


    def cycle(self, series_id):
        """Return the cycle number of a given series.

        :param series_id: The ID of a series in the Truchas output file. Valid
            values are in the range ``[1, TruchasData.num_series()]``.

        :type series_id: int

        :return: The cycle number at the given series_id
        :rtype: int
        """
        return self._series(series_id).attrs["cycle"]


    def series_id(self, cycle):
        """Return the series id of a given cycle. If not found, returns None.

        :param cycle: The number of a Truchas cycle.
        :type cycle: int

        :return: The series ID matching a given cycle. If there is no series
            exactly matching the cycle, returns None.

        :rtype: int or None
        """
        nseries = self.num_series()
        for sid in range(1,nseries+1):
            if cycle == self.cycle(sid):
                return sid
        return None


    def centroids(self):
        """
        :return: A ``[ncell,3]`` field of cell centroids.
        :rtype: :class:`numpy.ndarray`
        """
        if self._centroid is None:
            x = self.node_coordinates()
            cnode = self._cell_node_map()
            self._centroid = np.array([np.average(x[_cell_nodes(cn),:], axis=0) for cn in cnode])
        return self._centroid


    def volumes(self):
        """
        :return: A length ``ncell`` field of cell volumes.
        :rtype: :class:`numpy.ndarray`
        """
        if self._volume is None:
            x = self.node_coordinates()
            cnode = self._cell_node_map()
            self._volume = np.array([_cell_volume(x[_cell_nodes(cn),:]) for cn in cnode])
        return self._volume


    def node_coordinates(self):
        """
        :return: A ``[nnode,3]`` field of node coordinates.
        :rtype: :class:`numpy.ndarray`
        """
        if self._node_coordinates is None:
            self._node_coordinates = \
                self._root["Meshes/DEFAULT/Nodal Coordinates"][:][self._nodemap]
        return self._node_coordinates


    def region(self, *region_blockids):
        """Return a list of bools indicating which cells are part of a given
        set of block ids. Input is an arbitrary number of block id integers.

        :param region_blockids: An arbitrary number of block id integers. This
            function takes an arbitrary number of arguments, not a list. E.g.,

            .. code-block:: python

                TruchasData.region(1)
                TruchasData.region(1, 2, 3)

        :type region_blockids: int

        :return: A list of bools indicating which cells are part of a given
            set of block ids. E.g., ``[True, False, False, True, ...]``. This is
            useful for masking off operations, for example:

            .. code-block:: python

                mask = TruchasData.region(1, 2)
                f1[mask] = f2[mask]
                f1[mask] = 0

        :rtype: :class:`numpy.ndarray` of bools
        """
        cell_in_region = np.array([bid in region_blockids
                                   for bid in self.blockid()])
        # for c,cn in zip(cell_in_region,cnode):
        #     if c: node_in_region[cn] = True
        return cell_in_region


    def region_node(self, *region_blockids):
        """Return a list of bools indicating which nodes are part of a given
        set of block ids. Input is an arbitrary number of block id integers.
        A node is considered part of the given block set if any one of the cells
        adjacent to it is part of one block in that set.

        :param region_blockids: An arbitrary number of block id integers. This
            function takes an arbitrary number of arguments, not a list. E.g.,

            .. code-block:: python

                TruchasData.region(1)
                TruchasData.region(1, 2, 3)

        :type region_blockids: int

        :return: A list of bools indicating which nodes are part of a given
            set of block ids. E.g., ``[True, False, False, True, ...]``. This is
            useful for masking off operations, for example:

            .. code-block:: python

                mask = TruchasData.region(1, 2)
                f1[mask] = f2[mask]
                f1[mask] = 0

        :rtype: :class:`numpy.ndarray` of bools
        """
        region_cell = self.region(*region_blockids)
        cnode = self._cell_node_map()

        region_node = np.zeros(self.nnode, dtype=bool)
        for cn, rc in zip(cnode, region_cell):
            if rc:
                region_node[cn] = True

        return region_node


    def blockid(self):
        """
        :return: A ncell length field of integers, holding the Exodus block ID
            at for each cell.
        :rtype: :class:`numpy.ndarray`
        """
        if self._blockid is None:
            self._blockid = self._root["Simulations/MAIN/Non-series Data/BLOCKID"][:][self._cellmap] \
                if "BLOCKID" in self._root["Simulations/MAIN/Non-series Data/"] \
                else []
        return self._blockid


    def field_names(self, series_id=1):
        """Return a list of field names present in the first series.


        :param series_id: The ID of a series in the Truchas output file. Valid
            values are in the range ``[1, TruchasData.num_series()]``. Defaults
            to 1. Usually the fields present does not change over the course
            of a simulation.
        :type series_id: int, optional

        :result: A list of field names present in the requested series.
        :rtype: list of strings
        """
        return list(self._series(series_id).keys())


    def num_species(self):
        """
        :return: The number of species in this file.
        :rtype: int
        """
        return self._root["Simulations/MAIN"].attrs["NUM_SPECIES"]


    def num_series(self):
        """
        :return: The number of series in this file.
        :rtype: int
        """
        return len(self._root["Simulations/MAIN/Series Data/"])


    def _series(self, n):
        """Return the requested series as an h5 type."""
        return self._root["Simulations/MAIN/Series Data/Series {:d}".format(n)]


    def write_restart(self, outfile, series_id):
        """Write a restart file from given Truchas data object and series ID.

        :param outfile: The name of the restart file to be generated.
        :type outfile: str

        :param series_id: The ID of a series in the Truchas output file. Valid
            values are in the range ``[1, TruchasData.num_series()]``.

        :type series_id: int
        """

        fw = fortran_write.FortranWrite(outfile)
        fields = self.field_names(series_id)

        # HEADER SEGMENT
        # file format magic number
        fw.write_str("{:8s}".format("TRF-3"))

        # feature list
        features = self._feature_list(fields)
        fw.write_i4x0(len(features))
        for f in features:
            fw.write_str("{:32s}".format(f))

        # simulation specification (free use -- not using any now)
        fw.write_i4x0(0)

        # global data
        fw.write_r8x0(self.time(series_id))
        fw.write_r8x0(self.time_step(series_id))
        fw.write_i4x0(self.cycle(series_id))
        fw.write_i4x0(self.ncell)
        fw.write_i4x0(self.nnode)

        # MESH SEGMENT
        for cn in self._cell_node_map().transpose():
            fw.write_i4x1(cn+1)

        bid = self.blockid()
        if len(bid) > 0:
            fw.write_i4x0(1)
            fw.write_i4x1(bid)
        else:
            fw.write_i4x0(0)

        for x in self.node_coordinates().transpose():
            fw.write_r8x1(x)

        # CORE DATA SEGMENT
        fw.write_r8x1(self.field(series_id, "Z_RHO"))
        fw.write_r8x1(self.field(series_id, "Z_TEMP"))
        fw.write_r8x1(self.field(series_id, "Z_ENTHALPY"))

        if "fluid_flow" in features:
            fw.write_r8x1(self.field(series_id, "Z_P"))

            for v in self.field(series_id, "Z_VC").transpose():
                fw.write_r8x1(v)

            for v in self.field(series_id, "Face_Vel").transpose():
                fw.write_r8x1(v)
        else:
            dummy = np.zeros(self.ncell)
            for i in range(10):
                fw.write_r8x1(dummy)

        if "VOF" in fields:
            vof = self.field(series_id, "VOF").transpose()
            fw.write_i4x0(vof.shape[0])
            for vf in vof:
                fw.write_r8x1(vf)
        else:
            # single phase problem
            fw.write_i4x0(1)
            fw.write_r8x1(np.ones(self.ncell))

        # LEGACY SOLID MECHANICS SEGMENT
        if "legacy_solid_mechanics" in features:
            for n in range(12):
                for field in ("TOTAL_STRAIN_", "ELASTIC_STRESS_", "PLASTIC_STRAIN_"):
                    name = field + "{:02d}".format(n+1)
                    for t in self.field(series_id, name).transpose():
                        fw.write_r8x1(t)

                name = "PLASTIC_STRAIN_RATE_{:02d}".format(n+1)
                fw.write_r8x1(self.field(series_id, name))

            for t in self.field(series_id, "epsilon").transpose():
                fw.write_r8x1(t)

            for t in self.field(series_id, "sigma").transpose():
                fw.write_r8x1(t)

            for t in self.field(series_id, "e_plastic").transpose():
                fw.write_r8x1(t)

            fw.write_r8x1(self.field(series_id, "epsdot"))

            for t in self.field(series_id, "RHS").transpose():
                fw.write_r8x1(t)

            for t in self.field(series_id, "epstherm").transpose():
                fw.write_r8x1(t)

            for t in self.field(series_id, "epspc").transpose():
                fw.write_r8x1(t)

            for t in self.field(series_id, "Displacement").transpose():
                fw.write_r8x1(t)

        # SPECIES SEGMENT
        if "species" in features:
            nspecies = self.num_species()
            assert nspecies > 0
            fw.write_i4x0(nspecies)

            for n in range(nspecies):
                fw.write_r8x1(self.field(series_id, "phi" + str(n+1)))

        # MICROSTRUCTURE SEGMENT
        if "microstructure" in features:
            ustruc_map = self._series(series_id)["CP-USTRUC-MAP"][:]
            fw.write_i4x0(ustruc_map.size)
            fw.write_i4x1(ustruc_map)

            # get the number of components
            ncomp = 0
            while "CP-USTRUC-COMP-" + str(ncomp+1) in fields: ncomp += 1
            fw.write_i4x0(ncomp)

            for n in range(ncomp):
                name = "CP-USTRUC-COMP-" + str(n+1)
                ustruc_comp = self._series(series_id)[name][:]
                cid = self._series(series_id)[name].attrs["COMP-ID"]

                fw.write_i4x0(cid)
                fw.write_i4x0(ustruc_comp.shape[1])
                fw.write_i8x2(ustruc_comp)

        # JOULE HEAT SEGMENT
        if "joule_heat" in features:
            t = self.time(series_id)

            # scan through the EM simulations, in order, looking for
            # the last one whose TIME attribute value is <= the time
            # attribute for the requested series sequence number.
            em_fields = list(self._root["Simulations"].keys())
            n = 0; t_em_next = 0; name_next = ""
            while t_em_next <= t:
                t_em = t_em_next
                name = name_next

                n += 1
                name_next = "EM{:03d}".format(n)
                if name_next not in em_fields: break
                t_em_next = self._root["Simulations/" + name_next] \
                                        .attrs["TIME"]
            em_sim = self._root["Simulations/" + name + "/Non-series Data"]
            print("Using EM simulation {:s} (t = {:f})".format(name, t_em))

            fw.write_r8x0(em_sim["FREQ"][0])
            fw.write_r8x0(em_sim["UHFS"][0])

            coils = em_sim["COILS"][:]
            fw.write_i4x0(coils.shape[0])
            for c in coils:
                fw.write_r8x0(c[0])
                fw.write_r8x1(c[1:4])
                fw.write_r8x0(c[4])
                fw.write_r8x0(c[5])
                fw.write_i4x0(int(round(c[6])))

            array = em_sim["MU"][:]
            fw.write_i4x0(array.size)
            fw.write_r8x1(array)

            array = em_sim["SIGMA"][:]
            fw.write_i4x0(array.size)
            fw.write_r8x1(array)

            array = em_sim["JOULE"][:][self._cellmap]
            fw.write_i4x0(self.ncell)
            fw.write_r8x1(array)

        fw.close()


    def _feature_list(self, fields):
        """Returns a list of features present based on fields in the output."""
        # Only species is supported for mapped restarts.
        features = []
        if "Z_VC" in fields and not self.mapped: features.append("fluid_flow")

        if "TOTAL_STRAIN_01" in fields and not self.mapped:
            features.append("legacy_solid_mechanics")
        elif "epsilon" in fields and not self.mapped:
            # new solid mechanics doesn't yet contribute to the restart file
            features.append("solid_mechanics")

        if "phi1" in fields: features.append("species")
        if "Joule_P" in fields and not self.mapped: features.append("joule_heat")
        if "CP-USTRUC-MAP" in fields and not self.mapped: features.append("microstructure")
        return features


def _cell_nodes(cn):
    """Map node IDs from degenerate hex to tet, pyramid, wedge, or hex."""
    cn2 = (   cn[1:5] if cn[0] == cn[1] # TET
        else  cn[:5]  if cn[4] == cn[5] # PYR
        else  cn[:6]  if cn[4] == cn[7] # PYR
        else  cn)                       # HEX
    return cn2


def _cell_volume(x):
    """Return the volume for a cell, given its node coordinates."""
    vol = (  _tet_volume(x) if 4 == x.shape[0]
        else _pyr_volume(x) if 5 == x.shape[0]
        else _wed_volume(x) if 6 == x.shape[0]
        else _hex_volume(x))
    return vol


def _tet_volume(x):
    return npla.det([x[1,:]-x[0,:], x[2,:]-x[0,:], x[3,:]-x[0,:]]) / 6


def _pyr_volume(x):
    vol = (_tet_volume(x[[0,1,3,4],:])
        +  _tet_volume(x[[1,2,0,4],:])
        +  _tet_volume(x[[2,3,1,4],:])
        +  _tet_volume(x[[3,0,2,4],:])) / 2
    return vol


def _wed_volume(x):
    vol = (_tet_volume(x[[0,1,2,3],:])
        +  _tet_volume(x[[4,3,5,1],:])
        +  _tet_volume(x[[1,2,3,5],:])
        +  _tet_volume(x[[1,2,0,4],:])
        +  _tet_volume(x[[3,5,4,0],:])
        +  _tet_volume(x[[0,2,5,4],:])) / 2
    return vol


def _hex_volume(x):
    vol = (_tet_volume(x[[0,1,3,4],:])
        +  _tet_volume(x[[1,2,0,5],:])
        +  _tet_volume(x[[2,3,1,6],:])
        +  _tet_volume(x[[3,0,2,7],:])
        +  _tet_volume(x[[4,7,5,0],:])
        +  _tet_volume(x[[5,4,6,1],:])
        +  _tet_volume(x[[6,5,7,2],:])
        +  _tet_volume(x[[7,6,4,3],:])
        +  _tet_volume(x[[0,2,7,5],:])
        +  _tet_volume(x[[1,3,4,6],:])) / 2
    return vol
