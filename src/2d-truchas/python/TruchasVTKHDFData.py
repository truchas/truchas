#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkCompositeDataSet
from vtkmodules.vtkFiltersParallel import vtkRemoveGhosts
from vtkmodules.vtkIOHDF import vtkHDFReader


class TruchasVTKHDFData:
    """Read temporal fields from a Truchas VTKHDF output file.

    The reader deliberately uses VTK's VTKHDF reader rather than depending on
    the file's internal HDF5 layout.  This supports both the serial and
    parallel VTKHDF layouts.

    :param filename: VTKHDF output filename.
    :param blockname: Optional mesh block name. Defaults to the first block.
    """

    def __init__(self, filename, blockname=None):
        self.filename = str(filename)
        self.blockname = blockname
        self._reader = vtkHDFReader()
        self._reader.SetFileName(self.filename)
        self._reader.UpdateInformation()
        self._num_steps = max(1, self._reader.GetNumberOfSteps())

    @property
    def num_steps(self):
        """Number of temporal steps in the file."""
        return self._num_steps

    def time(self, step):
        """Return the time associated with a zero-based temporal step."""
        self._select_step(step)
        return float(self._reader.GetTimeValue())

    def field(self, step, field_name, association="cell"):
        """Return a field from a zero-based temporal step as a NumPy array.

        Cell and point fields exclude VTK ghost entities and are ordered by
        their pedigree IDs, giving a partition-independent representation.
        """
        if association not in ("cell", "point", "field"):
            raise ValueError(f"unknown VTKHDF association {association!r}")

        if association == "field":
            # UG output stores temporal field data on the root object.  In
            # parallel this is a vtkPartitionedDataSet, while serial output
            # may be the unstructured grid itself.
            data = self._block(step).GetFieldData()
            array = data.GetArray(field_name)
            if array is not None:
                return np.asarray(vtk_to_numpy(array))

            # Retain support for older multiblock files whose field data was
            # exposed on a leaf by VTK.
            for block in self._leaf_blocks(self._block(step)):
                array = block.GetFieldData().GetArray(field_name)
                if array is not None:
                    return np.asarray(vtk_to_numpy(array))
            raise KeyError(f"field {field_name!r} not found in {self.filename}")

        arrays = []
        pedigree_ids = []
        for block in self._leaf_blocks(self._block(step)):
            if association == "cell":
                remove_ghosts = vtkRemoveGhosts()
                remove_ghosts.SetInputData(block)
                remove_ghosts.Update()
                block = remove_ghosts.GetOutput()
            data = {
                "cell": block.GetCellData(),
                "point": block.GetPointData(),
            }[association]
            array = data.GetArray(field_name)
            if array is None:
                raise KeyError(f"field {field_name!r} not found in {self.filename}")
            values = np.asarray(vtk_to_numpy(array))
            if association in ("cell", "point"):
                pedigree = data.GetPedigreeIds()
                if pedigree is None:
                    raise ValueError(
                        f"{association} data in {self.filename} has no pedigree IDs"
                    )
                ids = np.asarray(vtk_to_numpy(pedigree))
                if association == "point":
                    ghost_type = data.GetArray("vtkGhostType")
                    if ghost_type is None:
                        raise ValueError(
                            f"point data in {self.filename} has no vtkGhostType"
                        )
                    on_process = np.asarray(vtk_to_numpy(ghost_type)) == 0
                    values = values[on_process]
                    ids = ids[on_process]
                arrays.append(values)
                pedigree_ids.append(ids)
        result = np.concatenate(arrays) if len(arrays) > 1 else arrays[0]
        pedigree_ids = np.concatenate(pedigree_ids) \
            if len(pedigree_ids) > 1 else pedigree_ids[0]
        result = self._order_by_pedigree(result, pedigree_ids, association)
        return result

    def cell_centers(self, step):
        """Return cell geometric centroids in pedigree-ID order."""
        centers = []
        pedigree_ids = []
        for block in self._leaf_blocks(self._block(step)):
            remove_ghosts = vtkRemoveGhosts()
            remove_ghosts.SetInputData(block)
            remove_ghosts.Update()
            block = remove_ghosts.GetOutput()
            pedigree = block.GetCellData().GetPedigreeIds()
            if pedigree is None:
                raise ValueError(f"cell data in {self.filename} has no pedigree IDs")
            centers.append(self._cell_centers(block))
            pedigree_ids.append(np.asarray(vtk_to_numpy(pedigree)))

        centers = np.concatenate(centers) if len(centers) > 1 else centers[0]
        pedigree_ids = np.concatenate(pedigree_ids) \
            if len(pedigree_ids) > 1 else pedigree_ids[0]
        return self._order_by_pedigree(centers, pedigree_ids, "cell")

    def cell_vertices(self, step):
        """Return cell vertex coordinates in pedigree-ID order."""
        vertices = []
        pedigree_ids = []
        for block in self._leaf_blocks(self._block(step)):
            remove_ghosts = vtkRemoveGhosts()
            remove_ghosts.SetInputData(block)
            remove_ghosts.Update()
            block = remove_ghosts.GetOutput()
            pedigree = block.GetCellData().GetPedigreeIds()
            if pedigree is None:
                raise ValueError(f"cell data in {self.filename} has no pedigree IDs")
            for index in range(block.GetNumberOfCells()):
                points = block.GetCell(index).GetPoints().GetData()
                vertices.append(np.asarray(vtk_to_numpy(points)).copy())
            pedigree_ids.extend(vtk_to_numpy(pedigree))

        pedigree_ids = np.asarray(pedigree_ids)
        if len(np.unique(pedigree_ids)) != len(pedigree_ids):
            raise ValueError(f"cell pedigree IDs in {self.filename} are not unique")
        return [vertices[index] for index in np.argsort(pedigree_ids, kind="stable")]

    def _select_step(self, step):
        if not 0 <= step < self._num_steps:
            raise IndexError(f"VTKHDF step {step} outside [0, {self._num_steps})")
        self._reader.SetStep(step)
        self._reader.Update()

    def _block(self, step):
        self._select_step(step)
        output = self._reader.GetOutputDataObject(0)
        if not hasattr(output, "GetNumberOfBlocks"):
            if self.blockname is not None:
                raise KeyError(
                    f"block {self.blockname!r} requested from output with no named blocks"
                )
            return output

        for index in range(output.GetNumberOfBlocks()):
            block = output.GetBlock(index)
            if block is None:
                continue
            metadata = output.GetMetaData(index)
            name = metadata.Get(vtkCompositeDataSet.NAME()) \
                if metadata is not None and metadata.Has(vtkCompositeDataSet.NAME()) \
                else None
            if self.blockname is None or name == self.blockname:
                return block

        requested = self.blockname or "<first block>"
        raise KeyError(f"block {requested!r} not found in {self.filename}")

    @staticmethod
    def _leaf_blocks(block):
        if hasattr(block, "GetCellData"):
            return [block]
        if hasattr(block, "GetNumberOfPartitions"):
            leaves = []
            for index in range(block.GetNumberOfPartitions()):
                partition = block.GetPartition(index)
                if partition is not None:
                    leaves.extend(TruchasVTKHDFData._leaf_blocks(partition))
            return leaves
        if hasattr(block, "GetNumberOfPieces"):
            leaves = []
            for index in range(block.GetNumberOfPieces()):
                piece = block.GetPiece(index)
                if piece is not None:
                    leaves.extend(TruchasVTKHDFData._leaf_blocks(piece))
            return leaves
        if hasattr(block, "GetNumberOfBlocks"):
            leaves = []
            for index in range(block.GetNumberOfBlocks()):
                child = block.GetBlock(index)
                if child is not None:
                    leaves.extend(TruchasVTKHDFData._leaf_blocks(child))
            return leaves
        raise TypeError(f"unsupported VTKHDF block type {type(block).__name__}")

    def _order_by_pedigree(self, values, pedigree_ids, association):
        if len(np.unique(pedigree_ids)) != len(pedigree_ids):
            raise ValueError(
                f"{association} pedigree IDs in {self.filename} are not unique"
            )
        return values[np.argsort(pedigree_ids, kind="stable")]

    @staticmethod
    def _cell_centers(block):
        centers = np.empty((block.GetNumberOfCells(), 3))
        for index in range(block.GetNumberOfCells()):
            x = np.asarray(vtk_to_numpy(block.GetCell(index).GetPoints().GetData()))
            if len(x) < 3:
                raise ValueError("cannot compute the centroid of a cell with fewer than 3 nodes")
            if len(x) == 3:
                centers[index] = x.mean(axis=0)
                continue

            triangle_area = np.linalg.norm(
                np.cross(x[1:-1] - x[0], x[2:] - x[0]), axis=1
            ) / 2.0
            triangle_centers = (x[0] + x[1:-1] + x[2:]) / 3.0
            centers[index] = np.average(triangle_centers, axis=0, weights=triangle_area)
        return centers
