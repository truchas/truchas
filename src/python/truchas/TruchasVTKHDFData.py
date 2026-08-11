#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkCompositeDataSet
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
        """Return a field from a zero-based temporal step as a NumPy array."""
        if association not in ("cell", "point", "field"):
            raise ValueError(f"unknown VTKHDF association {association!r}")

        arrays = []
        for block in self._leaf_blocks(self._block(step)):
            data = {
                "cell": block.GetCellData(),
                "point": block.GetPointData(),
                "field": block.GetFieldData(),
            }[association]
            array = data.GetArray(field_name)
            if array is None:
                raise KeyError(f"field {field_name!r} not found in {self.filename}")
            arrays.append(np.asarray(vtk_to_numpy(array)))
        return np.concatenate(arrays) if len(arrays) > 1 else arrays[0]

    def _select_step(self, step):
        if not 0 <= step < self._num_steps:
            raise IndexError(f"VTKHDF step {step} outside [0, {self._num_steps})")
        self._reader.SetStep(step)
        self._reader.Update()

    def _block(self, step):
        self._select_step(step)
        output = self._reader.GetOutputDataObject(0)
        if not hasattr(output, "GetNumberOfBlocks"):
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
