#!/usr/bin/env python

import sys
import os

import unittest
import numpy

import Truchas
import TruchasTest

class RemapTest(TruchasTest.BaseTestCase):

  def serialize_data(self,cellmap,data):
    ordered=numpy.zeros(data.shape,data.dtype)
    if data.ndim == 1:
      idx=0
      while idx < data.size:
	o_idx=cellmap[idx]-1
	ordered[o_idx]=data[idx]
	idx=idx+1
    return ordered	
      




  def runTest(self):

    self.test_name = 'remap'

    # Same input file, notice pointing to the restart directory
    self.truchas.input=os.path.join(self.get_input_rootdir('restart'),'ds11.inp')

    # Root directory for output
    output_rootdir=self.get_output_rootdir(self.test_name)

    # Run with 4 processors
    output_4pe=self.build_output_directory(self.test_name,'out4pe')
    h5file_4pe=self.build_output_filename(self.test_name,'out4pe','ds11.h5')
    TruchasTest.verify_directory(output_4pe)
    self.clean_output_directory(self.test_name,'out4pe')

    self.truchas.nprocs = 4
    self.truchas.outdir=output_4pe
    self.truchas.run()

    # Run with 2 processors
    output_2pe=self.build_output_directory(self.test_name,'out2pe')
    h5file_2pe=self.build_output_filename(self.test_name,'out2pe','ds11.h5')
    TruchasTest.verify_directory(output_2pe)
    self.clean_output_directory(self.test_name,'out2pe')

    self.truchas.outdir=output_2pe
    self.truchas.nprocs=2
    self.truchas.run()

    # Compare the final output 
    h5_4pe = Truchas.TruchasOutput(h5file_4pe)
    h5_2pe = Truchas.TruchasOutput(h5file_2pe)

    # Get te cell maps
    cellmap_4pe = h5_4pe.get_simulation().cellmap()
    cellmap_2pe = h5_2pe.get_simulation().cellmap()

    # Grab the last temperature field, and loose tolerance
    cell_field='Z_TEMP'
    tol=1.0e-2
    num_series = len(h5_4pe.get_simulation().series_names)

    temp_4pe = h5_4pe.get_simulation().get_series_data(cell_field,id=num_series)
    temp_2pe = h5_2pe.get_simulation().get_series_data(cell_field,id=num_series)

    # Serialize the arrays
    temp_4pe_serial = self.serialize_data(cellmap_4pe,temp_4pe)
    temp_2pe_serial = self.serialize_data(cellmap_2pe,temp_2pe)

    # Compute the norms
    max = numpy.linalg.norm(temp_4pe_serial-temp_2pe_serial,numpy.inf)
    min = numpy.linalg.norm(temp_4pe_serial-temp_2pe_serial,-numpy.inf)
    l2 = numpy.linalg.norm(temp_4pe_serial-temp_2pe_serial)

    print '%s min %.4e max %.4e L2 %.4e' % (cell_field,min,max,l2)
    msg='Fail to satisfy the tol of %.4f'%tol
    self.failUnless(max < tol, msg)

if __name__ == '__main__':
  unittest.main()




    




