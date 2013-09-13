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

    # Grab the last series for each run
    series_4pe = h5_4pe.get_simulation().get_last_series()
    series_2pe = h5_2pe.get_simulation().get_last_series()

    # Grab the last temperature field, and loose tolerance
    cell_field='Z_TEMP'
    tol=1.0e-2

    # Serialize the arrays
    temp_4pe_serial = series_4pe.get_data(cell_field,serialize=True)
    temp_2pe_serial = series_2pe.get_data(cell_field,serialize=True)

    # Compute the norms
    max = numpy.linalg.norm(temp_4pe_serial-temp_2pe_serial,numpy.inf)
    min = numpy.linalg.norm(temp_4pe_serial-temp_2pe_serial,-numpy.inf)
    l2 = numpy.linalg.norm(temp_4pe_serial-temp_2pe_serial)

    print '%s min %.4e max %.4e L2 %.4e' % (cell_field,min,max,l2)
    msg='Fail to satisfy the tol of %.4f'%tol
    self.failUnless(max < tol, msg)

    cell_field='Grad_T'

    # Serialize the arrays
    gtemp_4pe_serial = series_4pe.get_data(cell_field,serialize=True)
    gtemp_2pe_serial = series_2pe.get_data(cell_field,serialize=True)
    gtemp_2pe_serial_max=numpy.linalg.norm(gtemp_2pe_serial,numpy.inf)

    # Compute the norms
    max = numpy.linalg.norm(gtemp_4pe_serial-gtemp_2pe_serial,numpy.inf)
    min = numpy.linalg.norm(gtemp_4pe_serial-gtemp_2pe_serial,-numpy.inf)
    l2 = numpy.linalg.norm(gtemp_4pe_serial-gtemp_2pe_serial)

    print 'max of 2pe = %1.9e %s min %.4e max %.4e L2 %.4e' % (gtemp_2pe_serial_max,cell_field,min,max,l2)
    msg='Fail (%1.9e) to satisfy the tol of %.4f'%(max/gtemp_2pe_serial_max,tol)
    self.failUnless(max/gtemp_2pe_serial_max < tol, msg)

if __name__ == '__main__':
  unittest.main()




    




