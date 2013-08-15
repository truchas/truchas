#!/usr/bin/env python

import sys
import os

import unittest

import Truchas
import TruchasTest

class RestartTest(TruchasTest.BaseTestCase):

  input = 'ds11.inp'
  num_procs = 4

  def test_simple(self):

    # Standard output location
    input_rootdir = TruchasTest.get_test_source_rootdir() + \
	            os.path.sep + 'restart'
    output_rootdir = TruchasTest.get_test_build_rootdir() + \
	             os.path.sep + 'restart'
    
    # Set the number of processors if this a parallel binary
    if self.truchas_is_parallel:
      self.truchas.nprocs = self.num_procs

    # First run
    outdir1=output_rootdir+os.path.sep+'out1'
    TruchasTest.verify_directory(outdir1)

    self.truchas.input = input_rootdir + os.path.sep + self.input
    self.truchas.outdir = outdir1
    self.truchas.h5file = outdir1 + os.path.sep + 'ds11.h5'
    self.truchas.stdout = outdir1 + os.path.sep + 'ds11.tty'
    self.truchas.stderr = outdir1 + os.path.sep + 'ds11.err'
    self.truchas.run()

    # Second run, direct output to another directory 
    outdir2=output_rootdir+os.path.sep+'out2'
    TruchasTest.verify_directory(outdir2)

    # Write a restart file
    restart_cycle=45
    restart_file=outdir2 + \
	         os.path.sep + 'restart.%i'%restart_cycle
    self.truchas.write_restart(restart_cycle,output=restart_file)

    # Now pick up the restart file
    self.truchas.outdir = outdir2
    self.truchas.h5file = outdir2 + os.path.sep + 'ds11.h5'
    self.truchas.stdout = outdir2 + os.path.sep + 'ds11.tty'
    self.truchas.stderr = outdir2 + os.path.sep + 'ds11.err'
    self.truchas.run()

    # Compare the final output 
    h51 = Truchas.TruchasOutput(outdir1 + os.path.sep + 'ds11.h5')
    sim1 = h51.get_simulation()
    h52 = Truchas.TruchasOutput(outdir2 + os.path.sep + 'ds11.h5')
    sim2 = h52.get_simulation()

    num_series1 = len(sim1.series_names)
    num_series2 = len(sim2.series_names)

    series_last1 = sim1.get_series(id=num_series1)
    series_last2 = sim2.get_series(id=num_series2)

    # Check all the dataset fields
    for ds in series_last1.datasets:
      if ds not in series_last2.datasets:
	self.fail('Missing dataset')
      d1 = series_last1.get_data(ds)
      d2 = series_last2.get_data(ds)
      if d1.all() != d2.all():
	self.fail('Mismatched dataset values')








     

 
     

if __name__ == '__main__':
  unittest.main()
 
  
