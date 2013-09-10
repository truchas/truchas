#!/usr/bin/env python

import sys
import os

import unittest

import Truchas
import TruchasTest

class VFRad2(TruchasTest.BaseTestCase):
  
  test_name = 'vfrad2'
  num_procs = 4 # with a parallel executable

  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    if self.truchas_is_parallel:
      self.truchas.nprocs=self.num_procs

  def runTest(self):
    '''Verifying invariant final temperature field with centered and shifted enclosures'''
    
    # First run with centered enclosure.
    infile = 'run-centered'
    self.truchas.input = os.path.join(self.get_input_rootdir(self.test_name),infile+'.inp')
    self.truchas.outdir = self.build_output_directory(self.test_name,infile+'_output')
    self.truchas.h5file = os.path.join(self.truchas.outdir,infile+'.h5')   
    self.truchas.run()
    
    centered_output = Truchas.TruchasOutput(self.truchas.h5file)
    
    # Second run with shifted enclosure.
    infile = 'run-shifted'
    self.truchas.input = os.path.join(self.get_input_rootdir(self.test_name),infile+'.inp')
    self.truchas.outdir = self.build_output_directory(self.test_name,infile+'_output')
    self.truchas.h5file = os.path.join(self.truchas.outdir,infile+'.h5')   
    self.truchas.run()
    
    shifted_output = Truchas.TruchasOutput(self.truchas.h5file)
    
    # Compare the final temperature fields -- should be identical
    Tcent = centered_output.get_simulation().find_series(id=6).get_data('Z_TEMP')
    Tshft =  shifted_output.get_simulation().find_series(id=6).get_data('Z_TEMP')
    error = max(abs(Tcent-Tshft))
    self.assertTrue(error <= 1.0e-13)


if __name__ == '__main__':
  unittest.main()
 
  
