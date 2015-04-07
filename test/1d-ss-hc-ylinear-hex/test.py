#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class test(TruchasTest.GoldenTestCase):

  test_name = '1d-ss-hc-ylinear-hex'
  num_procs = 4 # with a parallel executable
  
  # Override the default setUp, omitting the opening of the golden output
  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass() # This runs Truchas
    self.test_output = Truchas.TruchasOutput(self.get_output_file())

  def test_final_fields(self):
    '''verifying the final T and H fields'''
    tol = 1.0e-6
    mesh=self.test_output.get_mesh()
    centroids=mesh.centroids()
    Tref = 1.0 + centroids[:,1]
    T = self.test_output.get_simulation().get_last_series().get_data('Z_TEMP',serialize=False)
    error = max(abs(T -Tref)/Tref)
    print 'T error=%1.9e tol=%1.9e\n' %(error,tol)
    self.assertTrue( error <= tol )
    H = self.test_output.get_simulation().get_last_series().get_data('Z_ENTHALPY',serialize=False)
    error = max(abs(0.5*H -Tref)/Tref)
    print 'H error=%1.9e tol=%1.9e\n' %(error,tol)
    self.assertTrue( error <= tol )

if __name__ == '__main__':
  import unittest
  unittest.main()

