#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class DS9(TruchasTest.GoldenTestCase):

  test_name = 'ds9'
  num_procs = 4 # with a parallel executable
  
  # Override the default setUp, omitting the opening of the golden output
  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass() # This runs Truchas
    self.test_output = Truchas.TruchasOutput(self.get_output_file())

  def test_final_temp(self):
    '''DS9: verifying the final temperature field'''
    #FAILtol = 1.0e-10
    tol = 5.0e-4
    mesh=self.test_output.get_mesh()
    centroids=mesh.centroids()
    x=centroids[:,0]
    y=centroids[:,1]
    z=centroids[:,2]
    Tref = 9.0 + 6*x*y -x*x - y*y
    T = self.test_output.get_simulation().get_last_series().get_data('Z_TEMP')
    error = max(abs(T -Tref)/Tref)
    print 'error=%1.9e tol=%1.9e\n' %(error,tol)
    self.assertTrue( error <= tol )

if __name__ == '__main__':
  import unittest
  unittest.main()

