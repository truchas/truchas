#!/usr/bin/env python

import sys
import os

import unittest
import numpy

import Truchas
import TruchasTest

class StaticWall(TruchasTest.GoldenTestCase):

  test_name = 'static-wall'
  num_procs = 4 # with a parallel executable

  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())

  def get_test_field(self,field,cycle,serialize=True,region=None):
    return self.test_output.get_simulation().find_series(cycle=cycle).get_data(field,serialize,region=region)

  def runTest(self):
    '''Verify zero velocity'''

    tol = 1.0e-11
    
    v = self.get_test_field('Z_VC',cycle=1)
    vmag = numpy.empty(v.shape[0],v.dtype)
    for j in range(vmag.size):
      vmag[j] = (v[j,0]**2 + v[j,1]**2 + v[j,2]**2)**0.5
    error = max(vmag)
    if error > tol:
      print 'velocity max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'velocity max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)
    

if __name__ == '__main__':
  unittest.main()




    




