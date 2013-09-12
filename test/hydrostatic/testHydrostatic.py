#!/usr/bin/env python

import sys
import os

import unittest
import numpy

import Truchas
import TruchasTest

class Hydrostatic(TruchasTest.GoldenTestCase):

  test_name = 'hydrostatic'
  num_procs = 4 # with a parallel executable

  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())

  def get_test_field(self,field,cycle,region=None):
    return self.test_output.get_simulation().find_series(cycle=cycle).get_data(field,region)
    
  def test_velocity(self):
    '''Verifying zero velocity'''
    
    n = 20
    tol = 1.0e-9
    fail = 0
    
    test = self.get_test_field('Z_VC',cycle=n)
    for k in range(3):
      error = max(abs(test[:,k]))
      if (error > tol):
        fail += 1
        print '%s-velocity: max abs error = %8.2e: FAIL (tol=%8.2e)'%(('x','y','z')[k],error,tol)
      else:
        print '%s-velocity: max abs error = %8.2e: PASS (tol=%8.2e)'%(('x','y','z')[k],error,tol)
    self.assertTrue(fail == 0)

  def test_pressure(self):
    '''Verify hydrostatic pressure head'''
    
    n = 20
    tol = 4.0e-7

    test = self.get_test_field('Z_P',cycle=n)
    xc = self.test_output.get_mesh().centroids()
    
    # This is really lame ...
    top = 54  # cell 55 (1-based)
    bot =  6  # cell  7 (1-based)
    dpdz = (test[top]-test[bot]) / (xc[top,2]-xc[bot,2])
    
    gold = -9.81e3
    error = abs((dpdz - gold)/gold) 
    if (error > tol):
      print 'pressure gradient rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'pressure gradient rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)
    

if __name__ == '__main__':
  import unittest
  unittest.main()
