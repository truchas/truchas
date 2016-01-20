#!/usr/bin/env python

import sys
import os
import math

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

  def get_test_field(self,field,cycle,serialize=True,region=None):
    return self.test_output.get_simulation().find_series(cycle=cycle).get_data(field,serialize,region=region)
    
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

    # The centroids function doesn't serialize, so we don't want to here either.
    p = self.get_test_field('Z_P',cycle=n,serialize=False)
    vof = self.get_test_field('VOF', cycle=n,serialize=False)
    xc = self.test_output.get_mesh().centroids()
    
    n = 0
    pavg = 0.0
    zavg = 0.0
    for j in range(p.size):
      if vof[j,0] > 0.99:
        pavg += p[j]
        zavg += xc[j,2]
        n += 1
    pavg = pavg/n
    zavg = zavg/n
    error = 0.0
    dpdz = -9.81e3
    for j in range(p.size):
      if vof[j,0] > 0.99:
        error += ((p[j]-pavg) - dpdz*(xc[j,2]-zavg))**2
    error = math.sqrt(error/n)
    
    tol = 1.0e-4
    if (error > tol):
      print dpdz
      print 'fluid pressure l2 error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'fluid pressure l2 error = %8.2e: PASS (tol=%8.2e)'%(error,tol)
    

if __name__ == '__main__':
  import unittest
  unittest.main()
