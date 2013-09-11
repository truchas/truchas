#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class ShearConstrHeat(TruchasTest.GoldenTestCase):

  test_name = 'shear-constr-heat'
  num_procs = 4 # with a parallel executable
  
  # Override the default setUp, omitting the opening of the golden output
  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass() # This runs Truchas
    self.test_output = Truchas.TruchasOutput(self.get_output_file())

  def get_test_field(self,field,cycle,region=None):
    return self.test_output.get_simulation().find_series(cycle=cycle).get_data(field,region)

  def test_initial_stress(self):
    '''Verify the initial stress field'''
    abs_tol = 1.0e-8
    rel_tol = 1.0e-8
    fail = 0

    test = self.get_test_field('sigma',cycle=0)

    error = max(abs(test[:,0]))
    tol = abs_tol
    if error > tol:
      fail += 1
      print 'sigma_xx: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_xx: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs(test[:,1]))
    tol = abs_tol
    if error > tol:
      fail += 1
      print 'sigma_yy: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_yy: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    gold = -152.533333333
    error = max(abs((test[:,2]-gold)/gold))
    tol = rel_tol
    if error > tol:
      fail += 1
      print 'sigma_zz: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_zz: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs(test[:,3]))
    tol = abs_tol
    if error > tol:
      fail += 1
      print 'sigma_xy: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_xy: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    gold = 500.0
    error = max(abs((test[:,4]-gold)/gold))
    tol = rel_tol
    if error > tol:
      fail += 1
      print 'sigma_xz: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_xz: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs(test[:,5]))
    tol = abs_tol
    if error > tol:
      fail += 1
      print 'sigma_xx: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_xx: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    self.assertTrue(fail == 0)


  def test_initial_strain(self):
    '''Verify the initial strain field'''

    abs_tol = 1.0e-8
    rel_tol = 1.0e-8
    fail = 0

    test = self.get_test_field('epsilon',cycle=0)

    gold = 2.93333333333e-3
    error = max(abs((test[:,0]-gold)/gold))
    tol = rel_tol
    if error > tol:
      fail += 1
      print 'epsilon_xx: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'epsilon_xx: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    gold = 2.93333333333e-3
    error = max(abs((test[:,1]-gold)/gold))
    tol = rel_tol
    if error > tol:
      fail += 1
      print 'epsilon_yy: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'epsilon_yy: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs(test[:,2]))
    tol = abs_tol
    if error > tol:
      fail += 1
      print 'epsilon_zz: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'epsilon_zz: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs(test[:,3]))
    tol = abs_tol
    if error > tol:
      fail += 1
      print 'epsilon_xy: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'epsilon_xy: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    gold = 500.0 / 5.2e4
    error = max(abs((test[:,4]-gold)/gold))
    tol = rel_tol
    if error > tol:
      fail += 1
      print 'epsilon_xz: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'epsilon_xz: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs(test[:,5]))
    tol = abs_tol
    if error > tol:
      fail += 1
      print 'epsilon_xx: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'epsilon_xx: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    self.assertTrue(fail == 0)

  def test_initial_displacement(self):
    '''Verify initial displacement'''
    
    tol = 1.0e-10
    fail = 0
    
    # Analytic displacements
    exx = 2.93333333333e-3
    eyy = 2.93333333333e-3
    exz = 500.0 / 5.2e4
    [x,y,z] = self.test_output.get_mesh().coordinates()
    xgold = exx * x + (2*exz) * z
    ygold = eyy * y
    
    test = self.get_test_field('Displacement',cycle=0)
    
    error = max(abs(test[:,0]-xgold))
    if error > tol:
      fail += 1
      print 'x-displacement: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'x-displacement: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)
    
    error = max(abs(test[:,1]-ygold))
    if error > tol:
      fail += 1
      print 'y-displacement: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'y-displacement: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)
    
    error = max(abs(test[:,2]))
    if error > tol:
      fail += 1
      print 'z-displacement: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'z-displacement: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    self.assertTrue(fail == 0)


if __name__ == '__main__':
  import unittest
  unittest.main()

