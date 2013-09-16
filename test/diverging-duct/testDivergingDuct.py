#!/usr/bin/env python

import sys
import os

import unittest
import numpy

import Truchas
import TruchasTest

class DivergingDuct(TruchasTest.GoldenTestCase):

  test_name = 'diverging-duct'
  num_procs = 4 # with a parallel executable

  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())

  def get_test_field(self,field,cycle,serialize=True,region=None):
    return self.test_output.get_simulation().find_series(cycle=cycle).get_data(field,serialize,region=region)

  def test_pressure(self):
    '''Verify pressure field against expected analytic value'''
    
    n = 274
    tol = 5.0e-3
    # The centroids function does not serialize, so we don't want to here either.
    test = self.get_test_field('Z_P',cycle=n,serialize=False)
    
    xc = self.test_output.get_mesh().centroids()
    gold = 1.0 - 0.5/(1.0 + 0.025*xc[:,0])**2
    print max(gold), min(gold)
    error = max(abs((test-gold)/gold))
    if error > tol:
      print 'pressure: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'pressure: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

  def test_velocity(self):
    '''Verify velocity field against expected analytic value'''
    
    n = 274
    tol = 5.0e-3
    # The centroids function does not serialize, so we don't want to here either.
    test = self.get_test_field('Z_VC',cycle=n,serialize=False)
    
    xc = self.test_output.get_mesh().centroids()
    gold = 1.0/(1.0 + 0.025*xc[:,0])
    print max(gold), min(gold)
    error = max(abs((test[:,0]-gold))/gold)
    if error > tol:
      print 'velocity: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'velocity: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)


if __name__ == '__main__':
  import unittest
  unittest.main()
