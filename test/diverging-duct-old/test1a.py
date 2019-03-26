#!/usr/bin/env python

import sys
import os

import unittest
import numpy
import math

import Truchas
import TruchasTest

class DivergingDuct(TruchasTest.GoldenTestCase):

  test_name = 'diverging-duct-old-1a'
  num_procs = 4 # with a parallel executable

  # Override the default setUp, omitting the opening of the golden output
  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())

  def test_pressure(self):
    '''Verify pressure field against expected analytic value'''

    # The centroids function does not serialize, so we don't want to here either.
    test = self.test_output.get_simulation().find_series(id=2).get_data('Z_P',serialize=False)

    # Analytic pressure solution at cell centrioids
    cc = self.test_output.get_mesh().centroids()
    p = 2.5 - 2 / (cc[:,0]**2 + cc[:,1]**2)

    tol = 0.16
    error = max(abs((test-p)/p))
    if error > tol:
      print 'pressure: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'pressure: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

  def test_velocity(self):
    '''Verify velocity field against expected analytic value'''

    # The centroids function does not serialize, so we don't want to here either.
    test = self.test_output.get_simulation().find_series(id=2).get_data('Z_VC',serialize=False)

    # Analytic velocity solution at cell centroids
    cc = self.test_output.get_mesh().centroids()
    u = cc[:,0] / (cc[:,0]**2 + cc[:,1]**2)
    v = cc[:,1] / (cc[:,0]**2 + cc[:,1]**2)

    fail = 0

    tol = 0.015
    error = max(abs(test[:,0] - u)/max(abs(u)))
    if error > tol:
      fail += 1
      print 'x-velocity: max error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'x-velocity: max error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    tol = 0.18
    error = max(abs(test[:,1] - v)/max(abs(v)))
    if error > tol:
      fail += 1
      print 'y-velocity: max error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'y-velocity: max error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    tol = 1.0e-11
    error = max(abs(test[:,2]))
    if error > tol:
      fail += 1
      print 'z-velocity: max error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'z-velocity: max error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    self.assertTrue(fail == 0)

if __name__ == '__main__':
  import unittest
  unittest.main()
