#!/usr/bin/env python

import sys
import os

import unittest
import numpy
import math

import Truchas
import TruchasTest

class mytest(TruchasTest.GoldenTestCase):

  test_name = 'hydrostatic-old-1b'
  num_procs = 4 # with a parallel executable

  # Override the default setUp, omitting the opening of the golden output
  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())

  def pressure_test(self, id):

    # The centroids function does not serialize, so we don't want to here either.
    test = self.test_output.get_simulation().find_series(id).get_data('Z_P',serialize=False)
    time = self.test_output.get_simulation().find_series(id).time

    # Analytic pressure solution at cell centrioids
    cc = self.test_output.get_mesh().centroids()
    p = -math.sqrt(2)*(cc[:,1] + cc[:,2])
    
    # Error array, accounting for arbitrary constant
    d = (test-p) - numpy.mean(test-p)

    tol = 4.0e-9
    error = max(abs(d))
    if error > tol:
      print 'pressure at t=%8.2e: max error = %8.2e: FAIL (tol=%8.2e)'%(time,error,tol)
      self.assertTrue(False)
    else:
      print 'pressure at t=%8.2e: max error = %8.2e: PASS (tol=%8.2e)'%(time,error,tol)

  def test_pressure1(self):
    '''Verify initial pressure field'''
    self.pressure_test(1)

  def test_pressure2(self):
    '''Verify final pressure field'''
    self.pressure_test(2)

  def test_final_velocity(self):
    '''Verify final velocity field'''

    test = self.test_output.get_simulation().find_series(id=2).get_data('Z_VC',serialize=False)

    tol = 1.0e-13
    error = max(numpy.sqrt(test[:,0]**2 + test[:,1]**2 + test[:,2]**2))
    if error > tol:
      print 'max velocity = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'max velocity = %8.2e: PASS (tol=%8.2e)'%(error,tol)

if __name__ == '__main__':
  import unittest
  unittest.main()
