#!/usr/bin/env python

import sys
import os

import unittest
import numpy

import Truchas
import TruchasTest

class PorousDrag(TruchasTest.GoldenTestCase):

  test_name = 'porous-drag'
  num_procs = 4 # with a parallel executable

  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())

  def get_test_field(self,field,cycle,serialize=True,region=None):
    return self.test_output.get_simulation().find_series(cycle=cycle).get_data(field,serialize,region=region)

  def test_velocity(self):
    '''Verify velocity field'''

    n = 114
    fail = 0

    test = self.get_test_field('Z_VC',cycle=n)

    vx = 0.5
    tol = 5.0e-8
    error = max(abs((test[:,0]-vx)/vx))
    if (error > tol):
      fail += 1
      print 'x-velocity: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'x-velocity: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    tol = 5.0e-9
    for k in [1,2]:
      error = max(abs(test[:,k]))
      if (error > tol):
        fail += 1
        print '%s-velocity: max abs error = %8.2e: FAIL (tol=%8.2e)'%(('x','y','z')[k],error,tol)
      else:
        print '%s-velocity: max abs error = %8.2e: PASS (tol=%8.2e)'%(('x','y','z')[k],error,tol)

    self.assertTrue(fail == 0)


  def test_pressure(self):
    '''Verify pressure field'''

    n = 114
    tol = 5.0e-6

    # The centroids function doesn't serialize, so we don't want to here either.
    test = self.get_test_field('Z_P',cycle=n,serialize=False)
    xc = self.test_output.get_mesh().centroids()

    # Analytic linear pressure field
    length = 20
    Pin = 20000.0
    gold = Pin*(1.0-xc[:,0]/length)

    error = max(abs(test - gold))
    if (error > tol):
      print 'pressure: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'pressure: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)


if __name__ == '__main__':
  import unittest
  unittest.main()
