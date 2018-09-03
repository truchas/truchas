#!/usr/bin/env python

import sys
import os

import numpy
import math

import Truchas
import TruchasTest

class mytest(TruchasTest.GoldenTestCase):

  test_name = 'pipe-flow-old-3'
  num_procs = 4 # with a parallel executable

  # Override the default setUp, omitting the opening of the golden output
  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass() # This runs Truchas
    self.test_output = Truchas.TruchasOutput(self.get_output_file())

  def test_final_velocity(self):
    '''PIPE-FLOW-3: final velocity'''

    centroids = self.test_output.get_mesh().centroids()
    u = (1 - centroids[:,1]**2)/2/math.sqrt(2)
    w = u

    data = self.test_output.get_simulation().find_series(id=2).get_data('Z_VC',serialize=False)

    fail = 0

    tol = 1.5e-3
    error = max(abs(data[:,0] - u))
    if (error > tol):
      fail += 1
      print 'x-velocity: max error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'x-velocity: max error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    tol = 2.0e-11
    error = max(abs(data[:,1]))
    if (error > tol):
      fail += 1
      print 'y-velocity: max error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'y-velocity: max error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    tol = 1.5e-3
    error = max(abs(data[:,2] - w))
    if (error > tol):
      fail += 1
      print 'z-velocity: max error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'z-velocity: max error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    self.assertTrue(fail == 0)

if __name__ == '__main__':
  import unittest
  unittest.main()

