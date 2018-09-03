#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class mytest(TruchasTest.GoldenTestCase):

  test_name = 'pipe-flow-inviscid-1'
  num_procs = 4 # with a parallel executable

  # Override the default setUp, omitting the opening of the golden output
  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass() # This runs Truchas
    self.test_output = Truchas.TruchasOutput(self.get_output_file())

  def test_early_velocity(self):
    '''PIPE-FLOW-INVISCID-1: early velocity'''

    time = self.test_output.get_simulation().find_series(id=2).time
    data = self.test_output.get_simulation().find_series(id=2).get_data('Z_VC')
    print 'velocity at t=%11.5e'%(time)

    fail = 0

    tol = 1.0e-11
    error = max(abs((data[:,0] - 2*time)/(2*time)))
    if (error > tol):
      fail += 1
      print 'x-velocity: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'x-velocity: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    tol = 1.0e-12
    error = max(abs(data[:,1]))
    if (error > tol):
      fail += 1
      print 'y-velocity: max error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'y-velocity: max error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs(data[:,2]))
    if (error > tol):
      fail += 1
      print 'z-velocity: max error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'z-velocity: max error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    self.assertTrue(fail == 0)

  def test_final_velocity(self):
    '''PIPE-FLOW-INVISCID-1: final velocity'''

    time = self.test_output.get_simulation().find_series(id=3).time
    data = self.test_output.get_simulation().find_series(id=3).get_data('Z_VC')
    print 'velocity at t=%11.5e'%(time)

    fail = 0

    tol = 1.0e-12
    error = max(abs(data[:,0] - 2*time))
    if (error > tol):
      fail += 1
      print 'x-velocity: max error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'x-velocity: max error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    tol = 1.0e-13
    error = max(abs(data[:,1]))
    if (error > tol):
      fail += 1
      print 'y-velocity: max error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'y-velocity: max error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs(data[:,2]))
    if (error > tol):
      fail += 1
      print 'z-velocity: max error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'z-velocity: max error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    self.assertTrue(fail == 0)

if __name__ == '__main__':
  import unittest
  unittest.main()

