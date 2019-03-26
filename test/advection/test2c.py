#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class mytest(TruchasTest.GoldenTestCase):

  test_name = 'advection-2c'
  num_procs = 4 # with a parallel executable

  def test_final_vof(self):
    '''verifying volume fractions at final time'''

    time = self.test_output.get_simulation().find_series(id=2).time
    test = self.test_output.get_simulation().find_series(id=2).get_data('VOF')

    # From a special run to compute vof for final exact configuration
    gold = self.gold_output.get_simulation().find_series(id=1).get_data('VOF')

    fail = 0
    max_tol = 0.3
    max_error = numpy.amax(abs(test-gold))
    if max_error > max_tol:
      fail += 1
      print 'vof at t=%8.2e: max error = %8.2e: FAIL (tol=%8.2e)'%(time,max_error,max_tol)
    else:
      print 'vof at t=%8.2e: max error = %8.2e: PASS (tol=%8.2e)'%(time,max_error,max_tol)

    # max <= l2 <= sqrt(ncell)*max = sqrt(3)*31*max
    l2_tol = 4*max_tol
    l2_error = numpy.linalg.norm(test-gold)
    if l2_error > l2_tol:
      fail += 1
      print 'vof at t=%8.2e: l2 error = %8.2e: FAIL (tol=%8.2e)'%(time,l2_error,l2_tol)
    else:
      print 'vof at t=%8.2e: l2 error = %8.2e: PASS (tol=%8.2e)'%(time,l2_error,l2_tol)

    self.assertTrue(fail==0)

if __name__ == '__main__':
  import unittest
  unittest.main()

