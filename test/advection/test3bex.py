#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class mytest(TruchasTest.GoldenTestCase):

  test_name = 'advection-3b-exact'
  num_procs = 4 # with a parallel executable

  def test_final_vof(self):
    '''verifying volume fractions at final time'''
    tol_max = 0.5; tol_l2 = 0.05

    time = self.test_output.get_simulation().find_series(id=2).time
    test = self.test_output.get_simulation().find_series(id=2).get_data('VOF')[:,0]

    # From a special run to compute vof for final exact configuration
    gold = self.gold_output.get_simulation().find_series(id=1).get_data('VOF')[:,0]

    error_max = max(abs(test-gold))
    error_l2 = numpy.linalg.norm(test-gold) / numpy.sqrt(len(test))

    # max <= l2 <= sqrt(ncell)*max = sqrt(3)*31*max
    passfail = 'PASS' if error_max <= tol_max else 'FAIL'
    print 'vof at t=%8.2e: max error = %8.2e: %s (tol=%8.2e)' % (time, error_max, passfail, tol_max)

    passfail = 'PASS' if error_l2 <= tol_l2 else 'FAIL'
    print 'vof at t=%8.2e: l2 error = %8.2e: %s (tol=%8.2e)' % (time, error_l2, passfail, tol_l2)

    self.assertTrue(error_max <= tol_max)
    self.assertTrue(error_l2 <= tol_l2)

if __name__ == '__main__':
  import unittest
  unittest.main()
