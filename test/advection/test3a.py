#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class mytest(TruchasTest.GoldenTestCase):

  test_name = 'advection-3a'
  num_procs = 4 # with a parallel executable

  def test_final_vof(self):
    '''verifying volume fractions at final time'''
    tol = 5e-9

    time = self.test_output.get_simulation().find_series(id=2).time
    test = self.test_output.get_simulation().find_series(id=2).get_data('VOF')[:,0]
    gold = self.gold_output.get_simulation().find_series(id=2).get_data('VOF')[:,0]

    error = max(abs(test-gold))
    passfail = 'PASS' if error <= tol else 'FAIL'
    print 'vof at t=%8.2e: max error = %8.2e: %s (tol=%8.2e)' % (time, error, passfail, tol)

    self.assertTrue(error <= tol)

if __name__ == '__main__':
  import unittest
  unittest.main()
