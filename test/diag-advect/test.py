#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class DiagAdvect(TruchasTest.GoldenTestCase):

  test_name = 'diag-advect'
  num_procs = 4 # with a parallel executable

  def test_final_cycle_number(self):
    '''DIAG-ADVECT: checking the final cycle number'''
    test_series = self.test_output.get_simulation().find_series(id=6)
    gold_series = self.gold_output.get_simulation().find_series(id=6)
    self.assertTrue(test_series.cycle == gold_series.cycle)
  
  def test_final_fluid_frac(self):
    '''DIAG-ADVECT: verifying volume fractions at final time'''
    tol = 1.0e-6
    test = self.test_output.get_simulation().find_series(id=6).get_data('VOF')
    gold = self.gold_output.get_simulation().find_series(id=6).get_data('VOF')
    error = max(abs(test[:,0]-gold[:,0])) # comp 0 is circle
    print 'final vof max error=', error, '(tol=', tol, ')'
    self.assertTrue(error <= tol)

if __name__ == '__main__':
  import unittest
  unittest.main()

