#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class DS3(TruchasTest.GoldenTestCase):

  test_name = 'ds3'
  num_procs = 4 # with a parallel executable

  def test_final_cycle_number(self):
    '''DS3: checking the final cycle number'''
    test_series = self.test_output.get_simulation().find_series(id=2)
    gold_series = self.gold_output.get_simulation().find_series(id=2)
    self.assertTrue(test_series.cycle == gold_series.cycle)

  def test_final_temperature(self):
    '''DS3: verifying the final temperature field'''
    tol = 1.0e-10
    T    = self.test_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    Tref = self.gold_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    error = max(abs(T-Tref)/Tref)
    print 'final temp max rel error=', error, '(tol=', tol, ')'
    self.assertTrue(error <= tol)

  def test_final_concentration(self):
    '''DS3: verifying the final concentration field'''
    tol = 1.0e-10
    C    = self.test_output.get_simulation().find_series(id=2).get_data('phi1')
    Cref = self.gold_output.get_simulation().find_series(id=2).get_data('phi1')
    error = max(abs(C-Cref)/Cref)
    print 'final conc max rel error=', error, '(tol=', tol, ')'
    self.assertTrue(error <= tol)

if __name__ == '__main__':
  import unittest
  unittest.main()

