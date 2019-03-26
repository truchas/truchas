#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class HTVoid1(TruchasTest.GoldenTestCase):

  test_name = 'htvoid1'
  num_procs = 4 # with a parallel executable

  def test_final_cycle_number(self):
    '''HTVOID1: checking the final cycle number'''
    test_series = self.test_output.get_simulation().find_series(id=2)
    gold_series = self.gold_output.get_simulation().find_series(id=2)
    self.assertTrue(test_series.cycle == gold_series.cycle)

  def test_final_temp(self):
    '''HTVOID1: verifying the temperature field at final time'''
    tol = 1.0e-9
    test = self.test_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    gold = self.gold_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    error = max(abs(test-gold))
    print 'final temp max error=', error, '(tol=', tol, ')'
    self.assertTrue(error <= tol)

if __name__ == '__main__':
  import unittest
  unittest.main()
