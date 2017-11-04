#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class VFRad1Tet(TruchasTest.GoldenTestCase):

  test_name = 'vfrad1-tet'
  num_procs = 4 # with a parallel executable
  
  def test_final_cycle_number(self):
    '''VFRAD1Tet: checking the final cycle number'''
    test_series = self.test_output.get_simulation().find_series(id=2)
    gold_series = self.gold_output.get_simulation().find_series(id=2)
    self.assertTrue(test_series.cycle == gold_series.cycle)

  def test_final_temperature(self):
    '''VFRAD1Tet: verifying the final temperature field'''
    tol = 1.0e-5
    test = self.test_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    gold = self.gold_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    error = max(abs(test-gold)/gold)
    print 'final temp rel error=', error, '(tol=', tol, ')'
    self.assertTrue(error <= tol)

if __name__ == '__main__':
  import unittest
  unittest.main()

