#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class SloshingFlow(TruchasTest.GoldenTestCase):

  test_name = 'sloshing-flow'
  num_procs = 4 # with a parallel executable

  def test_final_cycle_number(self):
    '''SLOSHING FLOW: checking the final cycle number'''
    test_series = self.test_output.get_simulation().find_series(id=2)
    gold_series = self.gold_output.get_simulation().find_series(id=2)
    self.assertTrue(test_series.cycle == gold_series.cycle)
  
  def test_final_temp(self):
    '''SLOSHING FLOW: verifying the temperature field at final time'''
    tol = 1.0e-6
    test = self.test_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    gold = self.gold_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    error = max(abs(test-gold))
    print 'final temp max error=', error, '(tol=', tol, ')'
    self.assertTrue(error <= tol)

  def test_final_fluid_frac(self):
    '''SLOSHING FLOW: verifying the solid volume fraction at final time'''
    tol = 1.0e-6
    test = self.test_output.get_simulation().find_series(id=2).get_data('VOF')
    gold = self.gold_output.get_simulation().find_series(id=2).get_data('VOF')
    error = max(abs(test[:,0]-gold[:,0])) # comp 0 is fluid
    print 'final vof max error=', error, '(tol=', tol, ')'
    self.assertTrue(error <= tol)

if __name__ == '__main__':
  import unittest
  unittest.main()

