#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class DS11(TruchasTest.GoldenTestCase):

  test_name='ds11'
  num_procs=4

  gap_block_ids = [10]
  other_block_ids = [1, 2]

  true_cells = None
  gap_cells = None

  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())
    self.gold_output=Truchas.TruchasOutput(self.get_golden_output_file())
    self.test_sim = self.test_output.get_simulation()
    self.gold_sim = self.gold_output.get_simulation()
    self.test_reg = Truchas.TruchasRegion(self.test_sim,self.other_block_ids)
    self.gold_reg = Truchas.TruchasRegion(self.gold_sim,self.other_block_ids)
  
  def test_cycle_numbers(self):
    '''DS8: checking the cycle numbers'''
    for n in [2,4]:
      test_series = self.test_output.get_simulation().find_series(id=n)
      gold_series = self.gold_output.get_simulation().find_series(id=n)
      self.assertTrue(test_series.cycle == gold_series.cycle)

  def test_temperature(self):
    '''Verifying the temperature field'''
    tol = 1.0e-9
    fail = 0
    for n in [3,4]:
      gold=self.gold_output.get_simulation().find_series(id=n).get_data('Z_TEMP',region=self.gold_reg)
      test=self.test_output.get_simulation().find_series(id=n).get_data('Z_TEMP',region=self.test_reg)
      error=max(abs((test-gold)/gold))
      if (error > tol):
        fail += 1
        print 'Output %3d: max temp rel error = %8.2e: FAIL (tol=%8.2e)'%(n,error,tol)
      else:
        print 'Output %3d: max temp rel error = %8.2e: PASS (tol=%8.2e)'%(n,error,tol)
    self.assertTrue(fail == 0)


  def test_fluid_fraction(self):
    '''Verifying the fluid volume fraction field'''
    tol = 1.0e-8
    fail = 0
    for n in [3,4]:
      gold=self.gold_output.get_simulation().find_series(id=n).get_data('VOF')
      test=self.test_output.get_simulation().find_series(id=n).get_data('VOF')
      error=max(abs((test[:,2]-gold[:,2])))  # comp 2 is fluid
      if (error > tol):
        fail += 1
        print 'Output %3d: max vof error = %8.2e: FAIL (tol=%8.2e)'%(n,error,tol)
      else:
        print 'Output %3d: max vof error = %8.2e: PASS (tol=%8.2e)'%(n,error,tol)
    self.assertTrue(fail == 0)

if __name__ == '__main__':
  import unittest
  unittest.main()


