#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class SimpleGap(TruchasTest.GoldenTestCase):

  test_name = 'simple-gap'
  num_procs = 4 # with a parallel executable

  gap_block_ids = [3, 4]
  other_block_ids = [1, 2]

  true_cells = None
  gap_cells = None

  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())
    self.gold_output=Truchas.TruchasOutput(self.get_golden_output_file())
    self.test_sim = self.test_output.get_simulation()
    self.gap_cells  = Truchas.TruchasRegion(self.test_sim,self.gap_block_ids)
    self.true_cells = Truchas.TruchasRegion(self.test_sim,self.other_block_ids)

  def get_test_field(self,field,cycle,region=None):
    return self.test_output.get_simulation().find_series(cycle=cycle).get_data(field,region=region)

  def get_gold_field(self,field,cycle,region=None):
    return self.gold_output.get_simulation().find_series(cycle=cycle).get_data(field,region=region)

  def test_stress(self):
    '''Verify initial stress'''
    tol = 1.0e-6
    test = self.get_test_field('sigma',cycle=0,region=self.true_cells)
    gold = self.get_gold_field('sigma',cycle=0,region=self.true_cells)
    fail = 0
    for j in range(6):
      error = max(abs((test[:,j]-gold[:,j])/gold[:,j]))
      if error > tol:
        fail += 1
        print 'sigma%1d: max rel error = %8.2e: FAIL (tol=%8.2e)'%(j+1,error,tol)
      else:
        print 'sigma%1d: max rel error = %8.2e: PASS (tol=%8.2e)'%(j+1,error,tol)
    self.assertTrue(fail == 0)

  def test_strain(self):
    '''Verify initial strain'''
    tol = 1.0e-7
    test = self.get_test_field('epsilon',cycle=0,region=self.true_cells)
    gold = self.get_gold_field('epsilon',cycle=0,region=self.true_cells)
    fail = 0
    for j in range(6):
      error = max(abs((test[:,j]-gold[:,j])/gold[:,j]))
      if error > tol:
        fail += 1
        print 'epsilon%1d: max rel error = %8.2e: FAIL (tol=%8.2e)'%(j+1,error,tol)
      else:
        print 'epsilon%1d: max rel error = %8.2e: PASS (tol=%8.2e)'%(j+1,error,tol)
    self.assertTrue(fail == 0)

  def test_traction_7(self):
    '''Verify normal traction on sideset 7'''
    tol = 1.0e-6
    gold = self.get_gold_field('NTRAC_07',cycle=0)
    gold = numpy.ma.masked_equal(gold,0.0)
    test = self.get_test_field('NTRAC_07',cycle=0)
    error = numpy.ma.max(abs((gold-test)/gold))
    if error > tol:
      print 'normal traction 7: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'normal traction 7: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

  def test_traction_8(self):
    '''Verify normal traction on sideset 8'''
    tol = 1.0e-6
    gold = self.get_gold_field('NTRAC_08',cycle=0)
    gold = numpy.ma.masked_equal(gold,0.0)
    test = self.get_test_field('NTRAC_08',cycle=0)
    error = numpy.ma.max(abs((gold-test)/gold))
    if error > tol:
      print 'normal traction 7: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'normal traction 7: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)


if __name__ == '__main__':
  import unittest
  unittest.main()

