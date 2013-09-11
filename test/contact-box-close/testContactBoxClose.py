#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class ContactBoxClose(TruchasTest.GoldenTestCase):

  test_name = 'contact-box-close'
  num_procs = 4 # with a parallel executable

  gap_block_id = 3
  other_block_ids = [1, 2]

  true_cells = None
  gap_cells = None

  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())
    self.gold_output=Truchas.TruchasOutput(self.get_golden_output_file())
    self.test_sim = self.test_output.get_simulation()
    self.gap_cells  = Truchas.TruchasRegion(self.test_sim,self.gap_block_id)

    # TODO: get_data with region is currently broken for vector field data, so we go
    # ahead and do comparisons on all cells, including gap cells, instead of using
    # the true_cells region.  This only works because Truchas supplies appropriate
    # dummy values on gap cells.  Eventually we want to test only on real cells and
    # not rely on the current Truchas behavior.
    # self.true_cells = Truchas.TruchasRegion(self.test_sim,self.other_block_ids)

  def get_test_field(self,field,cycle,region=None):
    return self.test_output.get_simulation().find_series(cycle=cycle).get_data(field,region)

  def get_gold_field(self,field,cycle,region=None):
    return self.gold_output.get_simulation().find_series(cycle=cycle).get_data(field,region)

  # TODO: The SKIP_* tests below need to be enabled (delete the SKIP_) once
  # regions can be used with node-based fields.

  def test_open_stress(self):
    '''Verify stresses just before gap closing'''
    n = 15
    tol = 1.0e3 # a little sloppy
    test = self.get_test_field('sigma',cycle=n,region=self.true_cells)
    gold = self.get_gold_field('sigma',cycle=n,region=self.true_cells)
    fail = 0
    for j in range(6):
      error = max(abs(test[:,j]-gold[:,j]))
      if error > tol:
        fail += 1
        print 'Cycle %2d, sigma%1d: max abs error = %8.2e: FAIL (tol=%8.2e)'%(n,j+1,error,tol)
      else:
        print 'Cycle %2d, sigma%1d: max abs error = %8.2e: PASS (tol=%8.2e)'%(n,j+1,error,tol)
    self.assertTrue(fail == 0)

  def test_open_strain(self):
    '''Verify strains just before gap closing'''
    n = 15
    tol = 2.0e-8
    test = self.get_test_field('epsilon',cycle=n,region=self.true_cells)
    gold = self.get_gold_field('epsilon',cycle=n,region=self.true_cells)
    fail = 0
    for j in range(6):
      error = max(abs(test[:,j]-gold[:,j]))
      if error > tol:
        fail += 1
        print 'Cycle %2d, epsilon%1d: max abs error = %8.2e: FAIL (tol=%8.2e)'%(n,j+1,error,tol)
      else:
        print 'Cycle %2d, epsilon%1d: max abs error = %8.2e: PASS (tol=%8.2e)'%(n,j+1,error,tol)
    self.assertTrue(fail == 0)

  def SKIP_test_open_traction(self):
    '''Verify normal traction just before gap closing'''
    n = 15
    tol = 1.0e-4
    gold = 0.0
    test = self.get_test_field('NTRAC_02',cycle=n,region=self.gap_cells)
    error = max(abs(test-gold))
    if error > tol:
      print 'Cycle %2d: max abs error = %8.2e: FAIL (tol=%8.2e)'%(n,error,tol)
      self.assertTrue(False)
    else:
      print 'Cycle %2d: max abs error = %8.2e: PASS (tol=%8.2e)'%(n,error,tol)

  def test_close_stress(self):
    '''Verify stresses just after gap closing'''
    n = 18
    tol = 1.0e4 # sloppy
    test = self.get_test_field('sigma',cycle=n,region=self.true_cells)
    gold = self.get_gold_field('sigma',cycle=n,region=self.true_cells)
    fail = 0
    for j in range(6):
      error = max(abs(test[:,j]-gold[:,j]))
      if error > tol:
        fail += 1
        print 'Cycle %2d, sigma%1d: max abs error = %8.2e: FAIL (tol=%8.2e)'%(n,j+1,error,tol)
      else:
        print 'Cycle %2d, sigma%1d: max abs error = %8.2e: PASS (tol=%8.2e)'%(n,j+1,error,tol)
    self.assertTrue(fail == 0)

  def test_close_strain(self):
    '''Verify strains just after gap closing'''
    n = 18
    tol = 2.0e-8
    test = self.get_test_field('epsilon',cycle=n,region=self.true_cells)
    gold = self.get_gold_field('epsilon',cycle=n,region=self.true_cells)
    fail = 0
    for j in range(6):
      error = max(abs(test[:,j]-gold[:,j]))
      if error > tol:
        fail += 1
        print 'Cycle %2d, epsilon%1d: max abs error = %8.2e: FAIL (tol=%8.2e)'%(n,j+1,error,tol)
      else:
        print 'Cycle %2d, epsilon%1d: max abs error = %8.2e: PASS (tol=%8.2e)'%(n,j+1,error,tol)
    self.assertTrue(fail == 0)

  def test_closed_displacement(self):
    '''Verify displacement magnitude just after gap closing'''
    n = 18
    abs_tol = 1.0e-9
    rel_tol = 1.0e-4
    test = self.get_test_field('Displacement',cycle=n)
    gold = self.get_gold_field('Displacement',cycle=n)
    dtest = (test[:,0]**2 + test[:,1]**2 + test[:,2]**2)**0.5
    dgold = (gold[:,0]**2 + gold[:,1]**2 + gold[:,2]**2)**0.5
    error = max(abs(dtest-dgold)/(abs_tol + rel_tol*abs(dgold)))
    if error > 1.0:
      print 'Cycle %2d: displacement error = %8.2e: FAIL (tol=1.0)'%(n,error)
      self.assertTrue(False)
    else:
      print 'Cycle %2d: displacement error = %8.2e: PASS (tol=1.0)'%(n,error)

  def test_closed_temperature(self):
    '''Verify temperature just after gap closing'''
    n = 18
    tol = 1.0e-5
    test = self.get_test_field('Z_TEMP',cycle=n,region=self.true_cells)
    gold = self.get_gold_field('Z_TEMP',cycle=n,region=self.true_cells)
    error = max(abs((test-gold)/gold))
    if error > tol:
      print 'Cycle %2d: max rel error = %8.2e: FAIL (tol=%8.2e)'%(n,error,tol)
      self.assertTrue(False)
    else:
      print 'Cycle %2d: max rel error = %8.2e: PASS (tol=%8.2e)'%(n,error,tol)

  def test_initial_stress(self):
    '''Verify initial stresses'''
    n = 0
    tol = 1.0e2
    gold = [7.62665e6, 1.90666e6, 5.71999e6, 3.81333e6, 6.60487e6, 3.30244e6]
    test = self.get_test_field('sigma',cycle=n,region=self.true_cells)
    fail = 0
    for j in range(6):
      error = max(abs(test[:,j]-gold[j]))
      if error > tol:
        fail += 1
        print 'Cycle %2d, sigma%1d: max abs error = %8.2e: FAIL (tol=%8.2e)'%(n,j+1,error,tol)
      else:
        print 'Cycle %2d, sigma%1d: max abs error = %8.2e: PASS (tol=%8.2e)'%(n,j+1,error,tol)
    self.assertTrue(fail == 0)

  def test_initial_strain(self):
    '''Verify initial strains'''
    n = 0
    tol = 1.0e-8
    gold = [-1.46667e-4, -2.56666e-4, -1.83333e-4, 7.33332e-5, 1.27017e-4, 6.35085e-5]
    test = self.get_test_field('epsilon',cycle=n,region=self.true_cells)
    fail = 0
    for j in range(6):
      error = max(abs(test[:,j]-gold[j]))
      if error > tol:
        fail += 1
        print 'Cycle %2d, epsilon%1d: max abs error = %8.2e: FAIL (tol=%8.2e)'%(n,j+1,error,tol)
      else:
        print 'Cycle %2d, epsilon%1d: max abs error = %8.2e: PASS (tol=%8.2e)'%(n,j+1,error,tol)
    self.assertTrue(fail == 0)

  def SKIP_test_initial_traction(self):
    '''Verify initial normal traction'''
    n = 0
    tol = 1.0e-4
    gold = 0.0
    test = self.get_test_field('NTRAC_02',cycle=n,region=self.gap_cells)
    error = max(abs(test-gold))
    if error > tol:
      print 'Cycle %2d: max abs error = %8.2e: FAIL (tol=%8.2e)'%(n,error,tol)
      self.assertTrue(False)
    else:
      print 'Cycle %2d: max abs error = %8.2e: PASS (tol=%8.2e)'%(n,error,tol)

  def test_final_stress(self):
    '''Verify final stresses'''
    n = 42
    tol = 1.0
    gold = [-2.288e7, -5.720e6, -1.716e7, 0.0, 0.0, -9.90733e6]
    test = self.get_test_field('sigma',cycle=n,region=self.true_cells)
    fail = 0
    for j in range(6):
      error = max(abs(test[:,j]-gold[j]))
      if error > tol:
        fail += 1
        print 'Cycle %2d, sigma%1d: max abs error = %8.2e: FAIL (tol=%8.2e)'%(n,j+1,error,tol)
      else:
        print 'Cycle %2d, sigma%1d: max abs error = %8.2e: PASS (tol=%8.2e)'%(n,j+1,error,tol)
    self.assertTrue(fail == 0)

  def test_final_strain(self):
    '''Verify final strains'''
    n = 42
    tol = 1.0e-9
    gold = [0.0, 3.3e-4, 1.1e-4, 0.0, 0.0, -1.90526e-4]
    test = self.get_test_field('epsilon',cycle=n,region=self.true_cells)
    fail = 0
    for j in range(6):
      error = max(abs(test[:,j]-gold[j]))
      if error > tol:
        fail += 1
        print 'Cycle %2d, epsilon%1d: max abs error = %8.2e: FAIL (tol=%8.2e)'%(n,j+1,error,tol)
      else:
        print 'Cycle %2d, epsilon%1d: max abs error = %8.2e: PASS (tol=%8.2e)'%(n,j+1,error,tol)
    self.assertTrue(fail == 0)

  def SKIP_test_final_traction(self):
    '''Verify final normal traction'''
    n = 42
    tol = 1.0e-6
    gold = -2.288e7
    test = self.get_test_field('NTRAC_02',cycle=n,region=self.gap_cells)
    error = max(abs(test-gold))
    if error > tol:
      print 'Cycle %2d: max abs error = %8.2e: FAIL (tol=%8.2e)'%(n,error,tol)
      self.assertTrue(False)
    else:
      print 'Cycle %2d: max abs error = %8.2e: PASS (tol=%8.2e)'%(n,error,tol)

  def test_final_temperature(self):
    '''Verify final temperature'''
    n = 42
    tol = 1.0e-8
    test = self.get_test_field('Z_TEMP',cycle=n,region=self.true_cells)
    gold = 308.0
    error = max(abs((test-gold)/gold))
    if error > tol:
      print 'Cycle %2d: max rel error = %8.2e: FAIL (tol=%8.2e)'%(n,error,tol)
      self.assertTrue(False)
    else:
      print 'Cycle %2d: max rel error = %8.2e: PASS (tol=%8.2e)'%(n,error,tol)


if __name__ == '__main__':
  import unittest
  unittest.main()

