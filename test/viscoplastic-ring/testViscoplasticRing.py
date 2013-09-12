#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class ViscoplasticRing(TruchasTest.GoldenTestCase):

  test_name = 'viscoplastic-ring'
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
    
    # TODO: Original test tested NTRAC_04 over the entire mesh -- we probably
    # only want to do this over gap elements, right?
    #self.gap_cells  = Truchas.TruchasRegion(self.test_sim,self.gap_block_id)
    self.true_cells = Truchas.TruchasRegion(self.test_sim,self.other_block_ids)

  def get_test_field(self,field,cycle,region=None):
    return self.test_output.get_simulation().find_series(cycle=cycle).get_data(field,region)

  def get_gold_field(self,field,cycle,region=None):
    return self.gold_output.get_simulation().find_series(cycle=cycle).get_data(field,region)

  # TODO: The SKIP_* tests below need to be enabled (delete the SKIP_) once
  # regions can be used with node-based fields.

  def test_initial_stress(self):
    '''Verify initial stress'''
    n = 0
    tol = 0.1
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

  def test_initial_strain(self):
    '''Verify initial strain'''
    n = 0
    tol = 1.0e-10
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

  def test_initial_traction(self):
    '''Verify initial normal traction'''
    n = 0
    tol = 0.1
    test = self.get_test_field('NTRAC_04',cycle=n,region=self.gap_cells)
    gold = self.get_gold_field('NTRAC_04',cycle=n,region=self.gap_cells)
    error = max(abs(test-gold))
    if error > tol:
      print 'Cycle %2d: max abs error = %8.2e: FAIL (tol=%8.2e)'%(n,error,tol)
      self.assertTrue(False)
    else:
      print 'Cycle %2d: max abs error = %8.2e: PASS (tol=%8.2e)'%(n,error,tol)

  def test_final_stress(self):
    '''Verify final stress'''
    n = 5
    tol = 5.0e3 # loose
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

  def test_final_strain(self):
    '''Verify final strain'''
    n = 5
    tol = 5.0e-9
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
    
  def eps_eff(self,eps):
    return ((2.0/9.0)*((eps[:,0]-eps[:,1])**2 + (eps[:,1]-eps[:,2])**2 + (eps[:,2]-eps[:,0])**2) +
            (4.0/3.0)*(eps[:,3]**2 + eps[:,4]**2 + eps[:,5]**2))**0.5

  def test_final_plastic_strain(self):
    '''Verify final plastic strain'''
    n = 5
    tol = 5.0e-9
    test = self.eps_eff(self.get_test_field('e_plastic',cycle=n,region=self.true_cells))
    gold = self.eps_eff(self.get_gold_field('e_plastic',cycle=n,region=self.true_cells))
    error = max(abs(test-gold))
    if error > tol:
      print 'Cycle %2d: eps_plas: max abs error = %8.2e: FAIL (tol=%8.2e)'%(n,error,tol)
      self.assertTrue(False)
    else:
      print 'Cycle %2d: eps_plas: max abs error = %8.2e: PASS (tol=%8.2e)'%(n,error,tol)


if __name__ == '__main__':
  import unittest
  unittest.main()

