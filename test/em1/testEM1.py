#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class EM1(TruchasTest.GoldenTestCase):

  test_name='em1'

  num_procs = 4

  free_space_ids = 10
  free_space_region = None

  conducting_ids = [11, 12]
  conducting_region = None

  # Run before each test, must call the base class setUp
  # to ensure all the truchas binary instance is set up correctly
  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())
    self.gold_output=Truchas.TruchasOutput(self.get_golden_output_file())
    self.test_sim = self.test_output.get_simulation()
    self.test_simulation=self.test_output.get_simulation()
    self.free_space_region = Truchas.TruchasRegion(self.test_simulation,self.free_space_ids)
    self.conducting_region = Truchas.TruchasRegion(self.test_simulation,self.conducting_ids)

  def get_test_field(self,field,cycle,region=None):
    return self.test_output.get_simulation().find_series(cycle=cycle).get_data(field,region=region)

  def get_gold_field(self,field,cycle,region=None):
    return self.gold_output.get_simulation().find_series(cycle=cycle).get_data(field,region=region)

  
  def test_initial_joule_heat(self):
    '''Test the initial joule heat'''
   
    tol = 1.0e-8
    fail = 0

    # For cycle = 0
    test = self.get_test_field('Joule_P',cycle=0,region=self.conducting_region)
    gold = self.get_gold_field('Joule_P',cycle=0,region=self.conducting_region)
    error=max(abs((test-gold)/gold))
    if error > tol:
      fail += 1
      print 'initial joule heat: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'initial joule heat: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    # For cycles 1-5 should be identical to the initial joule heat
    for j in range(1,6):
      next = self.get_test_field('Joule_P',cycle=j,region=self.conducting_region)
      error=max(abs(next-test))
      if error != 0:
        fail += 1
        print 'unexpected change in the Joule heat for cycle %d: FAIL'%(j)

    self.assertTrue(fail==0)

  
  def xtest_scaled_joule_heat(self):
    '''Test the scaled joule heat'''
   
    n = 6
    tol = 1.0e-15
    fail = 0

    # For cycle = 6
    gold = 4 * self.get_test_field('Joule_P',cycle=0,region=self.conducting_region)
    test = self.get_test_field('Joule_P',cycle=n,region=self.conducting_region)
    error=max(abs((test-gold)/gold))
    if error > tol:
      fail += 1
      print 'cycle %2d: joule heat: max rel error = %8.2e: FAIL (tol=%8.2e)'%(n,error,tol)
    else:
      print 'cycle %2d: joule heat: max rel error = %8.2e: PASS (tol=%8.2e)'%(n,error,tol)

    # For cycles 7-10 should be identical to the cycle 6 joule heat
    for j in range(7,11):
      next = self.get_test_field('Joule_P',cycle=j,region=self.conducting_region)
      error=max(abs(next-test))
      if error != 0:
        fail += 1
        print 'unexpected change in the Joule heat for cycle %d: FAIL'%(j)

    self.assertTrue(fail==0)


  def xtest_zero_joule_heat(self):
    '''Verify the zero joule heat'''
    fail=0
    for j in range(11,16):
      test = self.get_test_field('Joule_P',cycle=j)
      error = max(abs(test))
      if error != 0.0:
	fail += 1
        print 'unexpected non-zero Joule heat for cycle %d: FAIL'%(j)
    self.assertTrue(fail==0)	
      

  def xtest_final_joule_heat(self):
    '''Verify the final joule heat'''

    n = 16
    tol = 1.0e-8
    fail = 0

    # For cycle = 16
    test = self.get_test_field('Joule_P',cycle=n,region=self.conducting_region)
    gold = self.get_gold_field('Joule_P',cycle=n,region=self.conducting_region)
    error=max(abs((test-gold)/gold))
    if error > tol:
      fail += 1
      print 'cycle %2d: joule heat: max rel error = %8.2e: FAIL (tol=%8.2e)'%(n,error,tol)
    else:
      print 'cycle %2d: joule heat: max rel error = %8.2e: PASS (tol=%8.2e)'%(n,error,tol)

    # For cycles 17-20 should be identical to the cycle 16 joule heat
    for j in range(17-21):
      next = self.get_test_field('Joule_P',cycle=j,region=self.conducting_region)
      error=max(abs(next-test))
      if error != 0:
        fail += 1
        print 'unexpected change in the Joule heat for cycle %d: FAIL'%(j)

    self.assertTrue(fail==0)


  def xtest_free_space_joule_heat(self):
    '''Verify no free-space Joule heat'''
    fail = 0
    for j in range(0,21):
      test = self.get_test_field('Joule_P',cycle=j,region=self.free_space_region)
      error = max(abs(test))
      if error != 0.0:
	fail += 1
	print '  unexpected non-zero free-space Joule heat for cycle %2d: FAIL'%(j)
    self.assertTrue(fail==0)

  def xtest_final_temperature(self):
    '''Verify the final temperature field'''
    tol = 1.0e-8
    test = self.get_test_field('Z_TEMP',cycle=20,region=self.conducting_region)
    gold = self.get_gold_field('Z_TEMP',cycle=20,region=self.conducting_region)
    error = max(abs((test-gold)/gold))
    if error > tol:
      print 'temperature: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'temperature: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)


if __name__ == '__main__':
  import unittest
  unittest.main()

