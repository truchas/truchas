#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class EM3(TruchasTest.GoldenTestCase):

  test_name = 'em3'
  num_procs = 4 # with a parallel executable
  
  def get_test_field(self,field,cycle,region=None):
    return self.test_output.get_simulation().find_series(cycle=cycle).get_data(field,region=region)

  def get_gold_field(self,field,cycle,region=None):
    return self.gold_output.get_simulation().find_series(cycle=cycle).get_data(field,region=region)

  def test_joule_heat(self):
    '''Verify the Joule heat fields'''
    tol = 1.0e-6
    test = self.get_test_field('Joule_P',cycle=0)
    gold = self.get_gold_field('Joule_P',cycle=0)
    error = max(abs((test-gold)/gold))
    if error > tol:
      print 'joule heat: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'joule heat: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)
    

  def test_final_temperature(self):
    '''Verify the final temperature field'''
    tol = 1.0e-6
    test = self.get_test_field('Z_TEMP',cycle=30)
    gold = self.get_gold_field('Z_TEMP',cycle=30)
    error = max(abs((test-gold)/gold))
    if error > tol:
      print 'temperature: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'temperature: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

if __name__ == '__main__':
  import unittest
  unittest.main()

