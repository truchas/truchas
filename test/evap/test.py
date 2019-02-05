#!/usr/bin/env python

import sys
import os

import numpy
import numpy.ma
import unittest

import Truchas
import TruchasTest

class EvapTest(TruchasTest.GoldenTestCase):

  test_name = 'evap'
  num_procs = 4 # with a parallel executable

  def get_test_field(self,field,id,region=None):
    return self.test_output.get_simulation().find_series(id=id).get_data(field,region=region)

  def get_gold_field(self,field,id,region=None):
    return self.gold_output.get_simulation().find_series(id=id).get_data(field,region=region)

  def test_temp1(self):
    '''EVAP: verifying the temperature field at t=5'''
    tol = 1.0e-7
    test = self.get_test_field('Z_TEMP',id=1)
    gold = self.get_gold_field('Z_TEMP',id=1)
    error = max(abs((test-gold)/gold))
    print 'temp max rel error =', error, '(tol =', tol, ')'
    self.assertTrue(error <= tol)

  def test_temp2(self):
    '''EVAP: verifying the temperature field at t=10'''
    tol = 1.0e-7
    test = self.get_test_field('Z_TEMP',id=2)
    gold = self.get_gold_field('Z_TEMP',id=2)
    error = max(abs((test-gold)/gold))
    print 'temp max rel error =', error, '(tol =', tol, ')'
    self.assertTrue(error <= tol)

  def test_temp3(self):
    '''EVAP: verifying the temperature field at t=15'''
    tol = 1.0e-7
    test = self.get_test_field('Z_TEMP',id=3)
    gold = self.get_gold_field('Z_TEMP',id=3)
    error = max(abs((test-gold)/gold))
    print 'temp max rel error =', error, '(tol =', tol, ')'
    self.assertTrue(error <= tol)

if __name__ == '__main__':
  import unittest
  unittest.main()
