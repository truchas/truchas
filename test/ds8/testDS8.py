#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class DS8(TruchasTest.GoldenTestCase):

  test_name = 'ds8'
  num_procs = 4 # with a parallel executable

  def test_early_temp(self):
    '''DS8: verifying the temperature field at early time'''
    # TODO: want to check that the final cycle number is 94
    tol = 1.0e-10
    test = self.test_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    gold = self.gold_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    error = max(abs(test-gold)/gold)
    self.assertTrue(error <= tol)

  def test_final_temp(self):
    '''DS8: verifying the temperature field at final time'''
    # TODO: want to check that the final cycle number is 290
    tol = 1.0e-10
    test = self.test_output.get_simulation().find_series(id=4).get_data('Z_TEMP')
    gold = self.gold_output.get_simulation().find_series(id=4).get_data('Z_TEMP')
    error = max(abs(test-gold)/gold)
    self.assertTrue(error <= tol)

  def test_early_velocity(self):
    '''DS8: verifying the velocity field at early time'''
    tol = 1.0e-10
    test = self.test_output.get_simulation().find_series(id=2).get_data('Z_VC')
    gold = self.gold_output.get_simulation().find_series(id=2).get_data('Z_VC')
    uerror = max(abs(test[:,0]-gold[:,0]))
    verror = max(abs(test[:,1]-gold[:,1]))
    werror = max(abs(test[:,2]-gold[:,2]))
    error = max(uerror,verror)
    self.assertTrue(error <= tol)

  def test_final_velocity(self):
    '''DS8: verifying the velocity field at final time'''
    tol = 1.0e-10
    test = self.test_output.get_simulation().find_series(id=4).get_data('Z_VC')
    gold = self.gold_output.get_simulation().find_series(id=4).get_data('Z_VC')
    uerror = max(abs(test[:,0]-gold[:,0]))
    verror = max(abs(test[:,1]-gold[:,1]))
    werror = max(abs(test[:,2]-gold[:,2]))
    error = max(uerror,verror)
    self.assertTrue(error <= tol)

if __name__ == '__main__':
  import unittest
  unittest.main()

