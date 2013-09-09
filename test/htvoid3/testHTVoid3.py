#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class HTVoid3(TruchasTest.GoldenTestCase):

  test_name = 'htvoid3'
  num_procs = 4 # with a parallel executable

  # TODO: want to check that the final cycle number is 276

  def test_final_temp(self):
    '''HTVOID3: verifying the temperature field at final time'''
    tol = 1.0e-6
    test = self.test_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    gold = self.gold_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    error = max(abs(test-gold))
    self.assertTrue(error <= tol)

  def test_final_solid_frac(self):
    '''HTVOID3: verifying the solid volume fraction at final time'''
    tol = 1.0e-6
    test = self.test_output.get_simulation().find_series(id=2).get_data('VOF')
    gold = self.gold_output.get_simulation().find_series(id=2).get_data('VOF')
    error = max(abs(test[:,1]-gold[:,1])) # comp 1 is solid
    self.assertTrue(error <= tol)

if __name__ == '__main__':
  import unittest
  unittest.main()

