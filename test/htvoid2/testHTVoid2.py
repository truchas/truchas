#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class HTVoid2(TruchasTest.GoldenTestCase):

  test_name = 'htvoid2'
  num_procs = 4 # with a parallel executable

  # TODO: want to check that the final cycle number is 49

  def test_final_temp(self):
    '''HTVOID2: verifying the temperature field at final time'''
    tol = 1.0e-6
    test = self.test_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    gold = self.gold_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    error = max(abs(test-gold))
    self.assertTrue(error <= tol)

if __name__ == '__main__':
  import unittest
  unittest.main()

