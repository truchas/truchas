#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class DS7(TruchasTest.GoldenTestCase):

  test_name = 'ds7'
  num_procs = 4 # with a parallel executable

  def test_early_temperature(self):
    '''DS7: verifying the early temperature field'''

    # TODO: want to check that the final cycle number is 21

    tol = 1.0e-10
    test = self.test_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    gold = self.gold_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    error = max(abs(test-gold)/gold)
    self.assertTrue(error <= tol)

  def test_final_temperature(self):
    '''DS7: verifying the final temperature field'''

    # TODO: want to check that the final cycle number is 48

    tol = 1.0e-10
    test = self.test_output.get_simulation().find_series(id=3).get_data('Z_TEMP')
    gold = self.gold_output.get_simulation().find_series(id=3).get_data('Z_TEMP')
    error = max(abs(test-gold)/gold)
    self.assertTrue(error <= tol)

if __name__ == '__main__':
  import unittest
  unittest.main()

