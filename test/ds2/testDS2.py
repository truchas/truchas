#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class DS2(TruchasTest.GoldenTestCase):

  test_name = 'ds2'
  num_procs = 4 # with a parallel executable

  def runTest(self):
    '''DS2: verifying the final temperature field'''

    # TODO: want to check that the final cycle number is 209

    tol = 1.0e-10
    T    = self.test_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    Tref = self.gold_output.get_simulation().find_series(id=2).get_data('Z_TEMP')
    error = max(abs(T-Tref)/Tref)
    self.assertTrue(error <= tol)

if __name__ == '__main__':
  import unittest
  unittest.main()

