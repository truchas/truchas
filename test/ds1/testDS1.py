#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class DS1(TruchasTest.GoldenTestCase):

  test_name = 'ds1'
  num_procs = 4 # with a parallel executable

  def runTest(self):
    '''DS1: verifying the final concentration field'''

    # TODO: want to check that the final cycle number is 209

    tol = 1.0e-10
    C    = self.test_output.get_simulation().find_series(id=2).get_data('phi1')
    Cref = self.gold_output.get_simulation().find_series(id=2).get_data('phi1')
    error = max(abs(C-Cref)/Cref)
    self.assertTrue(error <= tol)

if __name__ == '__main__':
  import unittest
  unittest.main()

