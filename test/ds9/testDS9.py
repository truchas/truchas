#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class DS9(TruchasTest.GoldenTestCase):

  test_name = 'ds9'
  num_procs = 4 # with a parallel executable
  
  # Override the default setUp, omitting the opening of the golden output
  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass() # This runs Truchas
    self.test_output = Truchas.TruchasOutput(self.get_output_file())

  def test_final_temp(self):
    '''DS9: verifying the final temperature field'''
    tol = 1.0e-10
    coord = self.test_output.get_mesh().coordinates()
    cnode = self.test_output.get_mesh().connectivity()
    print coord.shape
    print cnode.shape
    # compute cell centroids using the coord and cnode arrays.
    # TODO: want to check that the final cycle number is 67

if __name__ == '__main__':
  import unittest
  unittest.main()

