#!/usr/bin/env python

import sys
import os

import unittest

import Truchas
import TruchasTest

class mytest(TruchasTest.GoldenTestCase):

  test_name = 'advection-0'
  num_procs = 4 # with a parallel executable

  # Override the default setUp, omitting the opening of the golden output
  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())

  # Really nothing we can test here, other than it ran successfully

if __name__ == '__main__':
  import unittest
  unittest.main()
