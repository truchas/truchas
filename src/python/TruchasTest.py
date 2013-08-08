#!/usr/bin/env python

import os
import sys

import logging
import unittest


class TruchasTestCase(unittest.TestCase):

  _test_is_initialized = False

  @classmethod
  def setUpClass(cls):
    pass

  @classmethod
  def tearDownClass(cls):
    pass

  def setUp(self):

    if self._test_is_initialized is False:
      self.setUpClass()


  def tearDown(self):
    pass


