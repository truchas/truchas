#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class EM1(TruchasTest.GoldenTestCase):

  test_name='ds11'

  num_procs=4

  def runTest(self):
    '''EM1 Test compare temperature and VOF fields'''

    fails = 0

    field='Z_TEMP'
    tol=1.0e-9
    for i in range(4):
      gold_data=self.gold_output.get_simulation().find_series(id=4).get_data(field)
      test_data=self.test_output.get_simulation().find_series(id=4).get_data(field)
      error=max(abs(test_data-gold_data)/gold_data)
      if error > tol:
	fails+=1

    field='VOF' 
    tol=1.0e-8
    for i in range(4):
      id=i+1
      print 'Compare sequence %d' % id
      gold_data=self.gold_output.get_simulation().find_series(id=id).get_data(field)
      test_data=self.test_output.get_simulation().find_series(id=id).get_data(field)
      error=max(abs(test_data[:,2]-gold_data[:,2]))
      if error > tol:
	fails+=1

    self.assertTrue(fails == 0)	

if __name__ == '__main__':
  import unittest
  unittest.main()


