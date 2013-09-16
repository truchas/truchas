#!/usr/bin/env python

import sys
import os

import unittest
import numpy

import Truchas
import TruchasTest

class NaturalConvection(TruchasTest.GoldenTestCase):

  test_name = 'natural-convection'
  num_procs = 4 # with a parallel executable
  restart_file = 'restart.bin'

  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())

  def runTest(self):
    '''Verify velocity at probes'''
    
    tol = 0.015
    fail = 0
    
    data = self.test_output.get_simulation().get_probe('VhMax','VC').get_data()
    test = data[-1,2] # final horizontal velocity: order is (cycle#,time,vx,vy,vz)
    gold = 7.585e-5
    error = abs((test - gold)/gold)
    if (error > tol):
      fail += 1
      print 'max horizontal velocity rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'max horizontal velocity rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)
    
    data = self.test_output.get_simulation().get_probe('VvMax','VC').get_data()
    test = data[-1,4] # final vertical velocity: order is (cycle#,time,vx,vy,vz)
    gold = 7.685e-5
    error = abs((test - gold)/gold)
    if (error > tol):
      fail += 1
      print 'max vertical velocity rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'max vertical velocity rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)
    
    self.assertTrue(fail==0)
    

if __name__ == '__main__':
  import unittest
  unittest.main()
