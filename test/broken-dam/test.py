#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class BrokenDam(TruchasTest.GoldenTestCase):

  test_name = 'broken-dam'
  num_procs = 4 # with a parallel executable
  
  def test_final_cycle_number(self):
    '''BROKEN DAM: checking the final cycle number'''
    test_series = self.test_output.get_simulation().find_series(id=2)
    gold_series = self.gold_output.get_simulation().find_series(id=2)
    self.assertTrue(test_series.cycle == gold_series.cycle)

  def test_final_fluid_frac(self):
    '''BROKEN DAM: verifying the fluid volume fraction at final time'''
    time = self.test_output.get_simulation().find_series(id=2).time
    test = self.test_output.get_simulation().find_series(id=2).get_data('VOF')
    gold = self.gold_output.get_simulation().find_series(id=2).get_data('VOF')
    error = max(abs(test[:,0]-gold[:,0])) # comp 0 is fluid
    tol = 1.0e-9
    if error > tol:
      print 'vof at t=%8.2e: max error = %8.2e: FAIL (tol=%8.2e)'%(time,error,tol)
      self.assertTrue(False)
    else:
      print 'vof at t=%8.2e: max error = %8.2e: PASS (tol=%8.2e)'%(time,error,tol)

  def test_final_velocity(self):
    '''BROKEN DAM: verifying the velocity field at final time'''
    time = self.test_output.get_simulation().find_series(id=2).time
    test = self.test_output.get_simulation().find_series(id=2).get_data('Z_VC')
    gold = self.gold_output.get_simulation().find_series(id=2).get_data('Z_VC')
    uerror = max(abs(test[:,0]-gold[:,0]))
    verror = max(abs(test[:,1]-gold[:,1]))
    error = max(uerror,verror)
    tol = 1.0e-9
    if error > tol:
      print 'velocity at t=%8.2e: max error = %8.2e: FAIL (tol=%8.2e)'%(time,error,tol)
      self.assertTrue(False)
    else:
      print 'velocity at t=%8.2e: max error = %8.2e: PASS (tol=%8.2e)'%(time,error,tol)

if __name__ == '__main__':
  import unittest
  unittest.main()

