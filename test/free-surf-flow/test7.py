#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class mytest(TruchasTest.GoldenTestCase):

  test_name = 'free-surf-flow-7'
  num_procs = 4 # with a parallel executable

  def test_fields(self):
    success = True
    
    # Final time
    success &= self.vof_test(2, 4e-2)
    success &= self.pressure_test(2, 1e-14)
    success &= self.velocity_test(2, 1e-13)
    
    self.assertTrue(success)

  def vof_test(self, id, tol):
    time = self.test_output.get_simulation().find_series(id).time
    test = self.test_output.get_simulation().find_series(id).get_data('VOF')[:,0]
    gold = self.gold_output.get_simulation().find_series(1).get_data('VOF')[:,0]
    error = numpy.amax(abs(test-gold))
    return self.report('vof', time, error, tol)

  def pressure_test(self, id, tol):
    time = self.test_output.get_simulation().find_series(id).time
    test = self.test_output.get_simulation().find_series(id).get_data('Z_P')
    gold = self.gold_output.get_simulation().find_series(id).get_data('Z_P')
    error = numpy.amax(abs(test))
    return self.report('pressure', time, error, tol)
  
  def velocity_test(self, id, tol):
    time = self.test_output.get_simulation().find_series(id).time
    test = self.test_output.get_simulation().find_series(id).get_data('Z_VC')
    vof = self.test_output.get_simulation().find_series(id).get_data('VOF')[:,0]

    # the x-velocity is 1 in cells containing fluid
    uerror = max(abs(u - (2.0*time/3.0)) if vf > 0.0 else abs(u) for u,vf in zip(test[:,0],vof))
    verror = max(abs(test[:,1]))
    werror = max(abs(test[:,2]))

    success = self.report('x-velocity', time, uerror, tol)
    success &= self.report('y-velocity', time, verror, tol)
    success &= self.report('z-velocity', time, werror, tol)
    return success

  def report(self, var, time, error, tol):
    success = error <= tol
    status = 'PASS' if success else 'FAIL'
    print '%s: %s at t=%8.2e: max error=%8.2e (tol=%8.2e)'%(status,var,time,error,tol)
    return success

if __name__ == '__main__':
  import unittest
  unittest.main()

