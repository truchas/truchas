#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class mytest(TruchasTest.GoldenTestCase):

  test_name = 'free-surf-flow-8'
  num_procs = 4 # with a parallel executable

  def test_fields(self):
    success = True

    for i in range(1,3):
      success &= self.vof_test(i, 1e-13)
      success &= self.pressure_test(i, 1e-12)
      success &= self.velocity_test(i, 1e-12)

    self.assertTrue(success)

  def vof_test(self, id, tol):
    time = self.test_output.get_simulation().find_series(id).time
    test = self.test_output.get_simulation().find_series(id).get_data('VOF')[:,0]
    gold = self.gold_output.get_simulation().find_series(id).get_data('VOF')[:,0]
    error = numpy.amax(abs(test-gold))
    return self.report('vof', time, error, tol)

  def pressure_test(self, id, tol):
    time = self.test_output.get_simulation().find_series(id).time
    test = self.test_output.get_simulation().find_series(id).get_data('Z_P')
    gold = self.gold_output.get_simulation().find_series(id).get_data('Z_P')

    error = numpy.amax(abs(test-gold))
    success = self.report('pressure', time, error, tol)

    # ensure pressure is zero in the void
    vof = self.test_output.get_simulation().find_series(id).get_data('VOF')[:,0]
    void_error = max(abs(p) for p,vf in zip(test,vof) if vf == 0)
    success &= self.report('void-pressure', time, void_error, tol)

    return success

  def velocity_test(self, id, tol):
    time = self.test_output.get_simulation().find_series(id).time
    test = self.test_output.get_simulation().find_series(id).get_data('Z_VC')
    gold = self.gold_output.get_simulation().find_series(id).get_data('Z_VC')
    vof = self.test_output.get_simulation().find_series(id).get_data('VOF')[:,0]

    # compare with golden output
    uerror = max(abs(test[:,0] - gold[:,0]))
    verror = max(abs(test[:,1] - gold[:,1]))
    werror = max(abs(test[:,2] - gold[:,2]))

    success = self.report('x-velocity', time, uerror, tol)
    success &= self.report('y-velocity', time, verror, tol)
    success &= self.report('z-velocity', time, werror, tol)

    # the velocity is 0 in purely void cells
    error = max(max(abs(u)) for u,vf in zip(test,vof) if vf == 0)
    success &= self.report('void-velocity', time, error, tol)

    return success

  def report(self, var, time, error, tol):
    success = error <= tol
    status = 'PASS' if success else 'FAIL'
    print '%s: %s at t=%8.2e: max error=%8.2e (tol=%8.2e)' % (status,var,time,error,tol)
    return success

if __name__ == '__main__':
  import unittest
  unittest.main()
