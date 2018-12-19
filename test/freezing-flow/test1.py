#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class FreezingFlowTest(TruchasTest.GoldenTestCase):

  test_name = 'freezing-flow-1'
  num_procs = 4 # with a parallel executable

  def test_fields(self):
    success = True

    # first output dump
    success &= self.temp_test(2, 5e-7)
    success &= self.pressure_test(2, 1e-10)
    success &= self.velocity_test(2, 1e-10)
    success &= self.vof_test(2, 5e-7)

    # final time
    success &= self.temp_test(3, 5e-6)
    success &= self.pressure_test(3, 1e-10)
    success &= self.velocity_test(3, 1e-10)
    success &= self.vof_test(3, 1e-10)

    self.assertTrue(success)

  def temp_test(self, timeid, tol):
    time = self.test_output.get_simulation().find_series(id=timeid).time
    test = self.test_output.get_simulation().find_series(id=timeid).get_data('Z_TEMP')
    gold = self.gold_output.get_simulation().find_series(id=timeid).get_data('Z_TEMP')
    error = max(abs(test - gold))
    self.report('temperature', time, error, tol)
    return error <= tol

  def vof_test(self, timeid, tol):
    time = self.test_output.get_simulation().find_series(id=timeid).time
    test = self.test_output.get_simulation().find_series(id=timeid).get_data('VOF')[:,0]
    gold = self.gold_output.get_simulation().find_series(id=timeid).get_data('VOF')[:,0]
    error = max(abs(test - gold))
    self.report('VOF', time, error, tol)
    return error <= tol

  def pressure_test(self, timeid, tol):
    time = self.test_output.get_simulation().find_series(id=timeid).time
    test = self.test_output.get_simulation().find_series(id=timeid).get_data('Z_P')
    error = max(abs(test))
    self.report('pressure', time, error, tol)
    return error <= tol

  def velocity_test(self, timeid, tol):
    time = self.test_output.get_simulation().find_series(id=timeid).time
    test = self.test_output.get_simulation().find_series(id=timeid).get_data('Z_VC')

    uerror = max(abs(test[:,0]))
    verror = max(abs(test[:,1]))
    werror = max(abs(test[:,2]))

    self.report('uvel', time, uerror, tol)
    self.report('vvel', time, verror, tol)
    self.report('wvel', time, werror, tol)

    return uerror <= tol and verror <= tol and werror <= tol

  def report(self, var, time, err, tol):
    status = 'PASS' if err <= tol else 'FAIL'
    print '%s: %s at t=%8.2e max error=%8.2e (tol=%8.2e)' % (status, var, time, err, tol)

if __name__ == '__main__':
  import unittest
  unittest.main()
