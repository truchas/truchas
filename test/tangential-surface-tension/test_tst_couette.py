#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class TangentialSurfaceTension(TruchasTest.GoldenTestCase):

  test_name = 'tangential-surface-tension-couette'
  num_procs = 4 # with a parallel executable

  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())

  def test_final_temp(self):
    '''tangential-surface-tension: verifying the temperature field at early time'''
    tol = 1e-8
    # The centroids function does not serialize, so we don't want to here either.
    test = self.test_output.get_simulation().find_series(id=2).get_data('Z_TEMP', serialize=False)
    cc = self.test_output.get_mesh().centroids()
    minx = min(cc[:,0]) - 0.5
    tex = 1 + (cc[:,0] - minx)
    error = max(abs(test - tex))
    self.report('final temperature', error, tol)
    self.assertTrue(error <= tol)

  def test_final_pressure(self):
    '''tangential-surface-tension: verifying the temperature field at early time'''
    tol = 1e-8
    test = self.test_output.get_simulation().find_series(id=2).get_data('Z_P')
    error = max(abs(test))
    self.report('final pressure', error, tol)
    self.assertTrue(error <= tol)

  def test_final_velocity(self):
    '''tangential-surface-tension: verifying the velocity field at final time'''
    tol = 1e-7
    tolt = 1e-10
    # The centroids function does not serialize, so we don't want to here either.
    test = self.test_output.get_simulation().find_series(id=2).get_data('Z_VC', serialize=False)

    dsig_dx = -1.0
    viscosity = 20.0
    zc = self.test_output.get_mesh().centroids()[:,2]
    minz = min(zc) - 0.5

    uerror = max(abs(test[:,0] - dsig_dx / viscosity * (zc - minz)))
    verror = max(abs(test[:,1]))
    werror = max(abs(test[:,2]))

    self.report('final uvel', uerror, tol)
    self.report('final vvel', verror, tolt)
    self.report('final wvel', werror, tolt)

    self.assertTrue(uerror <= tol)
    self.assertTrue(verror <= tolt)
    self.assertTrue(werror <= tolt)

  def report(self, var, err, tol):
    status = 'PASS' if err <= tol else 'FAIL'
    print '%s max error=%8.2e: %s (tol=%8.2e)' % (var, err, status, tol)

if __name__ == '__main__':
  import unittest
  unittest.main()
