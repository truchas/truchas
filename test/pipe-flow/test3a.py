#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class mytest(TruchasTest.GoldenTestCase):

  test_name = 'pipe-flow-3a'
  num_procs = 4 # with a parallel executable

  # Override the default setUp, omitting the opening of the golden output
  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass() # This runs Truchas
    self.test_output = Truchas.TruchasOutput(self.get_output_file())
    self.test_sim = self.test_output.get_simulation()
    self.flow_region = Truchas.TruchasRegion(self.test_sim,[1])

  def test_fields(self):
    success = True
    success &= self.pressure_test(1,1e-13)
    success &= self.pressure_test(2,1e-11)
    success &= self.velocity_test(2,1.3e-3,1e-13)
    self.assertTrue(success)

  def velocity_test(self, id, tol1, tol2):
    time = self.test_sim.find_series(id).time
    data = self.test_sim.find_series(id).get_data('Z_VC',serialize=False,region=self.flow_region)

    centroids = self.test_output.get_mesh().centroids(region=self.flow_region)
    u = (1 - centroids[:,1]**2)/2

    error = max(abs(data[:,0] - u))
    success = self.report('x-velocity', time, error, tol1)
    error = max(abs(data[:,1]))
    success &= self.report('y-velocity', time, error, tol2)
    error = max(abs(data[:,2]))
    success &= self.report('y-velocity', time, error, tol2)
    return success

  def pressure_test(self, id, tol):
    data = self.test_sim.find_series(id).get_data('Z_P',serialize=False,region=self.flow_region)
    time = self.test_sim.find_series(id).time
    cc = self.test_output.get_mesh().centroids(region=self.flow_region)
    p = 3*(0.5 - cc[:,0])
    error = numpy.amax(abs(data-p))
    return self.report('pressure', time, error, tol)

  def report(self, var, time, error, tol):
    success = error <= tol
    status = 'PASS' if success else 'FAIL'
    print '%s: %s at t=%8.2e: max error=%8.2e (tol=%8.2e)'%(status,var,time,error,tol)
    return success

if __name__ == '__main__':
  import unittest
  unittest.main()

