#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class mytest(TruchasTest.GoldenTestCase):

  test_name = 'free-surf-flow-2'
  num_procs = 4 # with a parallel executable

  # Override the default setUp, omitting the opening of the golden output
  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass() # This runs Truchas
    self.test_output = Truchas.TruchasOutput(self.get_output_file())

  def test_fields(self):
    success = True
    
    # Target tolerances -- old flow passes with these
    # Initial conditions
    success &= self.vof_test(1, 1e-8)
    success &= self.pressure_test(1, 1e-14)
    success &= self.velocity_test(1, 1e-14)
    
    # Intermediate time
    success &= self.vof_test(2, 1e-8)
    success &= self.pressure_test(2, 1e-14)
    success &= self.velocity_test(2, 1e-14)
    
    # Final time
    success &= self.vof_test(3, 1e-8)
    success &= self.pressure_test(3, 1e-14)
    success &= self.velocity_test(3, 1e-14)
    
    self.assertTrue(success)

  def vof_test(self, id, tol):

    # The centroids function does not serialize, so we don't want to here either.
    time = self.test_output.get_simulation().find_series(id).time
    test = self.test_output.get_simulation().find_series(id).get_data('VOF',serialize=False)

    # Analytic vof solution at cell centroids
    cc = self.test_output.get_mesh().centroids()
    p = 2 - time
    vof = numpy.empty_like(test[:,0])
    for j in range(vof.size):
      x = cc[j,0]
      if x < p-0.1:
        vof[j] = 1
      elif x > p+0.1:
        vof[j] = 0
      else:
        vof[j] = 5*(p-(x-0.1))
    
    error = numpy.amax(abs(test[:,0]-vof))
    return self.report('vof', time, error, tol)

  def pressure_test(self, id, tol):
    time = self.test_output.get_simulation().find_series(id).time
    test = self.test_output.get_simulation().find_series(id).get_data('Z_P')
    error = numpy.amax(abs(test))
    return self.report('pressure', time, error, tol)
  
  def velocity_test(self, id, tol):
    time = self.test_output.get_simulation().find_series(id).time
    test = self.test_output.get_simulation().find_series(id).get_data('Z_VC')
    vof = self.test_output.get_simulation().find_series(id).get_data('VOF')[:,0]

    # the x-velocity is -1 in cells containing fluid
    uerror = max(abs(u + 1.) if vf > 0.0 else abs(u) for u,vf in zip(test[:,0],vof))
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

