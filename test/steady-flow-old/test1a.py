#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class mytest(TruchasTest.GoldenTestCase):

  test_name = 'steady-flow-old-1a'
  num_procs = 4 # with a parallel executable

  # Override the default setUp, omitting the opening of the golden output
  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass() # This runs Truchas
    self.test_output = Truchas.TruchasOutput(self.get_output_file())

  def vof_test(self, id, tol):

    # The centroids function does not serialize, so we don't want to here either.
    test = self.test_output.get_simulation().find_series(id).get_data('VOF',serialize=False)
    time = self.test_output.get_simulation().find_series(id).time

    # Analytic vof solution at cell centrioids
    cc = self.test_output.get_mesh().centroids()
    p = -4 + 4*time
    vof = numpy.empty_like(test[:,0])
    for j in range(vof.size):
      x = cc[j,0]
      if x < p-0.5:
        vof[j] = 1
      elif x > p+0.5:
        vof[j] = 0
      else:
        vof[j] = p-(x-0.5)
    
    error = numpy.amax(abs(test[:,0]-vof))
    if error > tol:
      print 'vof at t=%8.2e: max error = %8.2e: FAIL (tol=%8.2e)'%(time,error,tol)
      self.assertTrue(False)
    else:
      print 'vof at t=%8.2e: max error = %8.2e: PASS (tol=%8.2e)'%(time,error,tol)

  def test_vof_step(self):
    '''Verify vof after first step'''
    self.vof_test(2, 2e-13)

  def test_final_vof(self):
    '''Verify final vof'''
    self.vof_test(3, 2e-13)

  def test_final_velocity(self):
    '''Verify final velocity'''
    data = self.test_output.get_simulation().find_series(id=3).get_data('Z_VC')
    data[:,0] -= 4
    error = numpy.amax(abs(data))
    tol = 1.0e-13
    if error > tol:
      print 'velocity: max error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'velocity: max error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

  def pressure_test(self, id, tol):
    data = self.test_output.get_simulation().find_series(id).get_data('Z_P')
    time = self.test_output.get_simulation().find_series(id).time
    error = numpy.amax(abs(data))
    if error > tol:
      print 'pressure at t=%8.2e: max error = %8.2e: FAIL (tol=%8.2e)'%(time,error,tol)
      self.assertTrue(False)
    else:
      print 'pressure at t=%8.2e: max error = %8.2e: PASS (tol=%8.2e)'%(time,error,tol)

  def test_initial_pressure(self):
    '''Verify initial pressure'''
    self.pressure_test(1,1e-13)

  def test_final_pressure(self):
    '''Verify final pressure'''
    self.pressure_test(3,1e-13)

if __name__ == '__main__':
  import unittest
  unittest.main()

