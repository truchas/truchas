#!/usr/bin/env python

import sys
import os

import numpy
import math

import Truchas
import TruchasTest

class mytest(TruchasTest.GoldenTestCase):

  test_name = 'inviscid-pipe-flow-old-2d'
  num_procs = 4 # with a parallel executable

  # Override the default setUp, omitting the opening of the golden output
  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass() # This runs Truchas
    self.test_output = Truchas.TruchasOutput(self.get_output_file())
    self.test_sim = self.test_output.get_simulation()
    self.flow_region = Truchas.TruchasRegion(self.test_sim,[1])

  def velocity_test(self, id, tol):
    time = self.test_sim.find_series(id).time
    data = self.test_sim.find_series(id).get_data('Z_VC',region=self.flow_region)

    fail = 0

    vel = math.sqrt(2)*time
    error = max(abs((data[:,0] - vel)/vel))
    if (error > tol):
      fail += 1
      print 'x-velocity at t=%8.2e: max rel error = %8.2e: FAIL (tol=%8.2e)'%(time,error,tol)
    else:
      print 'x-velocity at t=%8.2e: max rel error = %8.2e: PASS (tol=%8.2e)'%(time,error,tol)

    error = max(abs((data[:,2] - vel)/vel))
    if (error > tol):
      fail += 1
      print 'y-velocity at t=%8.2e: max rel error = %8.2e: FAIL (tol=%8.2e)'%(time,error,tol)
    else:
      print 'y-velocity at t=%8.2e: max rel error = %8.2e: PASS (tol=%8.2e)'%(time,error,tol)

    error = max(abs(data[:,2]))
    if (error > tol):
      fail += 1
      print 'z-velocity at t=%8.2e: max error = %8.2e: FAIL (tol=%8.2e)'%(time,error,tol)
    else:
      print 'z-velocity at t=%8.2e: max error = %8.2e: PASS (tol=%8.2e)'%(time,error,tol)

    self.assertTrue(fail == 0)
    
  def test_velocity2(self):
    '''verify early velocity'''
    self.velocity_test(2,1e-14)

  def test_velocity3(self):
    '''verify final velocity'''
    self.velocity_test(3,1e-14)

  def pressure_test(self, id, tol):
    data = self.test_sim.find_series(id).get_data('Z_P',serialize=False,region=self.flow_region)
    time = self.test_sim.find_series(id).time
    cc = self.test_output.get_mesh().centroids(region=self.flow_region)
    p = 6*(0.5 - (cc[:,0]+cc[:,1])/math.sqrt(2))
    error = numpy.amax(abs(data-p))
    if error > tol:
      print 'pressure at t=%8.2e: max error = %8.2e: FAIL (tol=%8.2e)'%(time,error,tol)
      self.assertTrue(False)
    else:
      print 'pressure at t=%8.2e: max error = %8.2e: PASS (tol=%8.2e)'%(time,error,tol)

  def test_pressure1(self):
    '''verify initial pressure'''
    self.pressure_test(1,1e-13)

  def test_pressure2(self):
    '''verify final pressure'''
    self.pressure_test(3,1e-13)

if __name__ == '__main__':
  import unittest
  unittest.main()

