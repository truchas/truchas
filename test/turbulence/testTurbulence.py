#!/usr/bin/env python

import sys
import os

import unittest
import numpy

import Truchas
import TruchasTest

class Turbulence(TruchasTest.GoldenTestCase):

  test_name = 'turbulence'
  num_procs = 4 # with a parallel executable

  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())

    # Empirically determined velocities :( MAC
    # The pressure values are from the linear pressure profile
    # The min/max velocity values are essentially "golden" values
    self.umin =  632.876
    self.umax = 2062.680
    self.pmin = 25000.0
    self.pmax = 75000.0
    
  def test_umin_velocity(self):
    '''Verify velocity at umin probe'''
    umin = self.test_output.get_simulation().get_probe('umin','VC').get_data()
    umin = umin[-1,2] # final x-velocity: order is (cycle#,time,vx,vy,vz)
    error = abs((umin-self.umin)/self.umin)
    tol = 1.2e-6
    if (error > tol):
      print '@umin: x-velocity rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print '@umin: x-velocity rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)
    
  def test_umax_velocity(self):
    '''Verify velocity at umax probe'''
    umax = self.test_output.get_simulation().get_probe('umax','VC').get_data()
    umax = umax[-1,2] # final x-velocity: order is (cycle#,time,vx,vy,vz)
    error = abs((umax-self.umax)/self.umax)
    tol = 1.2e-6
    if (error > tol):
      print '@umax: x-velocity rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print '@umax: x-velocity rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)
    
  def test_pdown_pressure(self):
    '''Verify pressure at pdown probe'''
    pmin = self.test_output.get_simulation().get_probe('pdown','P').get_data()
    pmin = pmin[-1,2] # final pressure: order is (cycle#,time,p)
    error = abs((pmin-self.pmin)/self.pmin)
    tol = 1.0e-9
    if (error > tol):
      print '@pdown: pressure rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print '@pdown: pressure rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)
    
  def test_umax_pressure(self):
    '''Verify pressure at umax probe'''
    pmax = self.test_output.get_simulation().get_probe('umax','P').get_data()
    pmax = pmax[-1,2] # final pressure: order is (cycle#,time,p)
    error = abs((pmax-self.pmax)/self.pmax)
    tol = 1.0e-9
    if (error > tol):
      print '@umax: pressure rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print '@umax: pressure rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)
    

if __name__ == '__main__':
  import unittest
  unittest.main()
