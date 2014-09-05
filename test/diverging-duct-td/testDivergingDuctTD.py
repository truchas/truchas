#!/usr/bin/env python

import sys
import os

import unittest
import numpy

import Truchas
import TruchasTest

class DivergingDuctTD(TruchasTest.GoldenTestCase):

  test_name = 'diverging-duct-td'
  num_procs = 4 # with a parallel executable

  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())

  def get_test_field(self,field,cycle,serialize=True,region=None):
    return self.test_output.get_simulation().find_series(cycle=cycle).get_data(field,serialize,region=region)

  def get_test_time(self,cycle):
    return self.test_output.get_simulation().find_series(cycle=cycle).time

  def vin(self,t):
    if t <= 5.0:
      return 1.0
    elif t >= 15.0:
      return 0.5
    else:
      return 1.0 - 0.05*(t - 5.0)

  def velocity(self,x,t):
    return self.vin(t)/(1.0 + x/40.0)

  def pressure(self,x,t):
    return 0.68 + 0.5*self.vin(t)**2 * (0.64 - 1.0/(1.0 + x/40.0)**2)

  def test_pressure_0(self):
    '''Verify initial pressure field against expected analytic value'''
    
    n = 207 #141 #115
    tol = 3.0e-3
    # The centroids function does not serialize, so we don't want to here either.
    test = self.get_test_field('Z_P',cycle=n,serialize=False)
    time = self.get_test_time(cycle=n)
    
    xc = self.test_output.get_mesh().centroids()
    gold = self.pressure(xc[:,0],time)
    print max(gold), min(gold)
    error = max(abs((test-gold)/gold))
    if error > tol:
      print 'pressure: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'pressure: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

  def test_pressure_1(self):
    '''Verify final pressure field against expected analytic value'''
    
    n = 647 #361 #225
    tol = 1.0e-3
    # The centroids function does not serialize, so we don't want to here either.
    test = self.get_test_field('Z_P',cycle=n,serialize=False)
    time = self.get_test_time(cycle=n)
    
    xc = self.test_output.get_mesh().centroids()
    gold = self.pressure(xc[:,0],time)
    print max(gold), min(gold)
    error = max(abs((test-gold)/gold))
    if error > tol:
      print 'pressure: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'pressure: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

  def test_velocity_0(self):
    '''Verify initial velocity field against expected analytic value'''
    
    n = 207 #141 #115
    tol = 1.0e-3
    # The centroids function does not serialize, so we don't want to here either.
    test = self.get_test_field('Z_VC',cycle=n,serialize=False)
    time = self.get_test_time(cycle=n)
    
    xc = self.test_output.get_mesh().centroids()
    gold = self.velocity(xc[:,0],time)
    print max(gold), min(gold)
    error = max(abs((test[:,0]-gold)/gold))
    if error > tol:
      print 'velocity: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'velocity: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

  def test_velocity_1(self):
    '''Verify intermediate velocity field against expected analytic value'''
    
    n = 422 #248 #169
    tol = 3.0e-3
    # The centroids function does not serialize, so we don't want to here either.
    test = self.get_test_field('Z_VC',cycle=n,serialize=False)
    time = self.get_test_time(cycle=n)
    print time
    
    xc = self.test_output.get_mesh().centroids()
    gold = self.velocity(xc[:,0],time)
    print max(gold), min(gold)
    error = max(abs((test[:,0]-gold)/gold))
    if error > tol:
      print 'velocity: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'velocity: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

  def test_velocity_2(self):
    '''Verify final velocity field against expected analytic value'''
    
    n = 647 #361 #225
    tol = 1.0e-3
    # The centroids function does not serialize, so we don't want to here either.
    test = self.get_test_field('Z_VC',cycle=n,serialize=False)
    time = self.get_test_time(cycle=n)
    
    xc = self.test_output.get_mesh().centroids()
    gold = self.velocity(xc[:,0],time)
    print max(gold), min(gold)
    error = max(abs((test[:,0]-gold)/gold))
    if error > tol:
      print 'velocity: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'velocity: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)


if __name__ == '__main__':
  import unittest
  unittest.main()
