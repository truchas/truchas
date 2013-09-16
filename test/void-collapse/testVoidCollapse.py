#!/usr/bin/env python

import sys
import os

import unittest
import numpy

import Truchas
import TruchasTest

class VoidCollapse(TruchasTest.GoldenTestCase):

  test_name = 'void-collapse'
  num_procs = 4 # with a parallel executable

  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())

  def get_test_field(self,field,cycle,region=None):
    return self.test_output.get_simulation().find_series(cycle=cycle).get_data(field,region=region)

# NNC, Sept 2013.  I've kept comments by DAK about the test for reference.
# I've migrated the test as it was.  I don't understant it at all.

  def foo_test_velocity(self):
    '''Verify velocity field'''

    n = 113
    fail = 0

    test = self.get_test_field('Z_VC',cycle=n)

    vx = 0.5
    tol = 5.0e-8
    error = max(abs((test[:,0]-vx)/vx))
    if (error > tol):
      fail += 1
      print 'x-velocity: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'x-velocity: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    tol = 5.0e-10
    for k in [1,2]:
      error = max(abs(test[:,k]))
      if (error > tol):
        fail += 1
        print '%s-velocity: max abs error = %8.2e: FAIL (tol=%8.2e)'%(('x','y','z')[k],error,tol)
      else:
        print '%s-velocity: max abs error = %8.2e: PASS (tol=%8.2e)'%(('x','y','z')[k],error,tol)

    self.assertTrue(fail == 0)


#  This test seems to check some numerical discretization transient that I do not understand. DAK 8-24-10
  def test_pressure_1(self):
    '''Verify pressure field at cycle 1'''

    tol = 1.0e-3

    cmap = self.test_output.get_simulation().cellmap()
    print cmap

    p = self.get_test_field('Z_P',cycle=1)
    print p

    vof = self.get_test_field('VOF',cycle=1)
    test = self.get_test_field('Z_P',cycle=1)
    print test
    test = numpy.ma.masked_array(test, mask=vof[:,0]<0.01).compressed()
    print test

    #P_anal = [66,56,46,36,26,16,6,0.5]
    P_anal = [65.9934,55.9944,45.9954,35.9964,25.9974,15.9984,5.9994,0.49995]
    gold = numpy.array(P_anal)

    error = max(abs(test - gold))
    if (error > tol):
      print 'cycle 1 pressure: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'cycle 1 pressure: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

  def test_pressure_2(self):
    '''Verify pressure field at cycle 2'''

    tol = 1.0e-3

    vof = self.get_test_field('VOF',cycle=1)
    test = self.get_test_field('Z_P',cycle=2)
    test = numpy.ma.masked_array(test, mask=vof[:,0]<0.01).compressed()

    error = max(abs(test))
    if (error > tol):
      print 'cycle 2 pressure: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'cycle 2 pressure: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

  def test_void_fraction(self):
    '''Verify void volume fraction at cycle 30'''
    vof = self.get_test_field('VOF',cycle=30)
    test = vof[9,1] # void volume fraction in cell 10 (1-based)
    gold = 3.0e-4
    error = abs(test-gold)
    tol = 1.0e-7
    if (error > tol):
      print 'cycle 30 void fraction: abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'cycle 30 void fraction: abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

  def test_velocity(self):
    '''Verify x-velocity field at cycle 30'''
    tol = 1.0e-4
    test = self.get_test_field('Z_VC',cycle=30)
    gold = 0.9999
    error = max(abs(test[:,0]-gold))
    if (error > tol):
      print 'cycle 30 x-velocity: abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'cycle 30 x-velocity: abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)


if __name__ == '__main__':
  import unittest
  unittest.main()
