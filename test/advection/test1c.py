#!/usr/bin/env python

import sys
import os

import unittest
import numpy
import math

import Truchas
import TruchasTest

class mytest(TruchasTest.GoldenTestCase):

  test_name = 'advection-1c'
  num_procs = 4 # with a parallel executable

  # Override the default setUp, omitting the opening of the golden output
  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())

  def vof_test(self, id, tol):

    # The centroids function does not serialize, so we don't want to here either.
    test = self.test_output.get_simulation().find_series(id).get_data('VOF',serialize=False)
    time = self.test_output.get_simulation().find_series(id).time

    # Analytic vof solution at cell centrioids
    cc = self.test_output.get_mesh().centroids()
    p = -6 + 6*time
    vof = numpy.empty_like(test[:,0])
    for j in range(cc.shape[0]):
      x = (cc[j,0] + cc[j,2])/math.sqrt(2)
      if x < p-0.5:
        vof[j] = 1
      elif x > p+0.5:
        vof[j] = 0
      else:
        vof[j] = 0.5
    
    error = max(abs(test[:,0]-vof))
    if error > tol:
      print 'vof at t=%8.2e: max error = %8.2e: FAIL (tol=%8.2e)'%(time,error,tol)
      self.assertTrue(False)
    else:
      print 'vof at t=%8.2e: max error = %8.2e: PASS (tol=%8.2e)'%(time,error,tol)

  def test_vof1(self):
    '''Verify initial vof field'''
    self.vof_test(1, 1.0e-10)

  def test_vof2(self):
    '''Verify final vof field'''
    self.vof_test(2, 1.0e-10)

if __name__ == '__main__':
  import unittest
  unittest.main()
