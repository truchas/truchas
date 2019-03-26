#!/usr/bin/env python

import sys
import os

import unittest
import numpy

import Truchas
import TruchasTest

class NaturalConvTet(TruchasTest.GoldenTestCase):

  test_name = 'nat-conv-tet-old'
  num_procs = 4 # with a parallel executable
  restart_file = 'restart-tet.bin'

  def test_final_velocity(self):
    '''NaturalConvTet: verifying the final velocity field'''
    tol = 8.0e-5  # seeing an error of approx 3.7e-5
    test = self.test_output.get_simulation().find_series(id=2).get_data('Z_VC')
    gold = self.gold_output.get_simulation().find_series(id=3).get_data('Z_VC')
    error = max(abs((test[:,0]-gold[:,0]))) / max(abs(gold[:,0]))
    # Z velocity error
    fail = 0
    if error > tol:
      print 'max(|u - u_ref|) / max(u_ref) = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      fail = 1
    else:
      print 'max(|u - u_ref|) / max(u_ref) = %8.2e: PASS (tol=%8.2e)'%(error,tol)
    # Z velocity error
    error = max(abs((test[:,2]-gold[:,2]))) / max(abs(gold[:,2]))
    if error > tol:
      print 'max(|w - w_ref|) / max(w_ref) = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      fail = 1
    else:
      print 'max(|w - w_ref|) / max(w_ref) = %8.2e: PASS (tol=%8.2e)'%(error,tol)
    if fail:
      self.assertTrue(False)
    
if __name__ == '__main__':
  import unittest
  unittest.main()
