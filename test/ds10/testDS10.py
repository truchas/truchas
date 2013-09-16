#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class DS10(TruchasTest.GoldenTestCase):

  test_name = 'ds10'
  num_procs = 4 # with a parallel executable
  
  # Override the default setUp, omitting the opening of the golden output
  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass() # This runs Truchas
    self.test_output = Truchas.TruchasOutput(self.get_output_file())

  def get_test_field(self,field,cycle,serialize=True,region=None):
    return self.test_output.get_simulation().find_series(cycle=cycle).get_data(field,serialize,region=region)

  def test_final_temp(self):
    '''Verify the final temperature field'''
    
    tol = 1.0e-7
    
    coord = self.test_output.get_mesh().centroids()
    xcoord = coord[:,0]
    ycoord = coord[:,1]
    
    Tref = xcoord + ycoord
    for j in range(Tref.shape[0]):
      x = xcoord[j]
      y = ycoord[j]
      if x < 0:
        if y < 0:
          Tref[j] += 2.5
        else:
          Tref[j] += 2.7
      else:
        if y < 0:
          Tref[j] += 3.3
        else:
          Tref[j] += 3.5
    
    # The centroids function does not serialize the data, so we don't want to either.
    T = self.get_test_field('Z_TEMP',cycle=233,serialize=False)
    
    # Take care of gap cells so they don't enter into error calc
    for j in range(Tref.shape[0]):
      x = xcoord[j]
      y = ycoord[j]
      if abs(x) < 1.0e-3 or abs(y) < 1.0e-3:
        Tref[j] = T[j]
        
    error = max(abs(T-Tref)/Tref)
    if error > tol:
      print 'temperature max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'temperature max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

if __name__ == '__main__':
  import unittest
  unittest.main()

