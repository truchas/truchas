#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class EM1(TruchasTest.GoldenTestCase):

  test_name='em1'

  num_procs=1

  freespace_id = 10
  freespace_region = None

  conducting_ids = [11, 12]
  conducting_region = None

  # Useful routines to grab the field of interest
  # Python unittest skips these function because they do not
  # match test* or runTest patterns.
  def get_joule_test_data(self,cyc,region=None):
    return self.test_output.get_simulation().find_series(cycle=cyc).get_data('Joule_P',region=region)
  
  def get_joule_gold_data(self,cyc,region=None):
    return self.gold_output.get_simulation().find_series(cycle=cyc).get_data('Joule_P',region=region)

  # Run before each test, must call the base class setUp
  # to ensure all the truchas binary instance is set up correctly
  def setUp(self):
    TruchasTest.GoldenTestCase.setUp(self)
    self.test_simulation=self.test_output.get_simulation()
    self.freespace_region = Truchas.TruchasRegion(self.test_simulation,self.freespace_id)
    self.conducting_region = Truchas.TruchasRegion(self.test_simulation,self.conducting_ids)


  def test_initial_joule_heat(self):
    '''Test the initial joule heat'''
   
    tol = 1.0e-8
    fails=0

    # For cycle = 0
    q=self.get_joule_test_data(0,self.conducting_region)
    q_gold=self.get_joule_gold_data(0,self.conducting_region)
    error=max(abs(q-q_gold)/q_gold)
    if error > tol:
      fails+=1

    # For cycles 1-5
    q_ref=q_gold
    for cyc in range(1,6):
      q_next=self.get_joule_test_data(cyc,self.conducting_region)
      error=max(abs(q_next-q_ref)/q_ref)
      if error > tol:
        fails+=1

    self.assertTrue(fails == 0)	

  def test_scaled_joule_heat(self):
    '''Test scaled joule heat'''

    tol = 1.0e-15
    fails=0

    # For cycle = 6
    q_0=self.get_joule_test_data(0,self.conducting_region)
    q_ref=4*q_0
    q=self.get_joule_test_data(6,self.conducting_region)
    error=max(abs(q-q_ref)/q_ref)
    if error > tol:
      fails+=1

    # For cycles 7-10
    for cyc in range(7,11):
      qnext=self.get_joule_test_data(cyc,self.conducting_region)
      error=max(abs(qnext-q))
      if error > tol:
	fails+=1

    self.assertTrue(fails == 0)

  def test_zero_joule_heat(self):
    '''Test the zero joule heat'''

    fails=0

    # For cycles 11-15
    for cyc in range(11,16):
      q=self.get_joule_test_data(cyc)
      error = max(abs(q))
      if error != 0.0:
	fails+=1

    self.assertTrue(fails == 0)	
      

  def test_final_joule_heat(self):
    '''Test the final joule heat'''

    tol = 1.0e-8
    fails=0

    # For cycle 16
    q=self.get_joule_test_data(16,self.conducting_region)
    q_gold=self.get_joule_gold_data(16,self.conducting_region)
    error=max(abs(q-q_gold)/q_gold)
    if error > tol:
      fails+=1

    # For cycles 17-20
    q_ref=q
    for cyc in range(17,21):
      qnext=self.get_joule_test_data(cyc,self.conducting_region)
      error=max(abs(qnext-q_ref))
      if error != 0.0:
	fails+=1

    self.assertTrue(fails == 0)	

  def test_freespace_joule_heat(self):
    '''Test the freespace joule heat'''
    fails=0
    for cyc in range(0,21):
      q=self.get_joule_test_data(cyc,self.freespace_region)
      q_norm=numpy.linalg.norm(q,numpy.inf)
      if q_norm != 0.0:
	print 'cycle %d q_norm=%1.9e'%(cyc,q_norm)
	fails+=1
    self.assertTrue(fails == 0)

  def test_final_temperature(self):
    '''Test the final temperature'''

    tol = 1.0e-8

    T = self.test_simulation.find_series(cycle=20).get_data('Z_TEMP',self.conducting_region)
    gT = self.gold_output.get_simulation().find_series(cycle=20).get_data('Z_TEMP',self.conducting_region)
    error=max(abs(T-gT)/gT)

    self.assertFalse(error > tol)




if __name__ == '__main__':
  import unittest
  unittest.main()

