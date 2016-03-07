#!/usr/bin/env python

import sys
import os

import numpy
import numpy.ma
import unittest

import Truchas
import TruchasTest

class UStrucTest(TruchasTest.GoldenTestCase):

  test_name = 'ustruc2'
  num_procs = 4 # with a parallel executable
  restart_file = 'ustruc1.restart.149'

  def get_test_field(self,field,id,region=None):
    return self.test_output.get_simulation().find_series(id=id).get_data(field,region=region)

  def get_gold_field(self,field,id,region=None):
    return self.gold_output.get_simulation().find_series(id=id).get_data(field,region=region)

  def test_final_temp(self):
    '''USTRUC2: verifying the temperature field at the final time'''
    tol = 5.0e-5
    test = self.get_test_field('Z_TEMP',id=2)
    gold = self.get_gold_field('Z_TEMP',id=3)
    error = max(abs((test-gold)/gold))
    print 'final temp max rel error =', error, '(tol =', tol, ')'
    self.assertTrue(error <= tol)

  def test_final_G(self):
    '''USTRUC2: verifying the G field at the final time'''
    tol = 1.0e-4
    test = self.get_test_field('uStruc-G',id=2)
    gold = self.get_gold_field('uStruc-G',id=3)
    test = numpy.ma.masked_where(gold == 0, test)
    gold = numpy.ma.masked_where(gold == 0, gold)
    error = numpy.ma.max(abs((test-gold)/(200/tol + gold)))
    print 'final G max rel error =', error, '(tol =', tol, ')'
    self.assertTrue(error <= tol)

  def test_final_V(self):
    '''USTRUC2: verifying the V field at the final time'''
    tol = 1.0e-4
    test = self.get_test_field('uStruc-V',id=2)
    gold = self.get_gold_field('uStruc-V',id=3)
    test = numpy.ma.masked_where(gold == 0, test)
    gold = numpy.ma.masked_where(gold == 0, gold)
    error = numpy.ma.max(abs((test-gold)/(1.0e-4/tol + gold)))
    print 'final V max rel error =', error, '(tol =', tol, ')'
    self.assertTrue(error <= tol)

  def test_final_L1(self):
    '''USTRUC2: verifying the lambda1 field at the final time'''
    tol = 1.0e-4
    test = self.get_test_field('uStruc-gv1-lambda1',id=2)
    gold = self.get_gold_field('uStruc-gv1-lambda1',id=3)
    test = numpy.ma.masked_where(gold == 0, test)
    gold = numpy.ma.masked_where(gold == 0, gold)
    error = numpy.ma.max(abs((test-gold)/(4.0e-6/tol + gold)))
    print 'final lambda1 max rel error =', error, '(tol =', tol, ')'
    self.assertTrue(error <= tol)

  def test_final_L2(self):
    '''USTRUC2: verifying the lambda2 field at the final time'''
    tol = 1.0e-4
    test = self.get_test_field('uStruc-gv1-lambda2',id=2)
    gold = self.get_gold_field('uStruc-gv1-lambda2',id=3)
    test = numpy.ma.masked_where(gold == 0, test)
    gold = numpy.ma.masked_where(gold == 0, gold)
    error = numpy.ma.max(abs((test-gold)/(2.0e-6/tol + gold)))
    print 'final lambda2 max rel error =', error, '(tol =', tol, ')'
    self.assertTrue(error <= tol)

  def test_final_ustruc(self):
    '''USTRUC2: verifying the microstructure map at the final time'''
    test = self.get_test_field('uStruc-gv1-ustruc',id=2)
    gold = self.get_gold_field('uStruc-gv1-ustruc',id=3)
    self.assertTrue(numpy.all(test == gold))

if __name__ == '__main__':
  import unittest
  unittest.main()
