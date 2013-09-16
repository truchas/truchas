#!/usr/bin/env python

import sys
import os

import numpy

import Truchas
import TruchasTest

class TMPC(TruchasTest.GoldenTestCase):

  test_name = 'tm-pc'
  num_procs = 4 # with a parallel executable

  def get_test_field(self,field,cycle,region=None):
    return self.test_output.get_simulation().find_series(cycle=cycle).get_data(field,region=region)

  def get_gold_field(self,field,cycle,region=None):
    return self.gold_output.get_simulation().find_series(cycle=cycle).get_data(field,region=region)

  def test_final_temperature(self):
    '''Verify final temperature'''
    n = 216
    tol = 1.0e-6
    gold = self.get_gold_field('Z_TEMP',cycle=n)
    test = self.get_test_field('Z_TEMP',cycle=n)
    error = max(abs((test-gold)/gold))
    if error > tol:
      print 'temperature: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'temperature: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)
      
  def test_final_thermal_strain(self):
    '''Verify final thermal strain with value calculated from temperature'''

    n = 216
    ncells = 18
    cte1 = 2.2e-5
    cte2 = 2.1e-5
    cte3 = 2.0e-5
    T0 = 500.0
    Ti = 400.0
    Th = 375.0
    Tl = 350.0
    error = 0.0
    pcstrain1 = 1.15013007E-03
    pcstrain2 = -1.57897042E-03
    
    T_test = self.get_test_field('Z_TEMP',cycle=n)
    test = self.get_test_field('epstherm',cycle=n)
    gold = self.get_gold_field('epstherm',cycle=n)

#   Here we compare with a calculated solution where the material has not seen
#   any non-isothermal phase change, and a golden solution where it has.
#   The non-isothermal thermal strain is a simple explicit integration of
#   a nonlinear equation.

    print "   T, Calc, reference, error"
    for j in range(ncells):
      T = T_test[j]
      if T >= Ti:
          epsref = (T-T0) * cte1
          diff = abs((test[j,0] - epsref) / epsref)
          print "  > Ti"
      elif T >= Th:
          epsref = pcstrain1 + (T - Ti) * cte2
          diff = abs((test[j,0] - epsref) / epsref)
          print "  > Th"
      elif T >= Tl:
          epsref = gold[j,0]
          diff = abs((test[j,0] - epsref) / epsref)
          print "  > Tl"
      elif T < Tl:
          epsref = pcstrain2 + (T - Tl) * cte3
          diff = abs((test[j,0] - epsref) / epsref)
          print "  < Tl"

      print "  %13.7e %13.7e %13.7e %13.7e"%(T, test[j,0], epsref, diff)
      if diff > error: error = diff

    print "  Cycle %2d: maximum relative thermal strain error = %10.4e"%(n,error)

    self.assertTrue(error <= 1.0e-8)
 
  

  def test_final_stress(self):
    '''Verify the final stress field'''
    
    n = 216 # final cycle number
    abs_tol = 2.0e0
    rel_tol = 1.0e-8
    fail = 0

    epstherm = self.get_test_field('epstherm',cycle=n)
    l1 = 5.20e+10
    l2 = 2.60e+10

    XXref = -epstherm[:,0] * (2.0*l1 + 2.0*l2 - 2.0*l1**2/(l1 + 2.0*l2))
    YYref = XXref

    test = self.get_test_field('sigma',cycle=n)

    error = max(abs((test[:,0]-XXref)/XXref))
    tol = rel_tol
    if error > tol:
      fail += 1
      print 'sigma_xx: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_xx: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs((test[:,1]-YYref)/YYref))
    tol = rel_tol
    if error > tol:
      fail += 1
      print 'sigma_yy: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_yy: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs(test[:,2]))
    tol = abs_tol
    if error > tol:
      fail += 1
      print 'sigma_zz: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_zz: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs(test[:,3]))
    tol = abs_tol
    if error > tol:
      fail += 1
      print 'sigma_xy: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_xy: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs(test[:,4]))
    tol = abs_tol
    if error > tol:
      fail += 1
      print 'sigma_xz: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_xz: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs(test[:,5]))
    tol = abs_tol
    if error > tol:
      fail += 1
      print 'sigma_xx: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_xx: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    self.assertTrue(fail == 0)
  

  def test_final_strain(self):
    '''Verify the final strain field'''
    
    n = 216 # final cycle number
    abs_tol = 1.0e-10
    rel_tol = 1.0e-8
    fail = 0

    epstherm = self.get_test_field('epstherm',cycle=n)
    l1 = 5.20e+10
    l2 = 2.60e+10
    ZZref = epstherm[:,0] * (1.0 + 2.0*l1/(l1 + 2.0*l2))

    test = self.get_test_field('epsilon',cycle=n)

    error = max(abs(test[:,0]))
    tol = abs_tol
    if error > tol:
      fail += 1
      print 'sigma_xx: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_xx: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs(test[:,1]))
    tol = abs_tol
    if error > tol:
      fail += 1
      print 'sigma_yy: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_yy: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs((test[:,2]-ZZref)/ZZref))
    tol = rel_tol
    if error > tol:
      fail += 1
      print 'sigma_zz: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_zz: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs(test[:,3]))
    tol = abs_tol
    if error > tol:
      fail += 1
      print 'sigma_xy: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_xy: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs(test[:,4]))
    tol = abs_tol
    if error > tol:
      fail += 1
      print 'sigma_xz: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_xz: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    error = max(abs(test[:,5]))
    tol = abs_tol
    if error > tol:
      fail += 1
      print 'sigma_xx: max abs error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
    else:
      print 'sigma_xx: max abs error = %8.2e: PASS (tol=%8.2e)'%(error,tol)

    self.assertTrue(fail == 0)


if __name__ == '__main__':
  import unittest
  unittest.main()

