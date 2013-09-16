#!/usr/bin/env python

import sys
import os

import unittest
import numpy

import Truchas
import TruchasTest

class ChannelFlow(TruchasTest.GoldenTestCase):

  test_name = 'channel-flow'
  num_procs = 4 # with a parallel executable
  restart_file = 'restart.bin'

  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())

  def AnalyticSolution(self,Ncellx,Ncelly):

    # Analytic solution to balance between viscous stress and pressure gradient
    #      Assumes equal size cells along each axis
    Pinlet = 1000.
    Poutlet = 0.
    Length = 1.
    GradP = (Poutlet-Pinlet)/Length
    Rho = 1.
    Mu  = 1.
    Height   = 1.
    DeltaY = Height/Ncelly
    RHS = DeltaY*DeltaY*GradP/Mu
    Beta=-3.
    BCoeff=-2.
    Ux=numpy.array([RHS/Beta],'d')
    Gamma=numpy.array([],'d')

    #Forward Elimination Loop
    for j in range(1,Ncelly):
      if j==Ncelly-1:
        BCoeff=-3.
      Gamma.resize([j+1])
      Gamma[j]=1./Beta
      Beta=BCoeff-Gamma[j]
      Ux.resize([j+1])
      Ux[j]=(RHS-Ux[j-1])/Beta

    #Back Substitution Loop
    for i in range(1,Ncelly):
      j=Ncelly-i-1
      Ux[j]=Ux[j]-Gamma[j+1]*Ux[j+1]
    return Ux


  def get_test_field(self,field,cycle,region=None):
    return self.test_output.get_simulation().find_series(cycle=cycle).get_data(field,region=region)


  def test_velocity(self):
    '''Verify x-velocity field'''

    Ncellx = 2
    Ncelly = 7
    n = 51025
    tol = 1.0e-12

    test = self.get_test_field('Z_VC',cycle=n)
    test = test[::Ncellx,0] # regular ncellx by ncelly grid -- pulling x-vel from first column

    gold = self.AnalyticSolution(Ncellx,Ncelly)

    error = max(abs((test-gold)/gold))
    if error > tol:
      print 'x-velocity: max rel error = %8.2e: FAIL (tol=%8.2e)'%(error,tol)
      self.assertTrue(False)
    else:
      print 'x-velocity: max rel error = %8.2e: PASS (tol=%8.2e)'%(error,tol)


if __name__ == '__main__':
  import unittest
  unittest.main()
