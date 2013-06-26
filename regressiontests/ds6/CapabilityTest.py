#!/usr/bin/env python
if __name__=='__main__':

    import os, sys

    #developers : please specify the location of your truchas checkout directory 'truchasdir'
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    truchasdir  = os.path.abspath(thisdir + '/../../')

    errstring  = "\n\n    truchasdir is set to %s. No tools directory exists in this directory \n" %(truchasdir)
    errstring += "           Please alter your choice of truchasdir in your CapabilityTest script. \n"

    assert 'tools' in os.listdir(truchasdir), errstring

    extension   = '/tools/PythonPackages/TestingInfrastructure/'
    if os.path.isabs(truchasdir):
        testdir = os.path.abspath(truchasdir + extension)
    else:
        testdir = os.path.expanduser(truchasdir + extension)
    sys.path.append(testdir)
    sys.path.append(testdir+'/Runners')

    from RunThisCapability import RunThisCapability

    #developers : can specify clean=1 to remove all outputs and logdirs
    RunThisCapability(testdir,clean=0)

from TestCases import TruchasCapabilityTest

class CapabilityTest(TruchasCapabilityTest):
  "DS2 CapabilityTest"

  def setDataStores(self):
    "defines testdata and goldendata output directories needed by this Capability TestCase"

    self.testdirs   = ['ds6_output']
    self.goldendirs = ['ds6_golden', 'ds6_pgolden']

  def setDefinitions(self):
    self.reg = self.testdata[0].getRegion(['Default','all'])

  def testTemp(self):
    "Verifying the temperature field"

    tol = 1.0e-10

    self.logger.write("Verifying the temperature field at early and late times ...\n")
    
    fail = 0
    for n in [ 22, 50 ]:
      T = self.testdata[0].getField(field='T', cycle=n, region=self.reg)
      if 'serial' in self.basicrunspecs.parallel_env:
        Tref = self.goldendata[0].getField(field='T', cycle=n, region=self.reg)
      else:
        Tref = self.goldendata[1].getField(field='T', cycle=n, region=self.reg)
      error = max(abs(T.data-Tref.data)/Tref.data)
      self.logger.write("  Cycle %2d: max rel error = %10.4e" %(n,error))
      if error > tol:
          self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
          fail = fail + 1
      else:
          self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
    if fail: self.fail("  Error exceeds the tolerance")
        

  def testPhi1(self):
    "Verifying the phi1 field"

    tol = 1.0e-10

    self.logger.write("Verifying the phi1 field at early and late times ...\n")
    
    fail = 0
    for n in [ 22, 50 ]:
      C = self.testdata[0].getField(field='phi1', cycle=n, region=self.reg)
      if 'serial' in self.basicrunspecs.parallel_env:
        Cref = self.goldendata[0].getField(field='phi1', cycle=n, region=self.reg)
      else:
        Cref = self.goldendata[1].getField(field='phi1', cycle=n, region=self.reg)
      error = max(abs(C.data-Cref.data)/Cref.data)
      self.logger.write("  Cycle %2d: max rel error = %10.4e" %(n,error))
      if error > tol:
          self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
          fail = fail + 1
      else:
          self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
    if fail: self.fail("  Error exceeds the tolerance")
        

  def testPhi2(self):
    "Verifying the phi2 field"

    tol = 1.0e-10

    self.logger.write("Verifying the phi2 field at early and late times ...\n")
    
    fail = 0
    for n in [ 22, 50 ]:
      C = self.testdata[0].getField(field='phi2', cycle=n, region=self.reg)
      if 'serial' in self.basicrunspecs.parallel_env:
        Cref = self.goldendata[0].getField(field='phi2', cycle=n, region=self.reg)
      else:
        Cref = self.goldendata[1].getField(field='phi2', cycle=n, region=self.reg)
      error = max(abs(C.data-Cref.data)/Cref.data)
      self.logger.write("  Cycle %2d: max rel error = %10.4e" %(n,error))
      if error > tol:
          self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
          fail = fail + 1
      else:
          self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
    if fail: self.fail("  Error exceeds the tolerance")
        
