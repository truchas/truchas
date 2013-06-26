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

    self.testdirs   = ['ds5_output']
    self.goldendirs = ['ds5_golden', 'ds5_pgolden']

  def setDefinitions(self):
    self.reg = self.testdata[0].getRegion(['Default','all'])

  def testPhi1(self):
    "Verifying the phi1 field"

    tol = 1.0e-10

    self.logger.write("Verifying the phi1 field at early and late times ...\n")
    
    fail = 0
    for n in [ 24, 41 ]:
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
        

  def testPhi23(self):
    "Verifying the phi2 and phi3 fields"

    tol = 1.0e-14

    self.logger.write("Comparing the phi2 and phi3 fields against the reference phi1 ...\n")
    
    fail = 0
    for n in [ 24, 41 ]:
      Cref = self.testdata[0].getField(field='phi1', cycle=n, region=self.reg)
      C = self.testdata[0].getField(field='phi2', cycle=n, region=self.reg)
      error = max(abs((C.data/1.5)-Cref.data)/Cref.data)
      self.logger.write("  Cycle %2d: phi2 max rel error = %10.4e" %(n,error))
      if error > tol:
          self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
          fail = fail + 1
      else:
          self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
      C = self.testdata[0].getField(field='phi3', cycle=n, region=self.reg)
      error = max(abs((C.data/2.0)-Cref.data)/Cref.data)
      self.logger.write("  Cycle %2d: phi3 max rel error = %10.4e" %(n,error))
      if error > tol:
          self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
          fail = fail + 1
      else:
          self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
    if fail: self.fail("  Error exceeds the tolerance")
        
