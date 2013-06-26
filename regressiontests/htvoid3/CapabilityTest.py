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
  "HTVOID3 CapabilityTest"

  def setDataStores(self):
    "defines testdata and goldendata output directories needed by this Capability TestCase"

    self.testdirs   = ['htvoid3_output']
    self.goldendirs = ['htvoid3_golden', 'htvoid3_pgolden']

  def setDefinitions(self):
    self.reg = self.testdata[0].getRegion(['Default','all'])

  def testFinalSolidVolumeFraction(self):
    "Verifying the final solid volume fraction"

    tol = 1.0e-6
    n = 276 # should be the final cycle time

    self.logger.write("Verifying the final solid volume fraction...\n")
    
    V = self.testdata[0].getField(field='VOF0002', cycle=n, region=self.reg)
    if 'serial' in self.basicrunspecs.parallel_env:
      Vref = self.goldendata[0].getField(field='VOF0002', cycle=n, region=self.reg)
    else:
      Vref = self.goldendata[1].getField(field='VOF0002', cycle=n, region=self.reg)
    error = max(abs(V.data-Vref.data))
    self.logger.write("  Cycle %2d: maximum abs error = %10.4e\n" %(n,error))
    if error > tol:
        self.logger.write(": FAIL (tolerance = %12.7e)\n" %(tol))
        self.fail("  Error exceeds the tolerance")
    else:
        self.logger.write(": PASS (tolerance = %12.7e)\n" %(tol))
        

  def testFinalTemperature(self):
    "Verifying the final temperature field"

    tol = 1.0e-6
    n = 276 # should be the final cycle time

    self.logger.write("Verifying the final temperature field...\n")
    
    T = self.testdata[0].getField(field='T', cycle=n, region=self.reg)
    if 'serial' in self.basicrunspecs.parallel_env:
      Tref = self.goldendata[0].getField(field='T', cycle=n, region=self.reg)
    else:
      Tref = self.goldendata[1].getField(field='T', cycle=n, region=self.reg)
    error = max(abs(T.data-Tref.data))
    self.logger.write("  Cycle %2d: maximum abs error = %10.4e\n" %(n,error))
    if error > tol:
        self.logger.write(": FAIL (tolerance = %12.7e)\n" %(tol))
        self.fail("  Error exceeds the tolerance")
    else:
        self.logger.write(": PASS (tolerance = %12.7e)\n" %(tol))
        
