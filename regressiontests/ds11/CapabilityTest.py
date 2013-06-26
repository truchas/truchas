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

    self.testdirs   = ['ds11_output']
    self.goldendirs = ['ds11_golden', 'ds11_pgolden']

  def setDefinitions(self):
    self.reg = self.testdata[0].getRegion(['Default','all'])

  def testTemperature(self):
    "Verifying the temperature field"

    tol = 1.0e-9
    
    self.logger.write("Verifying the temperature field...\n")
    
    fail = 0
    for n in [45, 109]:
      T = self.testdata[0].getField(field='T', cycle=n, region=self.reg)
      if 'serial' in self.basicrunspecs.parallel_env:
        Tref = self.goldendata[0].getField(field='T', cycle=n, region=self.reg)
      else:
        Tref = self.goldendata[1].getField(field='T', cycle=n, region=self.reg)
      error = max(abs(T.data-Tref.data)/Tref.data)
      self.logger.write("  Cycle %2d: maximum relative error = %10.4e\n" %(n,error))
      if error > tol:
        self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
        fail += 1
      else:
        self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
    if fail: self.fail("  Error exceeds the tolerance")
        
  def testFluidVolumeFraction(self):
    "Verifying the fluid volume fraction"

    tol = 1.0e-8

    self.logger.write("Verifying the fluid volume fraction...\n")
    
    fail = 0
    for n in [45, 109]:
      u = self.testdata[0].getField(field='VOF0003', cycle=n, region=self.reg)
      if 'serial' in self.basicrunspecs.parallel_env:
        uref = self.goldendata[0].getField(field='VOF0003', cycle=n, region=self.reg)
      else:
        uref = self.goldendata[1].getField(field='VOF0003', cycle=n, region=self.reg)
      error = max(abs(u.data-uref.data))
      self.logger.write("  Cycle %2d: maximum error = %10.4e\n" %(n,error))
      if error > tol:
        self.logger.write(": FAIL (tolerance = %12.7e)\n" %(tol))
        fail += 1
      else:
        self.logger.write(": PASS (tolerance = %12.7e)\n" %(tol))
    if fail: self.fail("  Error exceeds the tolerance")
        
