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
  "DS8 CapabilityTest"

  def setDataStores(self):
    "defines testdata and goldendata output directories needed by this Capability TestCase"

    self.testdirs   = ['ds8_output']
    self.goldendirs = ['ds8_golden', 'ds8_pgolden']

  def setDefinitions(self):
    self.reg = self.testdata[0].getRegion(['Default','all'])

  def testTemp(self):
    "Verifying the temperature field"

    tol = 1.0e-10

    self.logger.write("Verifying the temperature field at early and late times ...\n")
    
    fail = 0
    for n in [ 94, 290 ]:
      T = self.testdata[0].getField(field='T', cycle=n, region=self.reg)
      if 'serial' in self.basicrunspecs.parallel_env:
        Tref = self.goldendata[0].getField(field='T', cycle=n, region=self.reg)
      else:
        Tref = self.goldendata[1].getField(field='T', cycle=n, region=self.reg)
      error = max(abs(T.data-Tref.data)/Tref.data)
      self.logger.write("  Cycle %2d: maximum relative error = %10.4e" %(n,error))
      if error > tol:
          self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
          fail = fail + 1
      else:
          self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
    if fail: self.fail("  Error exceeds the tolerance")
        
  def testVelocity(self):
    "Verifying the velocity field"

    tol = 1.0e-10

    self.logger.write("Verifying the velocity field at early and late times ...\n")
    
    fail = 0
    for n in [ 94, 290 ]:
      V = self.testdata[0].getField(field='Velocity', cycle=n, region=self.reg)
      if 'serial' in self.basicrunspecs.parallel_env:
        Vref = self.goldendata[0].getField(field='Velocity', cycle=n, region=self.reg)
      else:
        Vref = self.goldendata[1].getField(field='Velocity', cycle=n, region=self.reg)
      uerror = max(abs(V.data[:,0]-Vref.data[:,0]))
      verror = max(abs(V.data[:,1]-Vref.data[:,1]))
      werror = max(abs(V.data[:,2]-Vref.data[:,2]))
      error = max(uerror,verror)
      self.logger.write("  Cycle %2d: maximum error = %10.4e" %(n,error))
      if error > tol:
          self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
          fail = fail + 1
      else:
          self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
    if fail: self.fail("  Error exceeds the tolerance")
