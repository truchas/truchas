#!/usr/bin/env python
if __name__=='__main__':

    import os, sys

    #developers : please specify the location of your truchas checkout directory 'truchasdir'
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    truchasdir  = os.path.abspath(thisdir + '../../../../../../')

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
  "EM3 CapabilityTest" 

  def setDataStores(self):
    "defines testdata and goldendata output directories needed by this Capability TestCase"

    self.testdirs     = ['em3_output']
    self.goldendirs   = ['em3_golden']

  def testInitialJouleHeat(self):
    "Verifying the initial Joule heat field"

    self.logger.write("Verifying the initial Joule heat field...\n")
    q = self.testdata[0].getField(field='Joule_P', cycle=0)
    qref = self.goldendata[0].getField(field='Joule_P', cycle=0)
    error = max(abs(q.data-qref.data)/qref.data)
    self.logger.write("  Maximum relative error = %10.4e" %(error))
    tol = 1.0e-10
    if error > tol:
      self.logger.write(": FAIL (tolerance = %3.2e)\n" %(tol))
      self.fail("Incorrect initial Joule heat field")
    else:
      self.logger.write(": PASS (tolerance = %8.2e)\n" %(tol))

  def testFinalTemperature(self):
    "Verifying the final temperature field"

    self.logger.write("Verifying the final temperature field...\n")
    t = self.testdata[0].getField(field='T', cycle=20)
    tref = self.goldendata[0].getField(field='T', cycle=20)
    error = max(abs(t.data-tref.data)/tref.data)
    self.logger.write("  Maximum relative error = %10.4e" %(error))
    tol = 1.0e-10
    if error > tol:
      self.logger.write(": FAIL (tolerance = %3.2e)\n" %(tol))
      self.fail("Incorrect final temperature field" %(tol))
    else:
      self.logger.write(": PASS (tolerance = %8.2e)\n" %(tol))
        
