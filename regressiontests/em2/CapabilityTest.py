#!/usr/bin/env python
if __name__=='__main__':

    import os, sys

    #developers : please specify the location of your truchas checkout directory 'truchasdir'
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    truchasdir  = os.path.abspath(thisdir + '/../..')

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
    "EM2 CapabilityTest" 

    def setDataStores(self):
        "defines testdata and goldendata output directories needed by this Capability TestCase"

        self.testdirs   = ['em2_output']
        self.goldendirs = ['em2_golden']

    def testJouleHeat(self):
      "Verifying the Joule heat fields"
      
      self.logger.write("Verifying the Joule heat fields...\n")
      
      tol = 1.0e-8  # Tolerance for the max relative error norm.
      
      # Three Joule heat calculations, used for the following cycle groups.
      groups = [ range(0,5), range(5,16), range(16,21) ]
      
      fail = 0
      for g in groups:
        n = g[0]
        q = self.testdata[0].getField(field='Joule_P', cycle=n)
        qref = self.goldendata[0].getField(field='Joule_P', cycle=n)
        error = max(abs(q.data-qref.data)/qref.data)
        self.logger.write("  Cycle %2d: maximum relative error = %10.4e" %(n,error))
        if error > tol:
          fail = fail + 1
          self.logger.write(": FAIL (tolerance = %3.2e)\n" %(tol))
        else:
          self.logger.write(": PASS (tolerance = %8.2e)\n" %(tol))
        for j in g[1:]:
          qnext = self.testdata[0].getField(field='Joule_P', cycle=j)
          error = max(abs(qnext.data-q.data))
          if error != 0:
            fail = fail + 1
            self.logger.write("  Unexpected change in the Joule heat for cycle %d: FAIL\n" %(j))

      if fail: self.fail("  Incorrect Joule heat fields: %d failed tests" %(fail))
        
    def testFinalTemperature(self):
      "Verifying the final temperature field..."

      tol = 1.0e-8  # Tolerance for the max relative error norm.
      
      self.logger.write("Verifying the final temperature field...\n")
      t = self.testdata[0].getField(field='T', cycle=20)
      tref = self.goldendata[0].getField(field='T', cycle=20)
      error = max(abs(t.data-tref.data)/tref.data)
      self.logger.write("  Maximum relative error = %10.4e" %(error))
      if error > tol:
        self.logger.write(": FAIL (tolerance = %3.2e)\n" %(tol))
        self.fail("  Error exceeds the tolerance %10.4e" %(tol))
      else:
        self.logger.write(": PASS (tolerance = %8.2e)\n" %(tol))
        
