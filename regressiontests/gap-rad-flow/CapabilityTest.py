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
  "gap-rad-flow CapabilityTest"

  def setDataStores(self):
    "defines testdata and goldendata output directories needed by this Capability TestCase"

    self.testdirs   = ['gap-rad-flow_output']

  def setDefinitions(self):
    self.probe_le = self.testdata[0].getProbe(name='left end',field='TEMP')
    self.probe_gl = self.testdata[0].getProbe(name='gap left',field='TEMP')
    self.probe_re = self.testdata[0].getProbe(name='right end',field='TEMP')
    self.probe_gr = self.testdata[0].getProbe(name='gap right',field='TEMP')
    
  def testInitialProbeValues(self):
    "Verifying the initial probe values"
    
    fail = 0
    tol = 1.e-9 # use tight a tight tolerance for initial conditions
    self.logger.write('Checking initial probe values...\n')
    
    error = abs(self.probe_le.data[0,2] - 1.499998474)
    if error > tol:
        fail += 1
        self.logger.write('  Left-end probe error %8.2e exceeds tolerance %8.2e\n' %(error,tol))
    
    error = abs(self.probe_re.data[0,2] - 0.500001526)
    if error > tol:
        fail += 1
        self.logger.write('  Right-end probe error %8.2e exceeds tolerance %8.2e\n' %(error,tol))
    
    error = abs(self.probe_gl.data[0,2] - 1.498533630)
    if error > tol:
        fail += 1
        self.logger.write('  Gap-left probe error %8.2e exceeds tolerance %8.2e\n' %(error,tol))
    
    error = abs(self.probe_gr.data[0,2] - 0.5014663700)
    if error > tol:
        fail += 1
        self.logger.write('  Gap-right probe error %8.2e exceeds tolerance %8.2e\n' %(error,tol))
    
    if fail:
        self.fail('  FAIL: initial probe error exceeds tolerance\n')
    
  def testFinalProbeValues(self):
    "Verifying the final probe values"
    
    # Reference values computed from a 2-term asympotic expansion of the
    # exact solution in powers of 1/k, where k=100 is the thermal conductivity.
    
    fail = 0
    tol = 5.e-5
    n = self.probe_le.data.size/3 - 1
    self.logger.write('Checking final probe values...\n')
    
    error = abs(self.probe_le.data[n,2] - 1.371828885)
    if error > tol:
        fail += 1
        self.logger.write('  Left-end probe error %8.2e exceeds tolerance %8.2e\n' %(error,tol))
    
    error = abs(self.probe_re.data[n,2] - 0.628171115)
    if error > tol:
        fail += 1
        self.logger.write('  Right-end probe error %8.2e exceeds tolerance %8.2e\n' %(error,tol))
    
    error = abs(self.probe_gl.data[n,2] - 1.370837833)
    if error > tol:
        fail += 1
        self.logger.write('  Gap-left probe error %8.2e exceeds tolerance %8.2e\n' %(error,tol))
    
    error = abs(self.probe_gr.data[n,2] - 0.629162167)
    if error > tol:
        fail += 1
        self.logger.write('  Gap-right probe error %8.2e exceeds tolerance %8.2e\n' %(error,tol))
    
    if fail:
        self.fail('  FAIL: final probe error exceeds tolerance\n')
