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
  "DS1 CapabilityTest"

  def setDataStores(self):
    "defines testdata and goldendata output directories needed by this Capability TestCase"

    self.testdirs   = ['two-enclosures_output']
    self.goldendirs = ['two-enclosures_golden']

  def setTolerances(self):
    "define tolerances"

    self.tol['temperature l infinity error'] = 2e-3


  def setDefinitions(self):
    "Define the regions"
    self.reg = self.testdata[0].getRegion(['Default','all'])
    self.golden_reg = self.goldendata[0].getRegion(['Default','all'])

  def testHTencrad(self):
    "Verifying enclosure radiation boundary conditions computed using Chaparral"

    self.logger.write("Verifying Verifying enclosure radiation boundary conditions computed using Chaparral ... ")

    temp = self.testdata[0].getField(field='T',
                                     cycle=50,
                                     region=self.reg)
    
    temp_ref = self.goldendata[0].getField(field='T',
                                           cycle=50,
                   
                        region=self.golden_reg)
    

    TL2err = self.measure.linfError(temp.data,temp_ref.data)
    
    if TL2err > self.tol['temperature l infinity error']:
        self.logger.write(": FAIL (tolerance = %12.7e)\n" %(TL2err))
        self.fail("The temperature does not match the golden output.")
    else:
        self.logger.write(": PASS (tolerance = %12.7e)\n" %(TL2err))
        
