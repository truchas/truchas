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
  "Grain Growth CapabilityTest"

  def setTolerances(self):
    "Defines error tolerances"

    self.tol['R_error'] = 1.0e-02
    
  def setDataStores(self):
    "Defines testdata output directories needed by this Capability TestCase"

    self.testdirs   = ['grain_growth_output']

  def setDefinitions(self):
    "Define the radii at cell 1 (for t=1.5), cell 2 (for t=2.5) and cell 3 (for t=3.0)"

    self.radiusc1 = 6.54e-04
    self.radiusc2 = 3.79e-04
    self.radiusc3 = 8.14e-04

    self.c1       = self.testdata[0].getRegion(selector = ['ID','CELL','1'],
                                               desc     = 'Cell 1')
    self.c2       = self.testdata[0].getRegion(selector = ['ID','CELL','2'],
                                               desc     = 'Cell 2')        
    self.c3       = self.testdata[0].getRegion(selector = ['ID','CELL','3'],
                                               desc     = 'Cell 3')
    
  def testGrainRadius(self):
    "verifying grain growth radius"

    #get list of timesteps from Truchas testdata output that satisfy 1.5<=t<=3.5 
    these_timesteps = self.testdata[0].getTimeSteps(timerange = [1.5,3.5],
                                                    dtstep    = 0.5)

    for timestep in these_timesteps:

        thistime   = timestep.time
        index      = these_timesteps.index(timestep)

        if index == 0:
            rc1         = self.testdata[0].getField(field    = 'Grain_R',
                                                    time     = thistime,
                                                    region   = self.c1)
            
            err         = self.measure.percentError(rc1.data[0],self.radiusc1)
            #fail if err > R_error
            failmessage = 'In %s: R1 percent error = %10.3e is > %10.3e' %(self.methodName, err, self.tol['R_error'])
	    self.failIf(err>self.tol['R_error'],msg = failmessage)

        if index == 2:
            rc2         = self.testdata[0].getField(field    = 'Grain_R',
                                                    time     = thistime,
                                                    region   = self.c2)
            
            err         = self.measure.percentError(rc2.data[0],self.radiusc2)
            
            #fail if err > R_error
            failmessage = 'In %s: R2 percent error = %10.3e is > %10.3e' %(self.methodName, err, self.tol['R_error'])
	    self.failIf(err>self.tol['R_error'],msg = failmessage)
            

        if index == 3:
            rc3         = self.testdata[0].getField(field    = 'Grain_R',
                                                    time     = thistime,
                                                    region   = self.c3)

            err         = self.measure.percentError(rc3.data[0],self.radiusc3)
            #fail if err > R_error
            failmessage = 'In %s: R2 percent error = %10.3e is > %10.3e' %(self.methodName, err, self.tol['R_error'])
	    self.failIf(err>self.tol['R_error'],msg = failmessage)
