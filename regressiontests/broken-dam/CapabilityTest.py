#!/usr/bin/env python
if __name__=='__main__':

    import os, sys
    #specify the location of your truchas checkout directory 'truchasdir'
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    truchasdir  = os.path.abspath(thisdir + '../../../')

    errstring  = "\n\n    truchasdir is set to %s. \
                  No tools directory exists in this directory \n" %(truchasdir)
    errstring += \
       "           Please alter truchasdir in your CapabilityTest script. \n"
    
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

    def setDataStores(self):
        "defines all datastores to be used throughout this TestCase"

        self.testdirs   = ['broken-dam_output']
        self.goldendirs = ['broken-dam-golden']

    def setDefinitions(self):
        self.reg = self.testdata[0].getRegion(['Default','all'])
        
    def testVelocity(self):
        "Test velocity"
        
        self.logger.write("Verifying the final velocity field ...\n")
        
        tol = 1.0e-9
        n = 61
        
        V = self.testdata[0].getField(field='Velocity', cycle=n, region=self.reg)
        Vref = self.goldendata[0].getField(field='Velocity', cycle=n, region=self.reg)
        uerror = max(abs(V.data[:,0]-Vref.data[:,0]))
        verror = max(abs(V.data[:,1]-Vref.data[:,1]))
        #werror = max(abs(V.data[:,2]-Vref.data[:,2]))/max(abs(Vref.data[:,2]))
        error = max(uerror,verror)
        self.logger.write("  Cycle %3d: maximum error = %10.4e" %(n,error))
        if (error > tol):
          self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
          self.fail( "  Error exceeds the tolerance")
        else:
          self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
          
    def testVof(self):
        "test fluid VoF"
        
        self.logger.write("verifying the final water volume fractions ...\n")
        
        tol = 1.0e-9
        n = 61
        
        V = self.testdata[0].getField(field='VOF0001', cycle=n, region=self.reg)
        Vref = self.goldendata[0].getField(field='VOF0001', cycle=n, region=self.reg)
        error = max(abs(V.data-Vref.data))
        self.logger.write("  Cycle %3d: maximum error = %10.4e" %(n,error))
        if (error > tol):
          self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
          self.fail( "  Error exceeds the tolerance")
        else:
          self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
