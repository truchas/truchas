# This test was modified to avoid having the fluid run out of volume at the end
# The velocity was reduced from 1.0 to 0.9999 and out put adjusted accordingly

try:
    import numpy.oldnumeric as Numeric
except ImportError:
    import Numeric
except:
    raise

if __name__=='__main__':

    import os, sys
    #developers : please specify the location of your truchas checkout directory 'truchasdir'
    thisdir     = os.path.abspath(os.path.dirname(__file__)) 
    truchasdir  = os.path.abspath(thisdir + '../../../')

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

    def setTolerances(self):
        "define tolerances here"
        
        self.tol['vof2_error'] = 1.0e-07
        self.tol['VL1_error']  = 1.0e-04
        self.tol['P_error']    = 1.0e-03

    def setDataStores(self):
        "defines all datastores to be used throughout this TestCase"

        #define testdata output directories..these result from running a BasicRun TestCase
        self.testdirs     = ['void_collapse_output']

    def setDefinitions(self):
        "defines analytic expressions and regions needed by this Capability TestCase"

        #get region representing the entire mesh
        self.meshregion  = self.testdata[0].getRegion(selector=['Default'],
                                                      desc    = 'Region representing the entire Truchas mesh')

        self.fluid       = self.testdata[0].getRegion(selector= ['VAR','VOF0001',[0.01,1.0]],
                                                      desc    = 'Fluid region')

        #self.fluid.str()

#  This test seems to check some numerical discretization transient that I do not understand. DAK 8-24-10
    def testPressures1(self):
        "test Cycle 1 pressures"

        P1  = self.testdata[0].getField(field  = 'P',
                                       region = self.fluid,
                                       cycle  = 1)

        #P1.str()
        
#        P_anal = [66,56,46,36,26,16,6,0.5]
        P_anal = [65.9934,55.9944,45.9954,35.9964,25.9974,15.9984,5.9994,0.49995]
        P_anal = Numeric.array(P_anal)

        #P_anal.str()
        
	#get Linf error for P1
	Linferr    = self.measure.linfError(P1.data,P_anal)

        #fail if Linferr > P_error
        failmessage = 'In %s pressure Linf_error = %10.3e is > %10.3e' %(self.methodName, Linferr, self.tol['P_error'])

	self.failIf(Linferr > self.tol['P_error'],
                    msg = failmessage)
        
    def testPressures2(self):
        "test Cycle 2 pressures"

        P2  = self.testdata[0].getField(field  = 'P',
                                        region = self.fluid,
                                        cycle  = 2)

        #P2.str()
      
        
	#get Linf error for P2
	Linferr    = self.measure.linfError(P2.data)

        #fail if Linferr > P_error
        failmessage = 'In %s pressure Linf_error = %10.3e is > %10.3e' %(self.methodName, Linferr, self.tol['P_error'])

	self.failIf(Linferr > self.tol['P_error'],
                    msg = failmessage)        
        

    def testVof2(self):
	"tests Vof2 in final cell at cycle 30 "

        #get Vof2 field 
	Vof2  = self.testdata[0].getField(field  = 'VOF0002',
                                          region = self.meshregion,
                                          cycle  = 30)
        #Vof2.str()
        
	#log Vof2 field
	self.logger.writeField(Vof2,
                               region = self.meshregion,
                               desc   = 'Vof2 in the entire domain')

        # define error in residual void in cell 10
        err = Vof2.data[9] - 3.0e-4
        
        logmessage   = '\n Cell 10, cycle 30 Void Fraction = %10.3e \n' %(err)
        self.logger.write(logmessage)
	
	#fail if err > vof2_error
        failmessage = 'In %s: Residual Void = %10.3e is > %10.3e' %(
                       self.methodName, err, self.tol['vof2_error'])

	self.failIf(err > self.tol['vof2_error'],
                    msg = failmessage)


    def testXVelocities(self):
	"tests x direction velocity generation in entire domain"
        
 	#get velocities on the entire domain
	V     = self.testdata[0].getField(field  = 'Velocity',
                                          region = self.meshregion,
                                          cycle  = 30)

	Ucomp = V.data[:,0] - 0.9999

        #log V field
	self.logger.writeField(V,
                               region = self.meshregion,
                               desc   = 'Velocity at cycle 1')

	#get L1 error
	L1err    = self.measure.l1Error(Ucomp)

        #fail if L1err > VL1_error
        failmessage = 'In %s velocity L1_error = %10.3e is > %10.3e' %(self.methodName, L1err, self.tol['VL1_error'])

	self.failIf(L1err > self.tol['VL1_error'],
                    msg = failmessage)
