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

try:
   import numpy.oldnumeric as Numeric
except ImportError:
   import Numeric
except:
   raise

from TestCases import TruchasCapabilityTest

class CapabilityTest(TruchasCapabilityTest):
    "pure material solidification CapabilityTest" 

    def setTolerances(self):
        "define tolerances here"
        
        self.tol['SF_error'] = 2.0e-02

    def setDataStores(self):
        "defines testdata and goldendata output directories needed by this Capability TestCase"

        #define testdata output directories..these are created from running a BasicRun TestCase
        self.testdirs     = ['purematerial_solidification_output']

    def setDefinitions(self):
	"developer defined regions to be used throughout this Capability TestCase"

	#get region representing plane x=3.9 from testdata store
	self.myplane       = self.testdata[0].getRegion(selector = ['S','[3.9,0.0,0.0->3.9,3.9,0.1]'],
                                                        desc     = 'Region:x=3.9 plane')

    def testSolidFractionsOnPlane(self):
	"tests solid fractions along x=3.9 at various times 0.5<=t<=6.0"

	#get list of timesteps from Truchas testdata output that satisfy 0.5<=t<=6.0 at intervals=0.5
	these_timesteps  = self.testdata[0].getTimeSteps(timerange = [0.5,6.0],
                                                         dtstep    = 0.5)

        #obtain golden SFtotal (obtained from Lazaridis,Int.J.Heat Mass Transfer,Vol13,(1970),pp1459-1477)
        SFtotal_golden = Numeric.array([0.18,0.26,0.32,0.37,0.41,0.45,0.49,0.53,0.56,0.6,0.63,0.66],'d')

	for timestep in these_timesteps:

	    thistime   = timestep.time
            index      = these_timesteps.index(timestep)

            #get solidfraction field from testdata store on self.myplane region
	    SFrac      = self.testdata[0].getField(field  = 'VOF0002',
                                                   time   = thistime,
                                                   region = self.myplane)

	    #log SFrac field
	    description = 'Test solid fraction at time t =%5.2f' %(thistime)
	    self.logger.writeField(SFrac,
                                   region = self.myplane,
                                   desc   = description)


	    #calculate summation of SFrac on plane and normalize it  
	    SFtotal     = self.measure.totalValue(SFrac.data)/len(SFrac.data)

            message     = '\n Solid fraction along x=3.9 is %10.3e \n' %(SFtotal)
            self.logger.write(message)

            #calculate the l2error of SFtotal with golden output
            err         = self.measure.percentError(SFtotal,
                                                    SFtotal_golden[index])

	    #fail if err > SF_error
            failmessage = 'In %s: solid fraction sum = %10.3e is > %10.3e' %(self.methodName, err, self.tol['SF_error'])

	    self.failIf(err>self.tol['SF_error'],
                        msg = failmessage)

