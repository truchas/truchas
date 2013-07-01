if __name__=='__main__':

    import os, sys
    #developers : please specify the location of your truchas directory
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

    def setTolerances(self):
        "define tolerances here"
        
        self.tol['PD_error'] = 4.0 #10e-05
        self.tol['VL1_error']= 10e-06

    def setDataStores(self):
        "defines all datastores to be used throughout this TestCase"

        #define testdata directories...
        self.testdirs     = ['static_drop_output']

        #define goldendata directories...
        self.goldendirs   = ['static_drop_golden']

    def setDefinitions(self):
        "creates all analytic results and spatial regions needed by this Capability TestCase"

        #define region based on interface position
        self.interface     = self.testdata[0].getRegion(selector=['VAR','VOF0001',[0.99,1.0]],
                                                        desc    = 'Interface defined by 0.99<VOF0001<1.0')

        self.mybox         = self.testdata[0].getRegion(selector=['S','[0.0,0.0,0.0->1.0,1.0,1.0]'],
                                                        desc    = 'Box region:[0.0,0.0,0.0->1.0,1.0,1.0]')
        
        self.goldinterface = self.goldendata[0].getRegion(selector=['VAR','VOF0001',[0.99,1.0]],
                                                          desc    = 'Interface defined by 0.99<VOF0001<1.0 in golden output')
                                                                 
    def testPartialPressureDrop(self):
        "tests interface pressure drop using 'Partial' measure"
        pass

    def testMaxPressureDrop(self):
	"tests interface pressure drop using 'Max' measure"

	Pinterface     = self.testdata[0].getField(field  = 'P',
                                                   time   = 0.00105,
                                                   region = self.mybox)

        #log P field
	self.logger.writeField(Pinterface,
                               region = self.mybox,
                               desc   = 'Pressure in box region')

        PL2err         = self.measure.l2Error(Pinterface.data)

        if PL2err > self.tol['PD_error']:
            failmessage = '\n In %s: Pressure drop L2_error = %10.3e is > %10.3e \n' %(self.methodName, PL2err, self.tol['PD_error'])
            self.logger.write(failmessage)
            self.fail(msg = failmessage)
        
        #self.tol['PD_error']
 
    def testSpuriousVelocities(self):
	"tests spurious velocity generation at interface of drop"

 	#get velocities on the interface region
	V     = self.testdata[0].getField(field  = 'Velocity',
                                          time   = 0.00105,
                                          region = self.interface)
	Ucomp = V.data[:,0]
	Vcomp = V.data[:,1]

        #get golden velocities on interface region
        Vgold = self.goldendata[0].getField(field  = 'Velocity',
                                            time   = 0.00105,
                                            region = self.goldinterface)

        #log V field
	self.logger.writeField(V,
                               region = self.interface,
                               desc   = 'Velocity at circle interface at cycle 1')

        #log Vgold field
	self.logger.writeField(Vgold,
                               region = self.goldinterface,
                               desc   = 'Golden velocity at circle interface at cycle 1')

	#get L1, L2, Linf error
	L1err    = self.measure.l1Error(V.data,Vgold.data)
        L2err    = self.measure.l2Error(V.data)
        Linferr  = self.measure.linfError(V.data)

        #fail if L1err > VL1_error
        failmessage = 'Velocity L1_error = %10.3e is > %10.3e' %(L1err, self.tol['VL1_error'])

	self.failIf(L1err > self.tol['VL1_error'],
                    msg = failmessage)




