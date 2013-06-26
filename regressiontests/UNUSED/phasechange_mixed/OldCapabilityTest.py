if __name__=='__main__':

    import os, sys
    #developers : please specify the location of your truchas checkout 
    #             directory 'truchasdir'
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
        
        self.tol['error'] = 3e-07
        self.tol['Herror'] = 3e-06

    def setDataStores(self):
        "defines all datastores to be used throughout this TestCase"

        #define testdata output directories..these result from running a BasicRun TestCase
        self.testdirs = ['phasechange_mixed_output']
        self.goldendirs = ['phasechange_mixed_golden']

    def setDefinitions(self):
        "defines analytic expressions and regions needed by this Capability TestCase"

	#set values from input file
	self.timerange = [0,200]
	self.dtstep    = 1.0
	self.times     = [0,100,200]
	
        #get region representing the entire mesh
        self.meshregion = self.testdata[0].getRegion(selector=['Default'],
			                             desc='Region representing the entire Truchas mesh')
	
        #get region representing liquid phase Vof1 = 1
	self.liquid     = self.testdata[0].getRegion(selector=['VAR','VOF0001',[0.1,1.0]],
				                     desc = 'liquid region')
        #get region representing liquid phase Vof1 = 1
	self.solid      = self.testdata[0].getRegion(selector=['VAR','VOF0002',[0.1,1.0]],
				                     desc = 'solid region')


    def testEnthalpy(self):
	"tests enthalpy of each cell at each output time"

	timesteps = self.testdata[0].getTimeSteps(timerange = self.timerange,
			                          dtstep    = self.dtstep)

	for step in range(len(timesteps)): 
	    cycle = timesteps[step].cycle

            #compare enthalpy for the first cycle
            H = self.testdata[0].getField(field  = 'Z_ENTHALPY',
			                  cycle  = cycle,
				          region = self.meshregion)
	    self.logger.writeField(H,
			           region = self.meshregion,
			           desc   = 'Enthalpy of the system')
            Hgold = self.goldendata[0].getField(field  = 'Z_ENTHALPY',
			                      cycle  = cycle,
				              region = self.meshregion)

	    #calcaulte the error
	    err       = self.measure.l2Error(H.data,Hgold.data)
	    logmessage = '\n L2 Enthalpy error = %10.3e \n' % (err)
	    self.logger.write(logmessage)
		    
	    #fail if err > error
	    failmessage = 'In %s: Enthalpy L2_error = %10.3e is > %10.3e for cycle %5.3f ' % ( self.methodName, err, self.tol['Herror'], cycle)

	    self.failIf(err>self.tol['Herror'], msg = failmessage)
		

    def testTemperature(self):
	"tests temperature of each cell at each output time"

	timesteps = self.testdata[0].getTimeSteps(timerange = self.timerange,
			                          dtstep    = self.dtstep)

	for step in range(len(timesteps)): 
	    cycle = timesteps[step].cycle

            #compare enthalpy for the first cycle
            T = self.testdata[0].getField(field  = 'Z_TEMP',
			                  cycle  = cycle,
				          region = self.meshregion)
	    self.logger.writeField(T,
			           region = self.meshregion,
			           desc   = 'Temperature of the system')
            Tgold = self.goldendata[0].getField(field  = 'Z_TEMP',
			                      cycle  = cycle,
				              region = self.meshregion)
	    err       = self.measure.l2Error(T.data,Tgold.data)
	    logmessage = '\n L2 Temperature error = %10.3e \n' % (err)
	    self.logger.write(logmessage)

	    #fail if err > error
	    failmessage = 'In %s: Temperature L2_error = %10.3e is > %10.3e for cycle %5.3f ' % ( self.methodName, err, self.tol['error'], cycle)

	    self.failIf(err>self.tol['error'], msg = failmessage)


    def testVOF(self):
	"tests temperature of each cell at each output time"

	timesteps = self.testdata[0].getTimeSteps(timerange = self.timerange,
			                          dtstep    = self.dtstep)

	for step in range(len(timesteps)): 
	    cycle = timesteps[step].cycle

            #compare enthalpy for the first cycle
            VOF = self.testdata[0].getField(field  = 'VOF0001',
			                  cycle  = cycle,
				          region = self.meshregion)
	    self.logger.writeField(VOF,
			           region = self.meshregion,
			           desc   = 'VOF0001 of the system')
            VOFgold = self.goldendata[0].getField(field  = 'VOF0001',
			                          cycle  = cycle,
				                  region = self.meshregion)
	    err       = self.measure.l2Error(VOF.data,VOFgold.data)
	    logmessage = '\n L2 VOF error = %10.3e \n' % (err)
	    self.logger.write(logmessage)

	    #fail if err > error
	    failmessage = 'In %s: VOF L2_error = %10.3e is > %10.3e for cycle %5.3f ' % ( self.methodName, err, self.tol['error'], cycle)

	    self.failIf(err>self.tol['error'], msg = failmessage)
