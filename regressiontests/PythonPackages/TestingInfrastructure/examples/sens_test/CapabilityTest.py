if __name__=='__main__':

    import os, sys
    #developers : please specify the location of your truchas checkout 
    #             directory 'truchasdir'
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
        
        self.tol['error'] = 10e-07

    def setDataStores(self):
        "defines all datastores to be used throughout this TestCase"

        #define testdata output directories..these result from running a 
	#  BasicRun TestCase
        self.testdirs   = ['sens_test_output']
        self.goldendirs = ['sens_test_golden']

    def setDefinitions(self):
        "defines analytic expressions and regions needed by this Capability TestCase"

        #get region representing the entire mesh
        self.meshregion        = self.testdata[0].getRegion(
                                                            selector = ['Default'],
                                                            desc     = 'Region representing the entire Truchas mesh')
	
        #get region representing each sensitivity function location as defined in input
	#    first cell in the bottom left hand corner [0,0,0]
	self.bottomleftcorner  = self.testdata[0].getRegion(
                                                            selector = ['ID','CELL','1'],
				                            desc     = 'bottom left corner of liquid')
        
	#    last cell in the liquid region [0.49,0,0]
	self.bottomrightcorner = self.testdata[0].getRegion(
			                                    selector = ['ID','CELL','4'],
				                            desc     = 'bottom right corner of liquid')
        
	#    edge cell between mold1 and mold2 [1.1,0,0]
	self.betweenmolds      = self.testdata[0].getRegion(
			                                    selector = ['ID','CELL','9'],
				                            desc     = 'bottom corner between mold1 and mold2')


    def testTemperatureSensitivity(self):
	"tests temperature sensitivity"

	#set the number of variables to test
	#the input file defines 9 variables.  
	#test = 9
	#rather than just set the value, this extracts the number from Truchas
	tmp  = self.testdata[0].getField(field  = "nsens_design_variables",
                                         time   = 0.0,
                                         region = self.bottomleftcorner)
	test = tmp.data[0]
	
	#get list of timesteps from Truchas
	these_timesteps = self.testdata[0].getTimeSteps(timerange = [0,9000],
			                                dtstep    = 100)
	#loop over each time step and compare the values
	for timestep in these_timesteps:
	    
	    thistime = timestep.time

	    #loop over the location where temperature sentivity functions 
	    #are defined
	    for r in (self.bottomleftcorner, self.bottomrightcorner):

		#loop over each variable for the current location
  	        for i in xrange(test):

		    # increment field variable name
		    f = 'temp_sens_000'+str(i+1)

	            tsens = self.testdata[0].getField(field  = f, 
                                                      time   = thistime,
					              region = r)

	            description = 'Test tsens%s at time t=%5.3f' %(i,thistime)
	            self.logger.writeField(tsens,
                                           region = self.bottomleftcorner, 
                                           desc   = description)

	            #get golden tsnes field for the current variable
	            tsens_golden = self.goldendata[0].getField(field  = f, 
                                                               time   = thistime,
                                                               region = r)

	            #calculate the % error
	            err = self.measure.l1Error(tsens.data, tsens_golden.data)
	     
	            #fail if err > error
	            failmessage = 'In %s: tsens%s L1_error = %10.3e is > %10.3e for timestep %5.3f at %s' % (
                                   self.methodName, i, err, self.tol['error'],thistime, r.name)

	            self.failIf(err>self.tol['error'], msg = failmessage)
		
		
    def testEnthalpySensitivity(self):
	"tests enthalpy sensitivity"

	# set the number of variables to test
	#the input file defines 9 variables.  
	#test = 9
	#rather than just set the value, this extracts the number from Truchas
	tmp  = self.testdata[0].getField(field  = "nsens_design_variables",
                                         time   = 0.0,
                                         region = self.bottomleftcorner)
	test = tmp.data[0]
	
	#get list of timesteps from Truchas
	these_timesteps = self.testdata[0].getTimeSteps(timerange = [0,9000],
			                                dtstep    = 100)
	#loop over each time step and compare the values
	for timestep in these_timesteps:
	    
	    thistime = timestep.time

	    #loop over the location where temperature sentivity functions 
	    #   are defined
	    for r in (self.bottomleftcorner, self.bottomrightcorner, self.betweenmolds):

		#loop over each variable for the current location
	        for i in xrange(test):

		    # increment field variable name
		    f     = 'temp_sens_000'+str(i+1)

	            esens = self.testdata[0].getField(field  = f, 
			                              time   = thistime,
					              region = r)

	            description = 'Test esens%s at time t=%5.3f' %(i,thistime)
	            self.logger.writeField(esens,
		  	                   region = r, 
				           desc   = description)

	            #get golden esnes field for each variable
	            esens_golden = self.goldendata[0].getField(field  = f, 
                                                               time   = thistime,
                                                               region = r)

	            #calculate the % error
	            err = self.measure.l1Error(esens.data, esens_golden.data)
	     
	            #fail if err > error
	            failmessage = 'In %s: esens%s L1_error = %10.3e is > %10.3e for timestep %5.3f at %s' % (
                                   self.methodName, i, err, self.tol['error'],thistime, r.name)

	            self.failIf(err>self.tol['error'], msg = failmessage)
		
 
    def testVOF1Sensitivity(self):
	"tests VOF1 sensitivity"

	# set the number of variables to test
	#the input file defines 9 variables.  
	#test = 9
	#rather than just set the value, this extracts the number from Truchas
	tmp  = self.testdata[0].getField(field  = "nsens_design_variables",
                                         time   = 0.0,
                                         region = self.bottomleftcorner)
	test = tmp.data[0]
	
	#get list of timesteps from Truchas
	these_timesteps = self.testdata[0].getTimeSteps(timerange = [0,9000],
			                                dtstep    = 100)
	#loop over each time step and compare the values
	for timestep in these_timesteps:
	    
	    thistime = timestep.time

	    #loop over the location where temperature sentivity functions
	    #are defined
	    for r in (self.bottomleftcorner, self.bottomrightcorner):

		#loop over each variable for the current location
	        for i in xrange(test):

		    # increment field variable name
		    f    = 'volume_fraction_sens_000'+str(i+1)+'0001'

	            VOF1 = self.testdata[0].getField(field  = f, 
			                             time   = thistime,
					             region = r)

	            description = 'Test VOF1_%s at time t=%5.3f' %(i,thistime)
	            self.logger.writeField(VOF1,
		  	                   region = r, 
				           desc   = description)

	            # get golden tsnes field for each variable
	            VOF1_golden = self.goldendata[0].getField(field  = f, 
                                                              time   = thistime,
                                                              region = r)

	            #calculate the % error
	            err = self.measure.l1Error(VOF1.data, VOF1_golden.data)
	     
	            #fail if err > error
	            failmessage = 'In %s: VOF1_%s L1_error = %10.3e is > %10.3e for timestep %5.3f at %s' % (
                                   self.methodName, i, err, self.tol['error'],thistime,r.name)

	            self.failIf(err>self.tol['error'], msg = failmessage)
		
    def testVOF2Sensitivity(self):
	"tests VOF2 sensitivity"

	# set the number of variables to test
	#the input file defines 9 variables.  
	#test = 9
	#rather than just set the value, this extracts the number from Truchas
	tmp  = self.testdata[0].getField(field  = "nsens_design_variables",
                                         time   = 0.0,
                                         region = self.bottomleftcorner)
	test = tmp.data[0]
	
	#get list of timesteps from Truchas
	these_timesteps = self.testdata[0].getTimeSteps(timerange = [0,9000],
			                                dtstep    = 100)
	#loop over each time step and compare the values
	for timestep in these_timesteps:
	    
	    thistime = timestep.time

	    #loop over the location where temperature sentivity functions 
	    #are defined
	    for r in (self.bottomleftcorner, self.bottomrightcorner):
		    
		#loop over each variable for the current location
	        for i in xrange(test):
		    
		    # increment field variable name
		    f = 'volume_fraction_sens_000'+str(i+1)+'0002'

	            VOF2 = self.testdata[0].getField(field  = f, 
			                             time   = thistime,
					             region = r)

	            description = 'Test VOF2_%s at time t=%5.3f' %(i,thistime)
	            self.logger.writeField(VOF2,
		  	                   region = r, 
				           desc   = description)

	            #get golden tsnes field for each variable
	            VOF2_golden = self.goldendata[0].getField(field  = f, 
                                                              time   = thistime,
                                                              region = r)

	            #calculate the % error
	            err = self.measure.l1Error(VOF2.data, VOF2_golden.data)
	     
	            #fail if err > error
	            failmessage = 'In %s: VOF2_%s L1_error = %10.3e is > %10.3e for timestep %5.3f at %s' % (
                                   self.methodName, i, err, self.tol['error'],thistime, r.name)

	            self.failIf(err>self.tol['error'], msg = failmessage)

    
    def testUserDefSensitivity(self):
	"tests user defined function sensitivity"

	# set the number of variables to test
	#the input file defines 9 variables.  
	#test = 10 
	#rather than just set the value, this extracts the number from Truchas
	tmp  = self.testdata[0].getField(field  = "nsens_functions",
                                         time   = 0.0,
                                         region = self.bottomleftcorner)
	test = tmp.data[0]
	
	#get list of timesteps from Truchas
	these_timesteps = self.testdata[0].getTimeSteps(timerange = [0,9000],
			                                dtstep    = 100)

	#This function is only used in one location
	r= self.betweenmolds

	#loop over each time step and compare the values
	for timestep in these_timesteps:
	    
	        thistime = timestep.time

		#loop over each test variable
	        for i in xrange(test):
		    
		    #there are 10 defined variables.  
		    #this changes the string slightly
		    if i>8: f = 'sens_function_gradient_00'+str(i+1)
		    else:   f = 'sens_function_gradient_000'+str(i+1)

	            #get the value for this function
	            user = self.testdata[0].getField(field  = f, 
			                             time   = thistime,
					             region = r)

	            description = 'Test UserFunc_%s at time t=%5.3f' % (
				  i,thistime)
                    
	            self.logger.writeField(user,
		  	                   region = r, 
				           desc   = description)

	            #get golden field for the variable
	            user_golden = self.goldendata[0].getField(field  = f, 
                                                              time   = thistime,
                                                              region = r)

	            #calculate the % error
	            err = self.measure.l1Error(user.data, user_golden.data)
	     
	            #fail if err > error
	            failmessage = 'In %s: UserFunc_%s L1_error = %10.3e is > %10.3e for timestep %5.3f at %s' % (
                                   self.methodName, i, err, self.tol['error'],thistime, r.name)

	            self.failIf(err>self.tol['error'], msg = failmessage)


