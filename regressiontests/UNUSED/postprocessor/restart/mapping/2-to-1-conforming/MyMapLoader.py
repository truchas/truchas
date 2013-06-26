#!/usr/bin/env python
"""

 MyMapLoader

-----------------------------------------------------------------------------
   Purpose:
  
      Provides methods to create TestSuites containing:
      1) Truchas BasicRun TestCase(s)
      2) Truchas PostProcessorRun TestCase(s)
  
      by implementing RunLoader, PPLoader, PPVisualisationLoader, PPDiagnosticLoader
  
   Public Interface:
  
      T = MyMapLoader(dir,norerun,parallel_env,np,
                                     compiler,compile_mode,truchas_bin,
                                     truchas_exe)
      T.getSuite()
  
   Contains: MyMapLoader

   Unit Test Block
  
   Author: Sharen Cummins (scummins@lanl.gov)
           Erin Iesulauro Barker (eibarker@lanl.gov)
-----------------------------------------------------------------------------
"""

import unittest
import os, sys, string, re, imp, types, fnmatch, platform

thisdir     = os.path.abspath(os.path.dirname(__file__))
testingdir  = thisdir + '/../../../../../tools/PythonPackages/TestingInfrastructure'
parserdir   = thisdir + '/../../../../../tools/PythonPackages/TBrookParser'

if __name__=='__main__':
    sys.path.append(testingdir)
    sys.path.append(parserdir)
    
from Loaders             import RunLoader, PPLoader
from Interrogators       import DirectoryInterrogator

class MyMapLoader(unittest.TestLoader):
    """
    Overriding the unittest loader to allow us to search for input files and run 
    TruchasBasicRunTest, 
    TruchasPostProcessorRunTest,
    CheckPPGMV, CheckPPEnsight, CheckPPTecPlot,CheckPPVTK (depending on *.mac file)
    CheckPPDiagnostics (depending on *.mac file)
    TruchasRestartRunTest
    TruchasBasicRunFailTest
    """

    def __init__(self,
                 no_rerun     = 0,
                 srcdir       = os.path.abspath(thisdir + '/../../../../../regressiontests'),
                 currdir      = '.',
                 debug        = 0,
		 parallel_env = 'serial',
		 np           = 1,
		 compiler     = 'lahey',
		 compile_mode = 'debug',
		 version      = '2.3',
		 truchas_bin  = '.',
		 truchas_exe  = None):

        self.no_rerun     = no_rerun
	self.debug        = debug
	self.currdir      = currdir
	self.parallel_env = parallel_env
	self.np           = np
	self.compiler     = compiler
	self.compile_mode = compile_mode
	self.version      = version
	self.truchas_bin  = truchas_bin
	self.truchas_exe  = truchas_exe
        self.toolsdir     = testingdir + '/../../'
        self.toolsdir     = os.path.abspath(self.toolsdir)

        self.suite        = unittest.TestSuite()
        self.logdir       = None

    def getSuite(self):

        if self.no_rerun:
            "do not load in any PostProcessor and CheckPP* TestCase(s)"

            #just ascertain BasicRunSpecs or RestartRunSpecs from the most
            #previous BasicRun/RestartRun using the dates on the logs. directory

            X   = DirectoryInterrogator(self.currdir)
            
            try:
                try:
                    RunSpecs  = X.getBasicRunSpecs()
                except:
                    RunSpecs  = X.getRestartRunSpecs()
            finally:
                err  = X.nobasicrunerr
                err += '\n'
                err += X.norestartrunerr
                assert err

            self.logdir = RunSpecs.logdir

        else:
            "re-run will occur so recreate suite...."
            
	    #first load in all the BasicRun Tests - 
	    #when run these will create multiple Truchas output dirs

            datastores        = []
            restartdatastores = []
            X                 = RunLoader(datastores,
                                          restartdatastores,
                                          self.currdir,
					  self.debug,
                                          self.parallel_env,
                                          self.np,
                                          self.compiler,
                                          self.compile_mode,
					  self.version,
                                          self.truchas_bin,
                                          self.truchas_exe)
	
            RunSpecs                = X.runspecs
            self.suite, self.logdir = X.loadBasicRunTests(myinputfile='runA.inp')
            suite1, self.logdir     = X.loadBasicRunTests(myinputfile='runB.inp')
            self.suite.addTest(suite1)
            #suite1, self.logdir     = X.loadBasicRunFailTests()
            #self.suite.addTest(suite1)

	    #now load in any PostProcessor Tests - 
	    #when run these will create restart, gmv, ensight, diagnostic files
            
            Y      = PPLoader(self.currdir,
			      self.debug,
                              self.toolsdir,
                              RunSpecs)
                              

            suite2 = Y.loadPostProcessorTests()
            self.suite.addTest(suite2)

	    """
            #now check the results of these postprocessor tests

            Z      = PPVisualisationLoader(self.currdir,
                                           self.debug,
                                           RunSpecs)
                                                
            suite3 = Z.loadPPGMVChecker()
            self.suite.addTest(suite3)
            suite4 = Z.loadPPEnsightChecker()
            self.suite.addTest(suite4)
            suite5 = Z.loadPPTecPlotChecker()
            self.suite.addTest(suite5)
            suite6 = Z.loadPPVTKChecker()
            self.suite.addTest(suite6)


            A      = PPDiagnosticLoader(self.currdir,
                                        self.debug,
                                        RunSpecs)
            suite7 = A.loadPPDiagnosticChecker()
            self.suite.addTest(suite7)
	    """


	    #now load in any RestartRun Tests - 
	    #when run these will create multiple Truchas restart output dirs
            
            suite8, self.logdir = X.loadRestartRunTests(myinputfile='runAB.inp')
            self.suite.addTest(suite8)


if __name__=='__main__':

    from Runners  import TestRunner
    
    currdir       = os.getcwd()
    srcdir        = os.path.abspath(thisdir + '/../../../../../regressiontests')
    truchas_bin   = testingdir + '/../../../bin'
    truchas_bin   = os.path.abspath(truchas_bin)

    from PYTHONutils import uTestOpts
    
    # file, debug, output prefix, binary, ascii, clean
    dfltdir  = testingdir + '/TestCases/'
    opts = uTestOpts('dc', 
                     defaults = {'d' : False,
                                 'c' : False },
                     actions  = {'d' : 'store_true',
                                 'c' : 'store_true'},
                     dir      = dfltdir)
    (opt,args) = opts.parse_args()
    opts.header(__file__,None)

    try:

        if opt.c: #clean

	    print "Cleaning up from a past test run...\n"
	    cmd = 'rm *.inp'
	    os.system(cmd)
	    cmd = 'rm -r for_debugging'
	    os.system(cmd)
	    cmd = 'rm -r outputs*'
	    os.system(cmd)

	else:
		
            for exefile in os.listdir(truchas_bin):
         	if os.path.isfile(os.path.join(truchas_bin,exefile)) and exefile != 'truchas':
        	    vals  = string.splitfields(exefile,'.')
        	    if vals[0][2:] == string.lower(platform.system()):
                        
        	        parallel_env = vals[3]
        	        compiler     = vals[2]
        	        compile_mode = vals[4][:-2]
        		version      = vals[4][len(vals[4])-1] + '.' + vals[5]
        		np           = 1
        		no_rerun     = 0
        		debug        = 1
        
                        loader       = MyMapLoader(no_rerun,
                                                   srcdir,
                                                   currdir,
                                                   debug,
                                                   parallel_env,
                                                   np,
                                                   compiler,
                                                   compile_mode,
                                                   version,
                                                   truchas_bin)
                        loader.getSuite()
                    
                        fp      = sys.stdout
                        runner  = TestRunner(stream=fp,verbosity=2)
                        result  = runner.run( loader.suite )
        
    except:

	print "---> Test failed in some aspect <---"
	print "\nNature-Of-Error:",sys.exc_info()[0],"\n"
	if opt.d: raise
        
