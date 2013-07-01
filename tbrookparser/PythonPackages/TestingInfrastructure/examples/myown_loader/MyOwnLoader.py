#!/usr/bin/env python
"""

 MyOwnLoader

-----------------------------------------------------------------------------
   Purpose:
  
      Provides methods to create a customized TestSuite that in this case will test 
      extracting and plotting aborted residuals. Here the loader
      creates a TestSuite containing the following TestCases
      
      1. TruchasBasicRunFail
      2. GetAbortedResidualTest (lives in this directory)
      3. CheckPPGMV
  
   Public Interface:
  
      T = MyOwnLoader(dir,norerun,parallel_env,np,
                   compiler,compile_mode,version,
                   truchas_bin,truchas_exe)
      T.getSuite()
  
   Contains: MyOwnLoader

   Unit Test Block
  
   Author: Sharen Cummins (scummins@lanl.gov)
           Erin Iesulauro Barker (eibarker@lanl.gov)
-----------------------------------------------------------------------------
"""

import unittest
import os, sys, string, re, imp, types, fnmatch, platform

thisdir     = os.path.abspath(os.path.dirname(__file__))
testingdir  = thisdir + '/../../../TestingInfrastructure'
parserdir   = thisdir + '/../../../TBrookParser'

if __name__=='__main__':
    sys.path.append(testingdir)
    sys.path.append(parserdir)


from Loaders               import RunLoader, ModuleLoader, PPVisualisationLoader
from Interrogators         import DirectoryInterrogator

class MyOwnLoader(unittest.TestLoader):
    """
    Creation of a customized loader to create, extract and plot aborted
    residuals
    """

    def __init__(self,
                 no_rerun     = 0,
                 srcdir       = os.path.abspath(thisdir + '/../../../../regressiontests'),
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
        self.srcdir       = srcdir
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
            "do not load the TruchasBasicRunFail TestCase"

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

            self.logdir = string.split(RunSpecs.logdir,'/')[-1]

            #now load in the specific GetAbortedResidualTest that lives in this dir
            
            datastores        = []
            restartdatastores = []
            Y                 = ModuleLoader(datastores,
                                             restartdatastores,
                                             self.srcdir,
                                             self.currdir,
                                             RunSpecs,
                                             debug = self.debug)
            
            suite2 = Y.loadMyModuleTests(pattern='GetAbortedResidualTest')
            self.suite.addTest(suite2)
            
            #now check the GMV results of GetAbortedResidual TestCase 

            Z      = PPVisualisationLoader(self.currdir,
                                           self.debug,
                                           RunSpecs)
            
            suite3 = Z.loadPPGMVChecker(pattern='GetAbortedResidualTest.py')
            self.suite.addTest(suite3)

        else:
            "re-run will occur so recreate entire suite...."
            
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
            self.suite, self.logdir = X.loadBasicRunFailTests()

	    #now load in the specific GetAbortedResidualTest that lives in this dir
            
            Y      = ModuleLoader(datastores,
                                  restartdatastores,
                                  self.srcdir,
                                  self.currdir,
                                  RunSpecs,
                                  debug = self.debug)
            
            suite2 = Y.loadMyModuleTests(pattern='GetAbortedResidualTest')
            self.suite.addTest(suite2)

            #now check the GMV results of GetAbortedResidual TestCase 

            Z      = PPVisualisationLoader(self.currdir,
                                           self.debug,
                                           RunSpecs)
            
            suite3 = Z.loadPPGMVChecker(pattern='GetAbortedResidualTest.py')
            self.suite.addTest(suite3)


if __name__=='__main__':

    from Runners  import TestRunner
    
    currdir       = os.getcwd()
    srcdir        = os.path.abspath(thisdir + '/../../../../regressiontests'),
    truchas_bin   = testingdir + '/../../../bin'
    truchas_bin   = os.path.abspath(truchas_bin)

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

                loader       = MyOwnLoader(no_rerun,
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



