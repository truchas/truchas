#!/usr/bin/env python
"""
 ModuleLoader

-----------------------------------------------------------------------------
   Purpose:
  
      Provides methods to create TestSuites containing:
      1) Truchas Module (Capability) TestCase(s)
  
   Public Interface:
  
      T = ModuleLoader(datastores,restartdatastores,dir,toolsdir,BasicRunSpecs,RestartRunSpecs,debug)
      T.loadMyModuleTests(pattern)
  
   Contains:
      class ModuleLoader(unittest.TestLoader)
             __init__(datastores,restartdatastores,currdir,BasicRunSpecs,RestartRunSpecs,debug)
             loadMyModuleTests(pattern)
             getTestCaseNames(testCaseClass) (overloaded from unittest.TestLoader)
             loadTestsFromTestCase(testCaseClass)
             loadTestsFromModule(module)
             __rjustln(columns,s)
             __nonExistentModule()
             __moduleLoadDetails()
             
      class RunThisModuleTest
             
   Unit Test Block
  
   Author: Sharen Cummins (scummins@lanl.gov)
           Erin Iesulauro Barker (eibarker@lanl.gov)
-----------------------------------------------------------------------------
"""

import unittest
import imp
import os, sys, string, types, fnmatch, platform

if __name__=='__main__':
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    testingdir  = thisdir + '/../'
    sys.path.append(testingdir)
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)
    from Runners       import TestRunner, RunTimeSpecs

from Loggers           import TestLoggers
from TestCases         import TruchasCapabilityTest
from DataHandlers      import getDataCreator

class ModuleLoader(unittest.TestLoader):
    """
    Overriding the unittest TestLoader to allow us to search for 
    and load the following TestCase classes:

    -CapabilityTest
    """

    def __init__(self,datastores,
                 restartdatastores,
                 srcdir,
                 currdir,
                 BasicRunSpecs,
                 RestartRunSpecs=None,
                 debug=0): 

	self.currdir           = currdir
	self.srcdir            = srcdir
	self.BasicRunSpecs     = BasicRunSpecs
	self.RestartRunSpecs   = RestartRunSpecs
	self.macFilePat        = '*.mac'
	self.logdir            = 'for_debugging'
	#for formatting debug file
	self.column1           = 5
	self.column2           = 10
        self.datastores        = datastores
        self.restartdatastores = restartdatastores
        self.goldendatastores  = [] 
        self.debug             = debug

        if self.debug:
            if not os.path.isdir(os.path.abspath(os.path.join(self.currdir,self.logdir))):            
                os.mkdir(os.path.abspath(os.path.join(self.currdir,self.logdir)))

    def loadMyModuleTests(self,pattern='Capability'):
	"""
	Searches through the current directory looking for module name 
	containing "pattern" and then loads the test(s) into a suite.
	Typical patterns to match right now are Capability.
	"""
	self.suite            = unittest.TestSuite()
	self.pattern          = pattern
	self.testMethodPrefix = 'test'
        loadername            = self.pattern

        if self.debug:
            loadername            = os.path.abspath('%s/%s/%sLoader' %(self.currdir,self.logdir,self.pattern))
            
	self.logger           = TestLoggers(loadername,debugging=self.debug).getThem()
        cwd                   = os.getcwd()
        os.chdir(self.currdir)

	thismodule     = None
	for filename in os.listdir(self.currdir):
	    iend       = string.rfind(filename,".")
            if ".~" not in filename:
                modulename = filename[0:iend]
                if self.pattern == modulename:
                    thismodule  = modulename

	#check for module existence...

	if thismodule == None:
	    self.logger.debug(self.__nonExistentModule())
	    return self.suite

	self.file, self.pathname, self.desc = imp.find_module(thismodule,[self.currdir])
	self.testmod                        = imp.load_module(thismodule, 
							      self.file, 
							      self.pathname, 
							      self.desc)

	self.suite.addTest( self.loadTestsFromModule(self.testmod)) 

	self.logger.debug(self.__moduleLoadDetails())

        info  = '\n'
        tmp   = 'Suite to be run :' + str(self.suite)        
        info += self.__rjustln(self.column2, tmp)        

        self.logger.debug(info)

        os.chdir(cwd)

	return self.suite

    def getTestCaseNames(self, testCaseClass):
        """Return a sorted sequence of method names found within testCaseClass
        """
        testFnNames = filter(lambda n,p=self.testMethodPrefix: n[:len(p)] == p,
                             dir(testCaseClass))
        return testFnNames


    def loadTestsFromTestCase(self, testCaseClass):
        """Return a suite of all tests cases contained in testCaseClass"""
	suite          = unittest.TestSuite()
        numberoftests  = len(self.getTestCaseNames(testCaseClass))
        
	for i in range(len(self.getTestCaseNames(testCaseClass))):
            testnumber = i+1
            T          = testCaseClass(self.datastores,
                                       self.restartdatastores,
                                       self.goldendatastores,
                                       numberoftests,
                                       testnumber,
                                       self.BasicRunSpecs, 
                                       self.RestartRunSpecs, 
                                       methodName=self.getTestCaseNames(testCaseClass)[i])
            if testnumber == 1:
                T.setDataStores()
                Tdir = os.path.abspath(self.BasicRunSpecs.logdir + '/../')

                for goldir in T.goldendirs:
                    G  = getDataCreator(dir = os.path.join(self.srcdir,goldir))
                    self.goldendatastores.append(G)
                    
	    suite.addTest(T)

	return suite

    def loadTestsFromModule(self, module):
        """Overloaded from unittest:Return a suite of all tests cases contained in the given module"""

        tests = []
        for name in dir(module):
            obj = getattr(module, name)
            if (isinstance(obj, (type, types.ClassType)) and 
		issubclass(obj, TruchasCapabilityTest)):
                tests.append(self.loadTestsFromTestCase(obj))
	return unittest.TestSuite(tests)

    def __rjustln(self,column,s):

        tmp = s.rjust(column + len(s))
        tmp += '\n'

        return tmp

    def __nonExistentModule(self):

	info  = '\n'
	tmp   = '*No module exists'
        info += self.__rjustln(self.column1, tmp)        
	tmp   = 'Pattern         :' + self.pattern
        info += self.__rjustln(self.column1, tmp)        
	tmp   = 'Current dir     :' + self.currdir
        info += self.__rjustln(self.column1, tmp)        
	tmp   = 'No tests will be loaded and run from this directory.'
        info += self.__rjustln(self.column1, tmp)        
	info += '\n'
	return info

    def __moduleLoadDetails(self):

	info  = '\n'
	tmp   = '*Loaded Capability TestCase'
        info += self.__rjustln(self.column1, tmp)        
	tmp   = 'Filename        :' + self.pathname
        info += self.__rjustln(self.column1, tmp)        
	tmp   = 'Module          :' + str(self.testmod)
        info += self.__rjustln(self.column1, tmp)        
	tmp   = 'Basic RunSpecs : ' + self.BasicRunSpecs.str()
        info += self.__rjustln(self.column1, tmp)
	if self.RestartRunSpecs != None:
	    tmp   = 'Restart RunSpecs : ' + self.RestartRunSpecs.str()
            info += self.__rjustln(self.column2, tmp)        
	info += '\n'

	return info


if __name__=='__main__':

    test_type='CapabilityTest'
    
    currdir  = os.getcwd()
    toolsdir = currdir + '/../../../'
    toolsdir = os.path.abspath(toolsdir)
    truchas_bin       = testingdir + '/../../../bin'
    truchas_bin       = os.path.abspath(truchas_bin)

    datastores        = []
    restartdatastores = []

    from PYTHONutils import uTestOpts
    
    # file, debug, output prefix, binary, ascii, clean
    opts = uTestOpts('dc', 
                     defaults = {'d' : False,
                                 'c' : False},
                     actions  = {'d' : 'store_true',
                                 'c' : 'store_true'},
                     dir      = currdir)
    (opt,args) = opts.parse_args()
    opts.header(__file__,None)

    try:

        if opt.c: #clean

	    print "Cleaning up from a past test run...\n"
	    cmd = 'rm -r for_debugging'
	    os.system(cmd)
	    cmd = 'rm -r outputs*'
	    os.system(cmd)
	    cmd = 'rm *.inp'
	    os.system(cmd)

	else:
		
            for exefile in os.listdir(truchas_bin):
         	if os.path.isfile(os.path.join(truchas_bin,exefile)) and exefile != 'truchas':
        	    vals  = string.splitfields(exefile,'.')
        	    if vals[0][2:] == string.lower(platform.system()):
                        
        	        parallel_env   = vals[3]
        	        compiler       = vals[2]
        	        compile_mode   = vals[4][:-2]
        		np             = 1
        		debug          = 1
        
        	        logdir       = 'outputs'
        		for val in vals: logdir += '.'+val
                        logfile      = logdir+'/static_drop_logs/BasicRun.log'
                	LastRunSpecs = RunTimeSpecs(logdir,logfile,parallel_env,np,compiler,compile_mode)
                
        		srcdir   = '' # fix this later to be correct
        	        X        = ModuleLoader(datastores,restartdatastores,srcdir,currdir,LastRunSpecs,debug=1)
                	suite    = X.loadMyModuleTests(pattern=test_type)
                	runner   = TestRunner(verbosity=2)
                	result   = runner.run( suite )    

    except:

	print "---> Test failed in some aspect <---"
	print "\nNature-Of-Error:",sys.exc_info()[0],"\n"
	if opt.d: raise
