#!/usr/bin/env python
"""
 PPDiagnosticLoader

-----------------------------------------------------------------------------
   Purpose:
  
      Provides methods to create TestSuites containing:
      1) CheckPPDiagnostics TestCase(s)
  
   Public Interface:
  
      T = PPDiagnosticLoader(currdir,RunSpecs) 
      T.loadPPDiagnosticChecker
  
   Contains:
      class PPDiagnosticLoader(unittest.TestLoader)
             __init__(currdir,RunSpecs) 
             loadPPDiagnosticChecker(diagfiles,mymacfile)
             __rjustln(columns,s)
             __PPCheckerLoadDetails(macfile)
             
      class RunThisPPDiagTest
             
   Unit Test Block
  
   Author: Sharen Cummins (scummins@lanl.gov)
           Erin Iesulauro Barker (eibarker@lanl.gov)
-----------------------------------------------------------------------------
"""

import unittest
import os, sys, string, types, fnmatch, platform

if __name__=='__main__':
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    testingdir  = thisdir + '/../'
    sys.path.append(testingdir)
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

from Loggers           import TestLoggers
from TestCases         import CheckPPDiagnostics

class PPDiagnosticLoader(unittest.TestLoader):
    """
    Overriding the unittest TestLoader to allow us to search for 
    and load the following TestCase classes:

    -CheckPPDiagnostics (if stat.mac, query.mac or probe.mac file(s) exist(s))
    """

    def __init__(self,currdir,debug,RunSpecs): 

	self.currdir           = currdir
	self.debug             = debug
	self.logdir            = 'for_debugging'
	self.RunSpecs          = RunSpecs        
	#for formatting debug file
	self.column1           = 5
	self.column2           = 10

        if self.debug:
            if not os.path.isdir(os.path.abspath(os.path.join(self.currdir,self.logdir))):            
                os.mkdir(os.path.abspath(os.path.join(self.currdir,self.logdir)))
                

    def loadPPDiagnosticChecker(self,
                                diagfiles  = ['probe.mac','query.mac','stat.mac'],
                                mymacfile  = None):
	"""
	Searches through the current directory looking for any files containing 
	'mymacfile' or 'diagfiles' then loads the test(s) into a suite.
	"""
	suite       = unittest.TestSuite()
        loadername  = 'PPDiagnosticCheckLoader'
        
        if self.debug:
            loadername  = os.path.abspath('%s/%s/%s' %(self.currdir,self.logdir,'PPDiagnosticCheckLoader'))
            
	self.logger = TestLoggers(loadername,debugging=self.debug).getThem()

	if (mymacfile != None):
	    self.macfile  = mymacfile

	    self.ThisTest = CheckPPDiagnostics(self.currdir,
                                               self.RunSpecs) 

	    suite.addTest(self.ThisTest)

            self.logger.debug(self.__PPCheckerLoadDetails())
	    if self.debug: print self.ThisTest.str()


	else:
	    #perform search for macro files in this directory
	    files           = os.listdir(self.currdir)
	    for f in files:
		fullname = os.path.join(self.currdir, f)
                if f in diagfiles:
		    self.macfile  = f
		    self.ThisTest = CheckPPDiagnostics(self.currdir,
                                                       self.RunSpecs) 

		    suite.addTest( self.ThisTest )

                    self.logger.debug(self.__PPCheckerLoadDetails())
	            if self.debug: print self.ThisTest.str()

        info  = '\n'
	tmp   = 'Suite to be run          :' + str(suite)
        info += self.__rjustln(self.column2, tmp)        

        self.logger.debug(info)
        
        return suite


    def __rjustln(self,column,s):

        tmp = s.rjust(column + len(s))
        tmp += '\n'

        return tmp

    def __PPCheckerLoadDetails(self):

	info  = '\n'
	tmp   = '*Loaded CheckPPDiagnostics TestCase'
        info += self.__rjustln(self.column1, tmp)        
        tmp   = self.ThisTest.str()
        info += self.__rjustln(self.column1, tmp)        
        info += '\n'

	return info


    
if __name__=='__main__':

    currdir  = os.getcwd()
    toolsdir = currdir + '/../../../'
    toolsdir = os.path.abspath(toolsdir)

    truchas_bin       = testingdir + '/../../../bin'
    truchas_bin       = os.path.abspath(truchas_bin)

    from Runners import TestRunner, RunTimeSpecs        
    from PYTHONutils import uTestOpts
    
    # file, debug, output prefix, binary, ascii, clean
    dfltdir  = '../TestCases/static_drop/static_drop_golden/'
    dfltfile = 'static_drop'
    opts = uTestOpts('dc', 
                     defaults = {'d' : False,
                                 'c' : False},
                     actions  = {'d' : 'store_true',
                                 'c' : 'store_true'},
                     dir      = dfltdir)
    (opt,args) = opts.parse_args()
    opts.header(__file__,None)

    try:

        if opt.c: #clean

	    print "Cleaning up from a past test run...\n"
	    cmd = 'rm -r outputs*'
	    os.system(cmd)
	    cmd = 'rm -r for_debugging'
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
                        logfile      = logdir + '/static_drop_logs/BasicRun.log'
                        LastRunSpecs = RunTimeSpecs(logdir,logfile,parallel_env,np,compiler,compile_mode)
        
                        X       = PPDiagnosticLoader(currdir,debug,LastRunSpecs)
                        suite   = X.loadPPDiagnosticChecker()
                        runner  = TestRunner(verbosity=2)
                        result  = runner.run( suite )
            
    except:

	print "---> Test failed in some aspect <---"
	print "\nNature-Of-Error:",sys.exc_info()[0],"\n"
	if opt.d: raise
