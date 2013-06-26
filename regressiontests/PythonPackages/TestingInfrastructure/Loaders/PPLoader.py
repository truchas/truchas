#!/usr/bin/env python
"""
 PPLoader

-----------------------------------------------------------------------------
   Purpose:
  
      Provides methods to create TestSuites containing:
      1) Truchas PostProcessorRun TestCase(s)
  
   Public Interface:
  
      T = PPLoader(dir,toolsdir,BasicRunSpecs,RestartRunSpecs)
      T.loadPostProcessorTests(macfile)
  
   Contains:
      class PPLoader(unittest.TestLoader)
             __init__(currdir,toolsdir,BasicRunSpecs,RestartRunSpecs)
             loadPostProcessorTests(macfile)
             __rjustln(columns,s)
             __postProcessorLoadDetails(macfile)
             
      class RunThisPostProcessorTest
             
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
    from Runners           import TestRunner, RunTimeSpecs

from Loggers           import TestLoggers
from TestCases         import TruchasPostProcessorRunTest

class PPLoader(unittest.TestLoader):
    """
    Overriding the unittest TestLoader to allow us to search for 
    and load the following TestCase classes:

    -TruchasPostProcessorRunTest (if *.mac file exists) 
    """

    def __init__(self,currdir,debug,toolsdir,BasicRunSpecs,RestartRunSpecs=None): 

	self.currdir           = currdir
	self.debug             = debug
	self.postprodir        = toolsdir
	self.BasicRunSpecs     = BasicRunSpecs
	self.RestartRunSpecs   = RestartRunSpecs
	self.macFilePat        = '*.mac'
	self.logdir            = 'for_debugging'
	#for formatting debug file
	self.column1           = 5
	self.column2           = 10

        if self.debug:
            if not os.path.isdir(os.path.abspath(os.path.join(self.currdir,self.logdir))):            
                os.mkdir(os.path.abspath(os.path.join(self.currdir,self.logdir)))
                
    def loadPostProcessorTests(self,mymacfile=None):
	"""
	Searches through the current directory looking for any files containing 
	*.mac pattern then loads the test(s) into a suite.
	"""
	suite       = unittest.TestSuite()

        loadername  = 'PostProcessorLoader'
        if self.debug:
            loadername = os.path.abspath('%s/%s/PostProcessorLoader' %(self.currdir,self.logdir))

	self.logger = TestLoggers(loadername,debugging=self.debug).getThem()

	macfiles    = []

	if (mymacfile != None):
	    self.macfile = mymacfile

	    self.ThisPostProcessorRunTest = TruchasPostProcessorRunTest(self.currdir,
                                                                        self.postprodir,
                                                                        self.macfile,
                                                                        self.BasicRunSpecs) 

	    suite.addTest( self.ThisPostProcessorRunTest)

            self.logger.debug(self.__postProcessorLoadDetails(self.macfile))
	    if self.debug: print self.ThisPostProcessorRunTest.str(self.macfile)

	else:
	    #perform search for macro files in this directory
	    files        = os.listdir(self.currdir)
	    for f in files:
		fullname = os.path.join(self.currdir, f)
		if fnmatch.fnmatch(f, self.macFilePat):

		    self.macfile                  = f
		    self.ThisPostProcessorRunTest = TruchasPostProcessorRunTest(self.currdir,
                                                                                self.postprodir,
                                                                                self.macfile,
                                                                                self.BasicRunSpecs) 

		    suite.addTest( self.ThisPostProcessorRunTest )

                    self.logger.debug(self.__postProcessorLoadDetails(self.macfile))
	            if self.debug: print self.ThisPostProcessorRunTest.str(self.macfile)

        info  = '\n'
	tmp   = 'Suite to be run          :' + str(suite)
        info += self.__rjustln(self.column2, tmp)        

        self.logger.debug(info)
        
        return suite


    def __rjustln(self,column,s):

        tmp = s.rjust(column + len(s))
        tmp += '\n'

        return tmp


    def __postProcessorLoadDetails(self, macfile):

	info  = '\n'
	tmp   = '*Loaded TruchasPostProcessorRun TestCase'
        info += self.__rjustln(self.column1, tmp)        
        tmp   = self.ThisPostProcessorRunTest.str(macfile)
        info += self.__rjustln(self.column1, tmp)        
        info += '\n'

	return info


if __name__=='__main__':

    currdir  = os.getcwd()
    toolsdir = currdir + '/../../../'
    toolsdir = os.path.abspath(toolsdir)
    truchas_bin       = testingdir + '/../../../bin'
    truchas_bin       = os.path.abspath(truchas_bin)

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
	    cmd = 'rm -r outputs*'
	    os.system(cmd)
	    cmd = 'rm -r for_debugging'
	    os.system(cmd)

	else:
            # determine output dirs (need output to do PP test)
            output_dirs = []
            for dir in os.listdir(currdir):
        	if 'outputs.t-' in dir and os.path.isdir(dir): output_dirs.append(dir)
        
            # send error msg if no output dirs
            if len(output_dirs) == 0: 
        	print >> sys.stderr, ""
        	print >> sys.stderr, "ERROR!! No outputs.t-* directories exist."
        	print >> sys.stderr, "   Run another component test that generates the directories first."
        	print >> sys.stderr, "   Then try again."
        	print >> sys.stderr, ""
            
            # loop over avaialable output dirs
            for dir in output_dirs:
        	logdir  = dir
                logfile      = logdir+'/static_drop_logs/BasicRun.log'
        
                # get corresponding exe file
                for exefile in os.listdir(truchas_bin):
        	    if exefile in logdir:
        	        vals  = string.splitfields(exefile,'.')
        	        if vals[0][2:] == string.lower(platform.system()):
                        
        	            parallel_env   = vals[3]
        	            compiler       = vals[2]
        	            compile_mode   = vals[4][:-2]
        		    np             = 1
        		    debug          = 1
        
                	    LastRunSpecs = RunTimeSpecs(logdir,logfile,parallel_env,np,compiler,compile_mode)
        
                            T = LastRunSpecs.str()
                            print T
                
                	    X       = PPLoader(currdir,debug,toolsdir,LastRunSpecs)
                	    suite   = X.loadPostProcessorTests()
                	    runner  = TestRunner(verbosity=2)
                	    result  = runner.run( suite )
            
    except:

	print "---> Test failed in some aspect <---"
	print "\nNature-Of-Error:",sys.exc_info()[0],"\n"
	if opt.d: raise
