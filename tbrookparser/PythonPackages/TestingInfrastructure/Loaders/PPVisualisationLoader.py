#!/usr/bin/env python
"""
 PPVisualisationLoader

-----------------------------------------------------------------------------
   Purpose:
  
      Provides methods to create TestSuites containing:
      1) CheckPPGMV     TestCase(s)
      2) CheckPPEnsight TestCase(s)
      3) CheckPPTecPlot TestCase(s)
      3) CheckPPVTK     TestCase(s)
  
   Public Interface:
  
      T = PPVisualisationLoader(currdir,RunSpecs) 
      T.loadPPGMVChecker
      T.loadPPEnsightChecker
      T.loadPPTecPlotChecker
      T.loadPPVTKChecker
  
   Contains:
      class PPVisualisationLoader(unittest.TestLoader)
             __init__(currdir,RunSpecs) 
             loadPPGMVChecker
             loadPPEnsightChecker
             loadPPTecPlotChecker
             loadPPVTKChecker
             __loadPPGeneralVisChecker(testcase,ldrname,pattername,mymacfile)
             __rjustln(columns,s)
             __PPCheckerLoadDetails(macfile)
             
      class RunThisPPVisTest
             
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
from TestCases         import CheckPPGMV, CheckPPEnsight, CheckPPTecPlot, CheckPPVTK

class PPVisualisationLoader(unittest.TestLoader):
    """
    Overriding the unittest TestLoader to allow us to search for 
    and load the following TestCase classes:

    -CheckPPGMV     (if gmv.mac file exists)
    -CheckPPEnsight (if ensight.mac file exists)
    -CheckPPTecPlot (if tecplot.mac file exists)
    -CheckPPVTK     (if vtk.mac file exists)
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

    def loadPPGMVChecker(self,pattern = 'gmv.mac'):

        suite = unittest.TestSuite()

        if self.__machine('eclair') or self.__machine('rayo'):
        #if not self.__machine('milkyway') and not self.__machine('qsc') and not self.__flash():

            suite = self.__loadPPGeneralVisChecker(testcase   = CheckPPGMV,
                                                   ldrname    = 'PPGMVCheckLoader',
                                                   pattername = pattern)

        return suite

    def loadPPEnsightChecker(self,pattern = 'ensight.mac'):

        suite = unittest.TestSuite()

        if self.__machine('eclair') or self.__machine('rayo') or self.__machine('milkyway'):
        #if not self.__machine('qsc') and not self.__flash():

            suite = self.__loadPPGeneralVisChecker(testcase   = CheckPPEnsight,
                                                   ldrname    = 'PPEnsightCheckLoader',
                                                   pattername = pattern)
        return suite

    def loadPPTecPlotChecker(self,pattern = 'tecplot.mac'):

        suite = unittest.TestSuite()

        if self.__machine('eclair') or self.__machine('rayo'):
        #if not self.__machine('milkyway') and not self.__machine('qsc') and not self.__flash():
            
            suite = self.__loadPPGeneralVisChecker(testcase   = CheckPPTecPlot,
                                                   ldrname    = 'PPTecPlotCheckLoader',
                                                   pattername = pattern)
        return suite

    def loadPPVTKChecker(self,pattern = 'vtk.mac'):

        suite = unittest.TestSuite()

        #if self.__machine('eclair') or self.__machine('rayo'):
        #if not self.__machine('milkyway') and not self.__machine('qsc') and not self.__flash():
	if os.getenv('PARAVIEW_HOME'):
            
            suite = self.__loadPPGeneralVisChecker(testcase   = CheckPPVTK,
                                                   ldrname    = 'PPVTKCheckLoader',
                                                   pattername = pattern)
        return suite

    def __loadPPGeneralVisChecker(self,
                                  testcase   = CheckPPGMV,
                                  ldrname    = 'PPGMVCheckLoader',
                                  pattername = 'gmv.mac',
                                  mymacfile  = None):
	"""
	Searches through the current directory looking for any files containing 
	'mymacfile' or 'pattername' pattern then loads the test(s) into a suite.
	"""
	suite       = unittest.TestSuite()
        loadername  = ldrname
        
        if self.debug:
            loadername  = os.path.abspath('%s/%s/%s' %(self.currdir,self.logdir,ldrname))
            
	self.logger = TestLoggers(loadername,debugging=self.debug).getThem()
        
	if (mymacfile != None):
	    self.macfile  = mymacfile

	    self.ThisTest = testcase(self.currdir,
                                     self.RunSpecs) 

	    suite.addTest(self.ThisTest)

            self.logger.debug(self.__PPCheckerLoadDetails())
	    if self.debug: print self.ThisTest.str(self.macfile)


	else:
	    #perform search for macro files in this directory
	    files           = os.listdir(self.currdir)

	    for f in files:
		fullname = os.path.join(self.currdir, f)
		if fnmatch.fnmatch(f, pattername):

		    self.macfile = f
		    self.ThisTest = testcase(self.currdir,
                                             self.RunSpecs) 
		    suite.addTest( self.ThisTest )

                    self.logger.debug(self.__PPCheckerLoadDetails())
	            if self.debug:
                        print self.ThisTest.str()

        info  = '\n'
	tmp   = 'Suite to be run          :' + str(suite)
        info += self.__rjustln(self.column2, tmp)        

        self.logger.debug(info)
        
        return suite

    def __machine(self,machine):

        try:
            check = machine in os.environ['HOST']
        except EnvironmentError:
            check = machine in os.environ['HOSTNAME']
        except:
            check = 0

        return check

    def __flash(self):

        try:
            check = 'flash' in os.environ['HOSTNAME'] or 'ffe' in os.environ['HOSTNAME']
        except EnvironmentError:
            check = 'flash' in os.environ['HOST'] or 'ffe' in os.environ['HOST']    
        except:
            check = 0
            
        return check

    def __rjustln(self,column,s):

        tmp = s.rjust(column + len(s))
        tmp += '\n'

        return tmp

    def __PPCheckerLoadDetails(self):

	info  = '\n'
	tmp   = '*Loaded CheckPPVisualisation TestCase'
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
            from Runners import TestRunner, RunTimeSpecs         
            suite        = unittest.TestSuite()
        	
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
        
           	        X       = PPVisualisationLoader(currdir,debug,LastRunSpecs)
                	suite2  = X.loadPPGMVChecker()
                        suite.addTest(suite2)
                        suite3  = X.loadPPTecPlotChecker()
                        suite.addTest(suite3)
                        suite4  = X.loadPPEnsightChecker()
                        suite.addTest(suite4)
                        suite5  = X.loadPPVTKChecker()
                        suite.addTest(suite5)
                
                	runner  = TestRunner(verbosity=2)
                	result  = runner.run( suite )
      
    except:

	print "---> Test failed in some aspect <---"
	print "\nNature-Of-Error:",sys.exc_info()[0],"\n"
	if opt.d: raise
 
