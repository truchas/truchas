#!/usr/bin/env python
"""
CheckPPDiagnostics

-----------------------------------------------------------------------------
   Purpose:
  
      TestCase to check diagnostics created from a postprocessor macro file
      These diagnostics would be created from the query, stat and probe commands
  
   Public Interface:
  
      T = CheckPPDiagnostics(unittest.TestCase,TruchasBaseTest)
      T.setUp()
      T.shortDescription()
      T.tearDown()
      T.testDiagnosticsCreated()
      T.str()
  
   Contains:  
      class CheckPPDiagnostics
            __init__(dir,RunSpecs)
            setUp()
            shortDescription()
            tearDown()
            testDiagnosticsCreated()
            str()
  
   Unit Test Block
  
   Author(s): Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""

import unittest
import os, sys, string, re, fnmatch, platform, shutil

if __name__ == '__main__':
    print "\n for component test in %s \n" %(__file__)
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    testingdir  = os.path.abspath(thisdir + '/../')
    sys.path.append(testingdir)
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

from Runners            import TestRunner, RunTimeSpecs
from TruchasBaseTest    import TruchasBaseTest

class CheckPPDiagnostics(unittest.TestCase,TruchasBaseTest):

    def __init__(self,dir,RunSpecs):

	unittest.TestCase.__init__(self,methodName='testDiagnosticsCreated')
	self.dir          = dir
	#specifications from a previous Truchas run
	self.runspecs     = RunSpecs
	#for formatting log file
	self.column1      = 0
	self.column2      = 5	
	self.logdir       = self.runspecs.logdir
        self.logname      = 'CheckPPDiagnostics.log'
        self.basenamelog  = string.split(self.runspecs.logfile,'/')[1]
        self.logfile      = os.path.join(self.logdir,self.basenamelog,self.logname)
	self.outfile      = 'checkppdiagnostics.out'
        self.diagfilePat  = '.txt'
        self.diagfiles    = []

    def setUp(self):
	"defines diagnostic file location" 
	
	os.chdir(self.dir)

        L      = string.split(self.logdir,'/')
        logdir = L[-1]
        
        if logdir not in os.listdir(self.dir):
            os.mkdir(logdir)

        wdir = os.path.join(self.dir,self.logdir)
	os.chdir(wdir)

        if self.basenamelog not in os.listdir(wdir):
            os.mkdir(self.basenamelog)

        logfl          = os.path.join(self.basenamelog,self.logname)

	self.watcher   = open(logfl,'w')
        
	self.watcher.write('In CheckPPDiagnostics TestCase')
 	self.watcher.write('\n')
 	self.watcher.write('\n')

	tmp     = '*Setting up'
	str     = self.rjustln(self.column1, tmp)
	self.watcher.write(str)
 	self.watcher.write('\n')

    def shortDescription(self):
	"modified unittest short description to allow class name and input file name"

        classname = self.__class__.__name__
        classname = re.sub('Test[ers]*$', '', classname)
        docname   = '['+classname +']' + ' ' 
        return  docname.ljust(40) + ': ' + str(unittest.TestCase.shortDescription(self))

    def testDiagnosticsCreated(self):
	"tests diagnostics file created from a postprocessor run" 

        success       = 0
        for file in os.listdir(os.getcwd()):
            if self.diagfilePat in file:
                success       = 1
                self.diagfiles.append(file)

        tmp    = 'ERROR!! No diagnostic file created by Truchas postprocessor run'
        errstr = '\n \n' + tmp + '\n'

        if not success:
            self.watcher.write('\n')
	    self.watcher.write(errstr)
            self.watcher.write('\n')
            raise ValueError(errstr)


    def tearDown(self):
	"ensures working directory returned to original directory before test started" 

        tmp     = '*Tearing down'
        str     = self.rjustln(self.column1, tmp)
	self.watcher.write(str)

        self.watcher.write('\n') 
        self.watcher.close()

        os.chdir(self.dir)
        
        """
        if self.diagfile != None:

            tmp     = 'Moving %s to %s/ ' %(self.diagfile, self.logdir)
            str     = self.rjustln(self.column2, tmp)
            self.watcher.write(str)
            
            shutil.move(self.diagfile,self.logdir)
        """

    def str(self,column=15):
        "provides formatted info about diagnostic file created"

        info  = '\n'
	tmp   = 'Diagnostic file(s) : ' + str(self.diagfiles)
        info += self.rjustln(column, tmp)
	tmp   = 'Log dir            : ' + self.runspecs.logdir
        info += self.rjustln(column, tmp)

        return info

if __name__=='__main__':
    currdir     = os.getcwd()
    TestRunner  = unittest.TextTestRunner

    #determine runtime directories
    output_dirs = []
    for dir in os.listdir(currdir):
	if 'outputs.t-' in dir and os.path.isdir(dir): output_dirs.append(dir)
	    
    # loop over available output directories
    for dir in output_dirs:
        logdir  = dir
        logfile = logdir+'/static_drop_logs/BasicRun.log'

        # parse directory name to get runtime options
	vals = string.splitfields(dir,'.')
	if vals[1][2:] == string.lower(platform.system()):
            parallel_env = vals[4] #'serial'
            np           = 1
            compiler     = vals[3] #'lahey'
	    compile_mode = vals[5][:-2] #'opt'
            LastRunSpecs = RunTimeSpecs(logdir,logfile,parallel_env,np,compiler,compile_mode)
            
            suite        = unittest.TestSuite()
            T            = CheckPPDiagnostics(currdir,LastRunSpecs) 
            suite.addTest(T)
            runner       = TestRunner(verbosity=2)
            result       = runner.run( suite )    

            print T.str()

    



 
