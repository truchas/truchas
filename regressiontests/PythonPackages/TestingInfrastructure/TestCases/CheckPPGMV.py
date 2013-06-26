#!/usr/bin/env python
"""
CheckPPGMV

-----------------------------------------------------------------------------
   Purpose:
  
      TestCase to check validity of GMV file (gmv.bin, gmv.ascii) created from
      gmv.mac file 
  
   Public Interface:
  
      T = CheckPPGMV(unittest.TestCase,TruchasBaseTest)
      T.setUp()
      T.shortDescription()
      T.testGMVFile()
      T.tearDown()
      T.str()
  
   Contains:  
      class CheckPPGMV
            __init__(dir,RunSpecs)
            setUp()
            shortDescription()
            testGMVFile()
            tearDown()
            str()
  
   Unit Test Block
  
   Author(s): Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""

import unittest
import os, sys, string, re, fnmatch, shutil, platform

if __name__ == '__main__':
    print "\n for component test in %s \n" %(__file__)
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    testingdir  = os.path.abspath(thisdir + '/../')
    sys.path.append(testingdir)
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

from Runners            import TestRunner, RunTimeSpecs
from TruchasBaseTest    import TruchasBaseTest

class CheckPPGMV(unittest.TestCase,TruchasBaseTest):

    def __init__(self,dir,RunSpecs):

	unittest.TestCase.__init__(self,methodName='testGMVFile')
	self.dir          = dir
	#specifications from a previous Truchas run
	self.runspecs     = RunSpecs        
	#for formatting log file
	self.column1      = 0
	self.column2      = 5	
	self.logdir       = self.runspecs.logdir
        self.logname      = 'CheckPPGMV.log'
        self.basenamelog  = string.split(self.runspecs.logfile,'/')[1]
	self.logfile      = os.path.join(self.logdir,self.basenamelog,self.logname)
	self.outfile      = 'checkppgmv.out'
        self.gmvbinfile   = 'gmv*.bin'
        self.gmvasciifile = 'gmv*.ascii'
        self.gmvfiles     = []
        self.curfile      = None

    def setUp(self):
	"defines gmv file location" 
	
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
        
	self.watcher.write('In CheckPPGMV TestCase')
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

    def testGMVFile(self):
	"tests validity of GMV file created from a postprocessor run or from self.plotter"

	# write to watcher
	tmp     = '*testGMVFile'
	str     = self.rjustln(self.column1, tmp)
	self.watcher.write(str)
 	self.watcher.write('\n')

	# setup error status and message
        success = 0
        tmp     = 'ERROR!! No GMV file created by Truchas postprocessor run or from self.plotter'
        errstr  = '\n \n' + tmp + '\n'
	err     = ''

        for file in os.listdir(os.getcwd()):
            if fnmatch.fnmatch(file,self.gmvbinfile) or fnmatch.fnmatch(file,self.gmvasciifile):
                
                self.curfile = file
                self.gmvfiles.append(file)
	        self.watcher.write(self.str(self.column2))

		cmd = 'gmvbatch -i %s -s %s.rgb 1>stdout 2>stderr' % (file,file)
		os.system(cmd)

		# test for success
		rgbout = file + '.rgb'
		for rgbfile in os.listdir(os.getcwd()):
		    if rgbout in rgbfile:
                        success = 1
		        if rgbfile not in self.gmvfiles:
                            self.gmvfiles.append(rgbfile) 
		        tmp  = 'Output File     : '
			tmp += rgbfile
		        str  = self.rjustln(self.column2, tmp)
                        self.watcher.write(str)
                        
		# write result to watcher
		if success: 
		    tmp = 'Status          : PASSED'
		    str = self.rjustln(self.column2, tmp)
                    self.watcher.write(str)
 	            self.watcher.write('\n')
		else:       
		    tmp    = 'Status          : FAILED -> STDERR '
		    str    = self.rjustln(self.column2, tmp)
                    self.watcher.write(str)
		    # if failed, append stderr contents
		    stderr = open('stderr','r')
		    tmp    = stderr.readlines()
                    self.watcher.writelines(tmp)
		    stderr.close()
		    err   += 'ERROR!!  Failed to generate RGB file from %s \n' % file

        if not success:
            if len(err) > 0:
                errstr  = '\n \n' + err + '\n'
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

    def str(self,column=15):
        "provides formatted info about gmv file created"

        info  = '\n'
	tmp   = 'GMV file(s)     : ' + str(self.gmvfiles)
        info += self.rjustln(column, tmp)
	tmp   = 'Log dir         : ' + self.runspecs.logdir
        info += self.rjustln(column, tmp)

        return info


if __name__=='__main__':
    currdir     = os.getcwd()
    TestRunner  = unittest.TextTestRunner    

    #for testing purposes specify dummy runtime specs 
    for dir in os.listdir(currdir):
	if os.path.isdir(dir) and 'outputs' in dir:
            logdir  = dir
            logfile = logdir+'/static_drop_logs/BasicRun.log'

	    # parse directory name for run time specs
	    vals = string.splitfields(dir,'.')
	    if vals[1][2:] == string.lower(platform.system()):
	        parallel_env = vals[4]
	        np           = 1
	        compiler     = vals[3]
	        compile_mode = vals[5][:-2]
                LastRunSpecs = RunTimeSpecs(logdir,logfile,parallel_env,np,compiler,compile_mode)

                suite        = unittest.TestSuite()
                T            = CheckPPGMV(currdir,LastRunSpecs) 
                suite.addTest(T)
                runner       = TestRunner(verbosity=2)
                result       = runner.run( suite )    
                print T.str()




 
