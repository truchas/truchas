#!/usr/bin/env python
"""
CheckPPEnsight

-----------------------------------------------------------------------------
   Purpose:
  
      TestCase to check validity of EnSight file (*.CASE) created from
      ensight.mac file 
  
   Public Interface:
  
      T = CheckPPEnsight(unittest.TestCase,TruchasBaseTest)
      T.setUp()
      T.shortDescription()
      T.testEnSightFile()
      T.tearDown()
      T.str()
  
   Contains:  
      class CheckPPEnsight
            __init__(dir,RunSpecs)
            setUp()
            shortDescription()
            testEnSightFile()
            tearDown()
            str()
  
   Unit Test Block
  
   Author(s): Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""

import unittest
import os, sys, string, re, fnmatch, commands, platform

if __name__ == '__main__':
    print "\n for component test in %s \n" %(__file__)
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    testingdir  = os.path.abspath(thisdir + '/../')
    sys.path.append(testingdir)
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

from Runners            import TestRunner, RunTimeSpecs
from TruchasBaseTest    import TruchasBaseTest

class CheckPPEnsight(unittest.TestCase,TruchasBaseTest):

    def __init__(self,dir,RunSpecs):

	unittest.TestCase.__init__(self,methodName='testEnSightFile')
	self.dir                = dir
	#specifications from a previous Truchas run
	self.runspecs           = RunSpecs        
	#for formatting log file
	self.column1            = 0
	self.column2            = 5	
	self.logdir             = self.runspecs.logdir
        self.logname            = 'CheckPPEnsight.log'

        self.basenamelog        = None
        errstr                  = '/n/n In %s logfile is %s which does not have the correct form "_logs" /n/n ' %(__file__,self.runspecs.logfile)
        for word in string.split(self.runspecs.logfile,'/'): 
            if '_logs' in word:
                self.basenamelog = word
        assert self.basenamelog != None, errstr
        
        self.basename           = string.split(self.basenamelog,'_logs')[0]
	self.logfile            = os.path.join(self.logdir,self.basenamelog,self.logname) 
	self.outfile            = 'checkppensight.out'
        self.ensightcasefilePat = '.ensight.CASE'
	self.successlines       = ['Data verification SUCCESSFUL', 'Data format verification SUCCESSFUL'] 
	self.failurelines       = ['Verification of the data FAILED']
        self.ensightfile        = None

    def setUp(self):
	"defines ensight file location and sets the working dir" 
	
	os.chdir(self.dir)

        L      = string.split(self.logdir,'/')
        logdir = L[-1]
        
        if logdir not in os.listdir(self.dir):
            os.mkdir(logdir)

        wdir = os.path.join(self.dir,self.logdir)
	os.chdir(wdir)

        if self.basenamelog not in os.listdir(wdir):
            os.mkdir(self.basenamelog)
        
        logfl         = os.path.join(self.basenamelog,self.logname)

        self.watcher  = open(logfl,'w')
        
	self.watcher.write('In CheckPPEnsight TestCase')
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

    def testEnSightFile(self):
	"tests validity of EnSight case file created from a postprocessor run" 

        # write to watcher
        tmp = '*testEnSightFile'
        str = self.rjustln(self.column1, tmp)
        self.watcher.write(str)
        self.watcher.write('\n')

        # setup error status message
        success = 0
        tmp     = 'ERROR!! No EnSight case file created by Truchas postprocessor run'
        errstr  = '\n \n' + tmp + '\n'
        
	outputdir = os.path.join(os.getcwd(), self.basename + '_output')	
	 
	os.chdir(outputdir)
        for file in os.listdir(os.getcwd()):
            if self.ensightcasefilePat in file:
                self.ensightfile = file
                self.watcher.write(self.str(self.column2))
                
		# the os.system execution didn't work on milkyway,
		# could not deal with the '>&' redirect required by ens_checker
		# switched to commands instead
                cmd    = 'ens_checker ' + self.basename + '.ensight.CASE '
		output = commands.getoutput(cmd)

		fd     = open(self.outfile,'w')
		fd.write(output)
		fd.close()
                
                tmp = 'Output File     : ' + self.outfile
                str = self.rjustln(self.column2, tmp)
                self.watcher.write(str)
                outfl = os.path.join(os.getcwd(),self.outfile)

                # parse the ens_checker output
		stdout = open(self.outfile,'r')
		lines  = stdout.readlines()
		for line in lines:
		    for a in self.failurelines:
			if a in line:
			   success = 0
                           errstr  = '\n In %s \n' %(outfl)
                           errstr += '\n \n Ens_Checker found a problem: \n' + line #+ '\n'	
		    for b in self.successlines:
			if b in line:	
			   success = 1	
		    tmp = 'Status          : ' + line[:-1]
		    str = self.rjustln(self.column2, tmp)
		    self.watcher.write(str)
		stdout.close()

                # write to watcher and append output to logfile
		str = self.rjustln(self.column2,'Output from Ens_Checker: \n')
		self.watcher.write(str)
                
		stdout = open(self.outfile,'r')
                tmp    = stdout.readlines()
		for line in tmp:
                    self.watcher.write(self.rjustln(self.column2,'     '+line[:-1]))
                self.watcher.write('\n')
		stdout.close()

	os.chdir(self.dir)

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
        tmp     = 'Moving %s to %s/ ' %(self.outfile, self.logdir)
        str     = self.rjustln(self.column2, tmp)
        self.watcher.write(str)

        try:
            shutil.copytree(self.basename + '_output/',self.logdir)
        except:
            cmd = 'cp -r %s %s/' %(self.basename + '_output/', self.logdir)
            os.system(cmd)

        cmd     = 'mv %s %s/' %(self.logdir + '/'+ self.basename + '_output/' + self.outfile, self.logdir)
        os.system(cmd)
        """

    def str(self,column=15):
        "provides formatted info about ensight file created"

        info  = '\n'
	tmp   = 'EnSight file    : ' + str(self.ensightfile)
        info += self.rjustln(column, tmp)

        outputdir = self.basename + '_output'
        tmp   = 'Output dir      : ' + outputdir
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
                T            = CheckPPEnsight(currdir,LastRunSpecs) 
                suite.addTest(T)
                runner       = TestRunner(verbosity=2)
                result       = runner.run( suite )    
                print T.str()




 
