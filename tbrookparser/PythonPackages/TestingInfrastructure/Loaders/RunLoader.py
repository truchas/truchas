#!/usr/bin/env python
"""
 RunLoader

-----------------------------------------------------------------------------
   Purpose:
  
      Provides methods to create TestSuites containing:
      1) Truchas BasicRun TestCase(s)
      2) Truchas BasicRunFail TestCase(s)
      3) Truchas RestartRun TestCase(s)
  
   Public Interface:
  
      T = RunLoader(dir,parallel_env,np,compiler,compile_mode,
                   truchas_bin,truchas_exe)
      T.loadBasicRunTests(inpfile)
      T.loadBasicRunFailTests(inpfile)
      T.getBasicRunSpecs()
      T.loadRestartRunTests(inpfile,restartfile)
  
   Contains:
      class RunLoader(unittest.TestLoader)
            __init__(dir,datastores,restartdatastores,
                     parallel_env,np,compiler,compile_mode,
                     truchas_bin,truchas_exe)
            loadBasicRunTests(inpfile)
            loadBasicRunFailTests(inpfile)
            getBasicRunSpecs()
            loadRestartRunTests(inpfile,restartfile)
            __constructLogDirName()
            __rjustln(column,s)
            __basicRunLoadDetails(inpfile)
            __restartRunLoadDetails(inpfile, restrtfile)

      class RunTheseTruchasTests
            __init__(dir,parallel_env,np,compiler,compile_mode,
                           truchas_bin)

   Unit Test Block
  
   Author: Sharen Cummins (scummins@lanl.gov)
           Erin Iesulauro Barker (eibarker@lanl.gov)
-----------------------------------------------------------------------------
"""

import unittest
import os, sys, fnmatch, string, platform

if __name__=='__main__':
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    testingdir  = thisdir + '/../'
    sys.path.append(testingdir)
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)
    from Runners       import TestRunner

from Runners           import RunTimeSpecs 
from Loggers           import TestLoggers
from TestCases         import TruchasBasicRunTest, TruchasRestartRunTest, TruchasBaseTest, TruchasBasicRunFailTest

class RunLoader(unittest.TestLoader,TruchasBaseTest):
    """
    Overriding the unittest loader to allow us to search for 
    and load the following TestCase classes:
    TruchasBasicRunTest   (if *.inp file exists)
    TruchasRestartRunTest (if restart.mac and *.inp file exists).

    """

    def __init__(self,
                 datastores,
                 restartdatastores,
                 currdir      = '.',
		 debug        = 0,
		 parallel_env = 'serial',
		 np           = 1,
		 compiler     = 'lahey',
		 compile_mode = 'debug',
		 version      = '2.3',
		 truchas_bin  = '.',
		 truchas_exe  = None):

        self.datastores        = datastores
        self.restartdatastores = restartdatastores
	self.debug             = debug
	self.currdir           = currdir
	self.parallel_env      = parallel_env
	self.np                = np
	self.compiler          = compiler
	self.compile_mode      = compile_mode
	self.version           = version
	self.truchas_bin       = truchas_bin
	self.truchas_exe       = truchas_exe

	self.inpFilePat        = '*.inp'
	self.failFilePat       = '*_fail.inp'
	self.restartFilePat    = '*restart.mac'
        self.restartfile       = 'restart.bin'
	self.inpRestartFilePat = '*_restart.inp'

	self.logdir          = 'for_debugging'
	#for formatting debug file
	self.column1         = 5
	self.column2         = 10


        if self.debug:
            if not os.path.isdir(os.path.abspath(os.path.join(self.currdir,self.logdir))):            
                os.mkdir(os.path.abspath(os.path.join(self.currdir,self.logdir)))

        if self.truchas_exe == None:
            self.truchas_exe = self.constructTruchasExe()
            
        runspecslogdir  = self.constructLogDirName()
        runspecslogfile = ''

	self.runspecs   = RunTimeSpecs(runspecslogdir,
                                       runspecslogfile,
				       self.parallel_env,
				       self.np,
				       self.compiler,
				       self.compile_mode,
                                       self.version)           

    def loadBasicRunTests(self,myinputfile=None):
        """
	For the currdir use myinputfile to perform the TruchasBasicRun test.
        If myinputfile=None, search for any input files in the currdir, if 
        they exist then perform the TruchasBasicRun test
	"""

        suite       = unittest.TestSuite()
        loadername  = 'BasicRunLoader'

        if self.debug:
            loadername = os.path.abspath('%s/%s/BasicRunLoader' %(self.currdir,self.logdir))
            
	self.logger = TestLoggers(loadername,debugging=self.debug).getThem()

	if (myinputfile != None):
	    
	    self.ThisBasicRunTest = TruchasBasicRunTest(self.currdir,
							myinputfile,
							self.parallel_env,
							self.np,
							self.compiler,
							self.compile_mode,
							self.version,
							self.truchas_bin,
							self.truchas_exe) 

	    
	    suite.addTest( self.ThisBasicRunTest)
            self.datastores.append(self.ThisBasicRunTest.datastore)

            self.logger.debug(self.__basicRunLoadDetails(myinputfile))

            self.runspecs.logfile = self.ThisBasicRunTest.logfile
	    if self.debug: print self.ThisBasicRunTest.str(myinputfile)

	else:
	    #perform search for input files in this directory
	    dir   = [self.currdir]
	    files = os.listdir(self.currdir)
	
	    for f in files:
		fullname = os.path.join(self.currdir, f)
		if fnmatch.fnmatch(f, self.inpFilePat) and not('restart' in f) and not ('_fail' in f):

		    self.inpfile          = f
		    self.ThisBasicRunTest = TruchasBasicRunTest(self.currdir,
                                                                self.inpfile,
                                                                self.parallel_env,
                                                                self.np,
                                                                self.compiler,
                                                                self.compile_mode,
								self.version,
                                                                self.truchas_bin,
                                                                self.truchas_exe)


		    suite.addTest( self.ThisBasicRunTest)
                    self.datastores.append(self.ThisBasicRunTest.datastore)

                    self.logger.debug(self.__basicRunLoadDetails(self.inpfile))

                    self.runspecs.logfile = self.ThisBasicRunTest.logfile
	            if self.debug: print self.ThisBasicRunTest.str(self.inpfile)

        info  = '\n'
        tmp   = 'Suite to be run         : ' + str(suite)
        info += self.__rjustln(self.column2, tmp)
        
        self.logger.debug(info)


	return suite, self.runspecs.logdir 

    def loadBasicRunFailTests(self,myinputfile=None):
        """
	For the currdir use myinputfile to perform the TruchasBasicRunFail test.
        If myinputfile=None, search for any input files in the currdir, if 
        they exist then perform the TruchasBasicRunFail test
	"""

        suite       = unittest.TestSuite()

        loadername  = 'BasicRunFailLoader'
        if self.debug:
            loadername = os.path.abspath('%s/%s/BasicRunFailLoader' %(self.currdir,self.logdir))
        
	self.logger = TestLoggers(loadername,debugging=self.debug).getThem()

	if (myinputfile != None):
	    
	    self.ThisBasicRunFailTest = TruchasBasicRunFailTest(self.currdir,
						   	        myinputfile,
							        self.parallel_env,
							        self.np,
							        self.compiler,
							        self.compile_mode,
								self.version,
							        self.truchas_bin,
							        self.truchas_exe) 

	    
	    suite.addTest( self.ThisBasicRunFailTest)
            self.datastores.append(self.ThisBasicRunFailTest.datastore)

            self.logger.debug(self.__basicRunFailLoadDetails(myinputfile))

            self.runspecs.logfile = self.ThisBasicRunFailTest.logfile
	    if self.debug: print self.ThisBasicRunFailTest.str(myinputfile)

	else:
	    #perform search for input files in this directory
	    dir   = [self.currdir]
	    files = os.listdir(self.currdir)
	
	    for f in files:
		fullname = os.path.join(self.currdir, f)
		if fnmatch.fnmatch(f, self.failFilePat):

		    self.inpfile              = f
		    self.ThisBasicRunFailTest = TruchasBasicRunFailTest(self.currdir, 
				                                        self.inpfile, 
									self.parallel_env,
									self.np, 
									self.compiler, 
									self.compile_mode,
									self.version,
									self.truchas_bin, 
									self.truchas_exe)

		    suite.addTest( self.ThisBasicRunFailTest)
                    self.datastores.append(self.ThisBasicRunFailTest.datastore)

                    self.logger.debug(self.__basicRunFailLoadDetails(self.inpfile))

                    self.runspecs.logfile = self.ThisBasicRunFailTest.logfile
	            if self.debug: print self.ThisBasicRunFailTest.str(self.inpfile)

        info  = '\n'
        tmp   = 'Suite to be run         : ' + str(suite)
        info += self.__rjustln(self.column2, tmp)
        
        self.logger.debug(info)


	return suite, self.runspecs.logdir 

    def loadRestartRunTests(self,myinputfile=None,myrestartfile=None):
        """
	For the currdir use myinputfile to perform the TruchasRestartRun test.
        If myinputfile=None search for any restart.mac files, if they exist then 
	perform the TruchasRestartRun test with existing Truchas input file.
	"""

        suite        = unittest.TestSuite()

        loadername   = 'RestartRunLoader'

        if self.debug:
            loadername   = os.path.abspath('%s/%s/RestartRunLoader' %(self.currdir,self.logdir))
            
	self.logger  = TestLoggers(loadername,debugging=self.debug).getThem()


	#if (myinputfile != None) and (myrestartfile != None):
	if (myinputfile != None) or (myrestartfile != None):

	    if myinputfile   == None: myinputfile   = self.inpfile
	    if myrestartfile == None: myrestartfile = self.restartfile

	    self.ThisRestartRunTest = TruchasRestartRunTest(self.currdir,
                                                            myinputfile,
                                                            myrestartfile,
                                                            self.runspecs,
                                                            self.parallel_env,
                                                            self.np,
                                                            self.compiler,
                                                            self.compile_mode,
                                                            self.truchas_bin,
                                                            self.truchas_exe) 

            self.logger.debug(self.__restartRunLoadDetails(myinputfile,myrestartfile))
            self.restartdatastores.append(self.ThisRestartRunTest.datastore)

	    suite.addTest( self.ThisRestartRunTest)
	    if self.debug: print self.ThisRestartRunTest.str(myinputfile,myrestartfile)

	else:  
        
	    dir     = [self.currdir]
	    files   = os.listdir(self.currdir)
	
	    restart   = 0
            pprestart = 0
	    for f in files:
		fullname = os.path.join(self.currdir, f)
                if fnmatch.fnmatch(f, self.restartFilePat):
		    #restart.mac file exists in the directory...now check for input file
		    pprestart = 1
                    macrofile = os.path.abspath(self.currdir+'/'+f)
                if fnmatch.fnmatch(f,self.restartfile):
                    #restart.bin already exists in directory
                    restart   = 1
	    if pprestart:
		errstring        = '\n \n'
		errstring       += '      "%s" not in %s \n' %(self.restartfile,macrofile)
		errstring       += '\n    Ensure name of your restart file in your mac file is %s \n' %(self.restartfile)
                macfile          = open(macrofile,'r')
                S                = macfile.read()
                macfile.close()
		assert self.restartfile in S, errstring

            if pprestart or restart:
		#first check for *.inp file and use this as the first choice for input file
		for f2 in files:
		    if fnmatch.fnmatch(f2,self.inpFilePat) and not ('_fail' in f2):
			self.inpfile = f2
		#now check for *_restart.inp file, if it exists use this as the input file
		for f2 in files:
		    if fnmatch.fnmatch(f2,self.inpRestartFilePat):
			self.inpfile = f2

		errstring = '\n \n'
		errstring += '           No input file exists in %s. \n' %(self.currdir)
		errstring += '\n         Need input files of the form *.inp or *_restart.inp \n'
		assert len(self.inpfile) > 0, errstring

		self.ThisRestartRunTest  = TruchasRestartRunTest(self.currdir,
								 self.inpfile,
								 self.restartfile,
								 self.runspecs,
								 self.parallel_env,
								 self.np,
								 self.compiler,
								 self.compile_mode,
								 self.truchas_bin,
								 self.truchas_exe) 
								 
		suite.addTest( self.ThisRestartRunTest)
                self.restartdatastores.append(self.ThisRestartRunTest.datastore)

                self.logger.debug(self.__restartRunLoadDetails(self.inpfile,
			                                       self.restartfile))
	        if self.debug: print self.ThisRestartRunTest.str(self.inpfile,self.restartfile)

        info  = '\n'
        tmp   = 'Suite to be run         : ' + str(suite)
        info += self.__rjustln(self.column2, tmp)
        
        self.logger.debug(info)


	return suite, self.runspecs.logdir


    def __rjustln(self,column,s):

        tmp = s.rjust(column + len(s))
        tmp += '\n'

        return tmp

    def __basicRunLoadDetails(self, inpfile):

	info  = '\n'
	tmp   = '*Loaded TruchasBasicRun TestCase: ' + str(inpfile)
        info += self.__rjustln(self.column1, tmp)        
        tmp   = self.ThisBasicRunTest.str(inpfile)
        info += self.__rjustln(self.column1, tmp)        
        info += '\n'

	return info
	
    def __basicRunFailLoadDetails(self, inpfile):

	info  = '\n'
	tmp   = '*Loaded TruchasBasicRunFail TestCase: ' + str(inpfile)
        info += self.__rjustln(self.column1, tmp)        
        tmp   = self.ThisBasicRunFailTest.str(inpfile)
        info += self.__rjustln(self.column1, tmp)        
        info += '\n'

	return info
	

    def __restartRunLoadDetails(self, inpfile, restrtfile):

	info  = '\n'
	tmp   = '*Loaded TruchasRestartRun TestCase: ' + str(inpfile)
        info += self.__rjustln(self.column1, tmp)        
        tmp   = self.ThisRestartRunTest.str(inpfile, restrtfile)
        info += self.__rjustln(self.column1, tmp)        
        strng = tmp.rjust(self.column2 + len(tmp))
	info += '\n'

	return info
	

if __name__=='__main__':

    # use .inp files from TestCases
    # start by removing any .inp files in the curdir
    files = os.listdir('.')
    for f in files: 
	if fnmatch.fnmatch(f,'*.inp'): os.remove(f)
    # now copy .inp files from ../TestCases
    files = os.listdir('../TestCases')
    for f in files: 
	if fnmatch.fnmatch(f,'*.inp'): os.link('../TestCases/'+f,f)

    currdir           = os.getcwd()
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
        		version        = vals[4][len(vals[4])-1] + '.' + vals[5]
        		np             = 1
        		debug          = 1
        
                        X              = RunLoader(datastores,
        				           restartdatastores,
        					   currdir,
        					   debug,
        					   parallel_env,
        					   np,
        					   compiler,
        					   compile_mode,
        					   version,
        					   truchas_bin)
        		
                        suite1, logdir = X.loadBasicRunTests()
                        suite2, logdir = X.loadBasicRunFailTests()
                        suite3, logdir = X.loadRestartRunTests()
        
                        suite    = unittest.TestSuite()
                        suite.addTest(suite1)
                        suite.addTest(suite2)
        
                        fp       = sys.stdout
                        runner   = TestRunner(stream=fp,verbosity=2)
        
                        if X.restartfile in os.listdir(currdir):
                            #restart file exists ...
        		    #    so can run both Truchas BasicRun and Truchas RestartRun 
                            suite.addTest(suite3)
        		    
                        result   = runner.run( suite )
        
    except:

	print "---> Test failed in some aspect <---"
	print "\nNature-Of-Error:",sys.exc_info()[0],"\n"
	if opt.d: raise
