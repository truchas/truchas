#!/usr/bin/env python

"""
TruchasBasicRunTest

-----------------------------------------------------------------------------
   Purpose:
  
      TestCase to run a basic Truchas simulation
  
   Public Interface:
  
      T = TruchasBasicRunTest(unittest.TestCase,TruchasBaseTest)
      T.setUp()
      T.shortDescription()
      T.tearDown()
      T.testBasicRun()
      T.str()
  
   Contains:  
      class TruchasBasicRunTest
            __init__(dir,inpfile,parallel_env,
                     np,compiler,compile_mode,
                     truchas_bin,truchas_exe)
            setUp()
            shortDescription()
            tearDown()
            testBasicRun()
            str()
  
   Unit Test Block
  
   Author(s): Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""

import unittest
import os, sys, string, re, string, platform, shutil

if __name__=='__main__':
    print "\n for component test in %s \n" %(__file__)
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    testingdir  = thisdir + '/../'
    sys.path.append(testingdir)
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

from ProgramEnvironment import createRunTimeProgram
from Runners            import TestRunner
from TruchasBaseTest    import TruchasBaseTest
from DataHandlers       import getDataCreator

class TruchasBasicRunTest(unittest.TestCase,TruchasBaseTest):

    def __init__(self,dir,
		 inpfile,
		 parallel_env  = 'serial',
		 np            = 1,
		 compiler      = 'lahey',
		 compile_mode  = 'debug',
		 version       = '2.3',
		 truchas_bin   = '.',
		 truchas_exe   = None):
	
	unittest.TestCase.__init__(self,methodName='testBasicRun')
	self.inpfile      = inpfile
	self.dir          = dir
	self.parallel_env = parallel_env
	self.np           = np
	self.compile_mode = compile_mode
	self.compiler     = compiler
	self.version      = version
	self.truchas_bin  = truchas_bin
	self.truchas_exe  = truchas_exe

        #string to diagnose success of this Truchas run...
        self.success      = 'truchas terminated normally' 
	#for formatting log file
	self.column1      = 0
	self.column2      = 5
	self.column3      = 10

	#construct executable name
	if self.truchas_exe == None:
	    self.truchas_exe = self.constructTruchasExe()
        #given executable, extract parallel_env, compile_mode, compiler  
        else:
            self.deconstructTruchasExe()

	self.progname = '%s/%s' %(self.truchas_bin,self.truchas_exe)

	iend          = string.rfind(self.inpfile,".")
	self.basename = self.inpfile[0:iend]

	#construct log dir name and output files
	self.logdir      = self.constructLogDirName()
	self.logname     = 'BasicRun.log'
        self.basenamelog = self.basename +'_logs'
        self.logfile     = os.path.join(self.logdir,self.basenamelog,self.logname)
	self.outdir      = self.basename+'_output' 
	self.outfile     = self.basename+'.brout'

        dsdir            = os.path.join(self.dir,self.logdir,self.outdir)
        self.datastore   = getDataCreator(dir=dsdir)

    def setUp(self):
	"defines default runtime environment, input file name, working dir location" 

        os.chdir(self.dir)

        if self.logdir not in os.listdir(self.dir):
            os.mkdir(self.logdir)

        wdir = os.path.join(self.dir,self.logdir)
	os.chdir(wdir)

        if self.basenamelog not in os.listdir(wdir):
            os.mkdir(self.basenamelog)

        logfl          = os.path.join(self.basenamelog,self.logname)

	self.watcher   = open(logfl,'w')
	self.watcher.write('In BasicRun TestCase')
 	self.watcher.write('\n')
 	self.watcher.write('\n')

	tmp     = '*Setting up'
        str     = self.rjustln(self.column1, tmp)
	self.watcher.write(str)

	self.watcher.write(self.str(self.inpfile,self.column2))

    def shortDescription(self):
	"modified unittest short description to allow class name and input file name"

        classname = self.__class__.__name__
        classname = re.sub('Test[ers]*$', '', classname)
	inputname = self.inpfile
        docname   = '['+classname +']' + ' ' + inputname
        return  docname.ljust(40) + ': ' + str(unittest.TestCase.shortDescription(self))

    def tearDown(self):
	"ensures working directory returned to original directory before test started"

        self.datastore.dir = os.path.join(self.logdir,self.outdir)
	
	tmp     = '*Tearing down'
        str     = self.rjustln(self.column1, tmp)
	self.watcher.write(str)

        self.watcher.write('\n')
        self.watcher.close()

        os.chdir(self.dir)
        
    def testBasicRun(self):
	"runs a basic Truchas simulation" 
	"(i.e runs truchas -o:output_dir input_file)"

	tmp     = '*testMethod = testBasicRun'
        str     = self.rjustln(self.column1, tmp)
	self.watcher.write(str)
        self.watcher.write('\n')

	# runs and returns good output code
        inpfl  = '../%s' %(self.inpfile)
        
        if self.parallel_env == 'parallel':
            args = ['-np %d '%(self.np), '-t', '-o:%s'%(self.outdir), inpfl]
        else:
            args = ['-t', '-o:%s'%(self.outdir), inpfl]

	T        = createRunTimeProgram(self.progname, args, self.outfile)

	retcode = T.run()
        output  = file(self.outfile,'r').read()
	self.watcher.write(T.str())

	if retcode!=0:
            outfl = os.path.join(os.getcwd(),self.outfile)
            tmp   = 'ERROR!! Program %s returned code of %d.' %(self.truchas_exe,retcode)
            tmp  += '\nSee %s for details!' %(outfl)
            str   = self.rjustln(self.column2, tmp)
            self.watcher.write('\n')
	    self.watcher.write(str)
            self.watcher.write('\n') 

            errstr = '\n \n' + tmp + '\n'
	    raise ValueError(errstr)

        if not (self.success in output) :
            outfl = os.path.join(os.getcwd(),self.outfile)
	    tmp   = 'ERROR!! Program %s terminated prematurely. \nSee %s for details!' %(self.truchas_exe,outfl)
            str   = self.rjustln(self.column2, tmp)
            self.watcher.write('\n')
	    self.watcher.write(str)
            self.watcher.write('\n')

            errstr = '\n \n' + tmp + '\n'
	    raise ValueError(errstr)

        #now write out statistics associated with this run....

        tmp       = '%s_OUTPUT STATISTICS:' %(self.basename.upper())
        str  = self.rjustln(self.column2, tmp)
        self.watcher.write(str)
        stats     = self.createStatistics()
        str       = self.rjustln(self.column2,stats)
        self.watcher.write(str)
	self.watcher.write('\n')

    def str(self,inpfile,column=15):
	"provides formatted info about RunTime specs for a given inpfile"

	info  = '\n'
	tmp   = 'Input file         : ' + inpfile
        info += self.rjustln(column, tmp)
	tmp   = 'Output file        : ' + self.outfile
        info += self.rjustln(column, tmp)
	tmp   = 'Current dir        : ' + self.dir
        info += self.rjustln(column, tmp)
	tmp   = 'Output dir         : ' + self.outdir
        info += self.rjustln(column, tmp)
	tmp   = 'Parallel env       : ' + self.parallel_env
        info += self.rjustln(column, tmp)
	if self.parallel_env == 'parallel':
	    tmp   = 'np                 : ' + str(self.np)
            info += self.rjustln(column, tmp)
	tmp   = 'Compiler           : ' + self.compiler
        info += self.rjustln(column, tmp)
	tmp   = 'Compile mode       : ' + self.compile_mode
        info += self.rjustln(column, tmp)
	tmp   = 'Truchas bin dir    : ' + self.truchas_bin
        info += self.rjustln(column, tmp)
	tmp   = 'Truchas exe        : ' + str(self.truchas_exe)
        info += self.rjustln(column, tmp)
	tmp   = 'Log file           : ' + self.logfile 
        info += self.rjustln(column, tmp)
	info += '\n'

	return info

if __name__=='__main__':
    currdir     = os.getcwd()
    truchas_bin = testingdir + '/../../../bin'
    truchas_bin = os.path.abspath(truchas_bin)

    TestRunner  = unittest.TextTestRunner
    runner      = TestRunner(verbosity=2)

    from PYTHONutils import uTestOpts
    import Clean
    
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
	    Clean.Clean()

	else:
            for exefile in os.listdir(truchas_bin):
         	if os.path.isfile(os.path.join(truchas_bin,exefile)) and exefile != 'truchas':
        	    vals  = string.splitfields(exefile,'.')
        	    if vals[0][2:] == string.lower(platform.system()):
                        suite        = unittest.TestSuite()
                        
        	        parallel_env = vals[3]
        	        compiler     = vals[2]
        	        compile_mode = vals[4][:-2]
        		version      = vals[4][len(vals[4])-1] + '.' + vals[5]
                        
                        T = TruchasBasicRunTest(currdir,
        				        'static_drop.inp',
        					parallel_env,
        					1,
        					compiler,
        					compile_mode,
        					version,
        					truchas_bin)
                        print T.str(T.inpfile)
        
                        suite.addTest(T)
                        result = runner.run( suite )
        
    except:

	print "---> Test failed in some aspect <---"
	print "\nNature-Of-Error:",sys.exc_info()[0],"\n"
	if opt.d: raise
        
         
