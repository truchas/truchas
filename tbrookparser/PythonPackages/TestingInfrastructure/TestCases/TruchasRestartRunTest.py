#!/usr/bin/env python
"""
TruchasRestartRunTest

-----------------------------------------------------------------------------
   Purpose:
  
      TestCase to run a restarted Truchas simulation
  
   Public Interface:
  
      T = TruchasRestartRunTest(unittest.TestCase,TruchasBaseTest)
      T.setUp()
      T.shortDescription()
      T.tearDown()
      T.testBasicRun()
      T.str()
  
   Contains:  
      class TruchasRestartRunTest
            __init__(dir,inpfile,restartfile,parallel_env,
                     np,compiler,compile_mode,
                     truchas_bin,truchas_exe)
            setUp()
            shortDescription()
            tearDown()
            testRestartRun()
            str()
  
   Unit Test Block
  
   Author(s): Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""

import unittest
import os, sys, string, re, commands, string, platform, shutil

if __name__=='__main__':
    print "\n for component test in %s \n" %(__file__)
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    testingdir  = thisdir + '/../'
    sys.path.append(testingdir)
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

from ProgramEnvironment import createRunTimeProgram
from Runners            import TestRunner, RunTimeSpecs
from TruchasBaseTest    import TruchasBaseTest
from DataHandlers       import getDataCreator

class TruchasRestartRunTest(unittest.TestCase,TruchasBaseTest):

    def __init__(self,dir,
		 inpfile,
		 restartfile,
		 RunSpecs,
		 parallel_env  = 'serial',
		 np            = 1,
		 compiler      = 'lahey',
		 compile_mode  = 'debug',
                 truchas_bin   = '.',
		 truchas_exe   = None):

	unittest.TestCase.__init__(self,methodName='testRestartRun')
	self.inpfile      = inpfile
	self.restartfile  = restartfile
	self.dir          = dir
	self.parallel_env = parallel_env
	self.np           = np
	self.compile_mode = compile_mode
	self.compiler     = compiler
	self.truchas_bin  = truchas_bin
	self.truchas_exe  = truchas_exe
	#specifications from previous Truchas run
	self.prevrunspecs = RunSpecs
        self.version      = self.prevrunspecs.version
        self.success      = self.prevrunspecs.success
	#for formatting log file
	self.column1      = 0
	self.column2      = 5
	self.column3      = 10

	#construct executable name
	if self.truchas_exe == None:
	    self.truchas_exe = self.constructTruchasExe()
        #given executable extract parallel_env, compile_mode, compiler  
        else:
            self.deconstructTruchasExe()

	self.progname = '%s/%s' %(self.truchas_bin,self.truchas_exe)

	iend          = string.rfind(self.inpfile,".")
	self.basename = self.inpfile[0:iend] 

	self.logdir      = self.prevrunspecs.logdir
	self.logname     = 'RestartRun.log' 
        self.basenamelog = self.basename + '_logs' 
        self.logfile     = os.path.join(self.logdir,self.basenamelog,self.logname)

        if '_restart' in self.basename:
            r            = string.find(self.basename,'_restart')
            self.outdir  = self.basename[0:r]+'_restart_output'
        else:    
            self.outdir  = self.basename+'_restart_output'
	self.outfile     = self.basename+'.rrout'

        dsdir            = os.path.join(self.dir,self.logdir,self.outdir)
        self.datastore   = getDataCreator(dir=dsdir)

    def setUp(self):
	"defines default runtime environment, input file name" 
	
	os.chdir(self.dir)

        if self.logdir not in os.listdir(self.dir):
            os.mkdir(self.logdir)

        wdir = os.path.join(self.dir,self.logdir)
	os.chdir(wdir)            

        if self.basenamelog not in os.listdir(wdir):
            os.mkdir(self.basenamelog)

        logfl          = os.path.join(self.basenamelog,self.logname)        

	self.watcher   = open(logfl,'w')
	self.watcher.write('In RestartRun TestCase')
 	self.watcher.write('\n')
 	self.watcher.write('\n')

	tmp     = '*Setting up'
	str     = self.rjustln(self.column1, tmp)
	self.watcher.write(str)

	self.watcher.write(self.str(self.inpfile,self.restartfile,self.column2))
	
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

    def testRestartRun(self):
	"runs a Truchas restart simulation" 
	"(i.e runs truchas -o:output_dir -r:restartfile input_file)"
	
	tmp     = '*testMethod = testRestartRun'
	str     = self.rjustln(self.column1, tmp)
	self.watcher.write(str)
	self.watcher.write('\n')

	# runs and returns good output code
        inpfl    = '../%s' %(self.inpfile)
        if self.restartfile in os.listdir(os.getcwd()):
            restrtfl = '%s' %(self.restartfile)
        else:
            restrtfl = '../%s' %(self.restartfile)

        if self.parallel_env == 'parallel':
            args = ['-np %d '%(self.np), '-t', '-o:%s'%(self.outdir), '-r:%s'%(restrtfl), inpfl]
        else:
            args = ['-t', '-o:%s'%(self.outdir), '-r:%s'%(restrtfl), inpfl]
            
	T       = createRunTimeProgram(self.progname, args, self.outfile)

	retcode = T.run()
        output  = file(self.outfile,'r').read()
	self.watcher.write(T.str())

	if retcode!=0:
            outfl = os.path.join(os.getcwd(),self.outfile)
            tmp   = 'ERROR!! Program %s returned code of %d.' %(self.truchas_exe,retcode)
            tmp  += '\nSee %s for details!' %(outfl)
	    str   = self.rjustln(self.column2,tmp)
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

    def str(self,inpfile,restrtfile,column=15):
	"provides formatted info about RunTime specs for a given inpfile and restrtfile"

	info  = '\n'
	tmp   = 'Input file         : ' + inpfile
        info += self.rjustln(column, tmp)
        tmp   = 'Restart file       : ' + restrtfile
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
	tmp   = 'Previous BasicRunSpecs : ' + self.prevrunspecs.str()        
        info += self.rjustln(column, tmp)
	info += '\n'

	return info

if __name__=='__main__':
    suite       = unittest.TestSuite()
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
	    """
	    This is a Restart Test, therefore we need to access existing outputs.t-* directories.
	    To fully test, each executable should be used to restart each existing output.  This
	    results in n^2 analyses.  To accomplish this, first loop over the output directories.
	    Then loop over the executables in the truchas/bin directory.  

	    The above approach results in a test suite that will take more than the desired amount
	    of time for a component test.  Therefore, this example execuates useing the first output 
	    directory and the first executable only.  
	    
	    Sections of the following code have been marked for use with a full test of restart 
	    files and with the smaller component test run.
	     
	    Inspect test output to determine which output directory and truchas exe have been used.
	    """
		
            # create a list of existing output directories
            output_dirs = []
            for dir in os.listdir(currdir):
                if 'outputs.t-' in dir and os.path.isdir(dir): output_dirs.append(dir)
        
	    # Full test: loop over available output directories
            #for dir in output_dirs:
	    # Component test: use the first directory in the list
            for dir in output_dirs[0:1]:

                logdir  = dir
                logfile = logdir+'/static_drop_logs/BasicRun.log'    

		# Full test: loop over the exe files
                #for exefile in os.listdir(truchas_bin):
		# Component test: use just the first exe file
		for exefile in os.listdir(truchas_bin)[0:1]:
                    
                    suite   = unittest.TestSuite()
                        
                    vals    = string.splitfields(exefile,'.')
                    prevals = string.splitfields(logdir,'.')
                    if prevals[1][2:] == string.lower(platform.system()):
                        parallel_env = prevals[4]      #'serial'
                        np           = 1
                        compiler     = prevals[3]      #'lahey'
                        compile_mode = prevals[5][:-2] #'opt'
        		version      = prevals[5][len(prevals[5])-1] + '.' + prevals[6]
			# get version from RunTimeSpecs
                        LastRunSpecs = RunTimeSpecs(logdir,
        					    logfile,
        					    parallel_env,
        					    1,
        					    compiler,
        					    compile_mode,
        					    version)
        
        
                    parallel_env = vals[3]
                    compiler     = vals[2]
                    compile_mode = vals[4][:-2]
                                
                    T            = TruchasRestartRunTest(currdir,
        				                 'static_drop.inp',
        						 'restart.bin',
        						  LastRunSpecs,
        						  parallel_env,
        						  1,
        						  compiler,
        						  compile_mode,
        						  truchas_bin)
                    print T.str(T.inpfile,T.restartfile)
        
                    suite.addTest(T)
                    result = runner.run( suite )

    except:

	print "---> Test failed in some aspect <---"
	print "\nNature-Of-Error:",sys.exc_info()[0],"\n"
	if opt.d: raise


