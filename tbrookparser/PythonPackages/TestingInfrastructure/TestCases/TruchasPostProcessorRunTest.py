#!/usr/bin/env python
"""
TruchasPostProcessorRunTest

-----------------------------------------------------------------------------
   Purpose:
  
      TestCase to run a Truchas postprocessor macro file
  
   Public Interface:
  
      T = TruchasPostProcessorRunTest(unittest.TestCase,TruchasBaseTest)
      T.setUp()
      T.shortDescription()
      T.tearDown()
      T.testPostProcessorRun()
      T.str()
  
   Contains:  
      class TruchasPostProcessorRunTest
            __init__(dir,toolsdir,macfile,RunSpecs)
            setUp()
            shortDescription()
            tearDown()
            testPostProcessorRun()
            str()
  
   Unit Test Block
  
   Author(s): Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""

import unittest
import os, sys, string, re, platform, shutil

if __name__ == '__main__':
    print "\n for component test in %s \n" %(__file__)
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    testingdir  = os.path.abspath(thisdir + '/../')
    sys.path.append(testingdir)
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

from ProgramEnvironment import createRunTimeProgram
from Runners            import TestRunner, RunTimeSpecs
from TruchasBaseTest    import TruchasBaseTest
from processInput       import processInput

class TruchasPostProcessorRunTest(unittest.TestCase,TruchasBaseTest):

    def __init__(self,dir,
		 toolsdir,
		 macfile,
		 RunSpecs):

	unittest.TestCase.__init__(self,methodName='testPostProcessorRun')
	self.macfile      = macfile
	self.dir          = dir
	#directory where tools lives 
	self.toolsdir     = toolsdir 
	#specifications from a previous Truchas run
	self.runspecs     = RunSpecs
	self.driver       = 'python'
	self.progname     = '%s %s/scripts/TBrookParse.py' %(self.driver,self.toolsdir)
	#for formatting log file
	self.column1      = 0
	self.column2      = 5	
        self.fail         = 'Error'
        
	iend              = string.rfind(self.macfile,".")
	self.basename     = self.macfile[0:iend]

	self.logdir       = self.runspecs.logdir
        self.logname      = 'PostProcessor.log'
        self.basenamelog  = string.split(self.runspecs.logfile,'/')[1]
	self.logfile      = os.path.join(self.logdir,self.basenamelog,self.logname)
	self.outfile      = self.basename+'.ppout'

    def setUp(self):
	"defines default postprocessor runtime environment" 
	
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
	self.watcher.write('In PostProcessorRun TestCase')
 	self.watcher.write('\n')
 	self.watcher.write('\n')

	tmp     = '*Setting up'
	str     = self.rjustln(self.column1, tmp)
	self.watcher.write(str)

	self.watcher.write(self.str(self.macfile,self.column2))


    def shortDescription(self):
	"modified unittest short description to allow class name and input file name"

        classname = self.__class__.__name__
        classname = re.sub('Test[ers]*$', '', classname)
	inputname = self.macfile
        docname   = '['+classname +']' + ' ' + inputname
        return  docname.ljust(40) + ': ' + str(unittest.TestCase.shortDescription(self))

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

	shutil.move(self.outfile,self.logdir)

        # look for the resulting viz,diagnostic file and copy it to logdir
        for ext in ['.bin','.ascii','.txt']: 
	    filename = self.basename + ext
            if os.path.isfile(os.path.join(self.dir,filename)):
                tmp     = 'Copying %s to %s/ ' %(filename, self.logdir)
	        str     = self.rjustln(self.column2, tmp)
	        self.watcher.write(str)
                shutil.copy(filename,self.logdir)
        """

    def testPostProcessorRun(self):
	"runs a Truchas postprocessor simulation" 
	"(i.e runs processInput with macrofile as input)"

	tmp      = '*testMethod = testPostProcessorRun'
	str      = self.rjustln(self.column1, tmp)
	self.watcher.write(str)
	self.watcher.write('\n')

        macfl    = '../%s' %(self.macfile)

        inbuffer = '@' + macfl + ' q'
        fp       = open(self.outfile,'w')
        try:
            T    = processInput(inbuffer,fp)
        except:
            outfl = os.path.join(os.getcwd(),self.outfile)
            tmp   = 'ERROR!! Program processInput with inbuffer %s failed. See %s.' %(inbuffer,outfl)
	    str   = self.rjustln(self.column2,tmp)
            self.watcher.write('\n')
	    self.watcher.write(str)
	    self.watcher.write('\n')
            errstr = '\n \n' + tmp + '\n'
	    raise ValueError(errstr)
        fp.close()

        output = file(self.outfile,'r').read()
        
        if (self.fail in output) :
            outfl = os.path.join(os.getcwd(),self.outfile)
	    tmp   = 'ERROR!! Postprocessor terminated in error. See %s for details!' %(outfl)
            str   = self.rjustln(self.column2, tmp)
            self.watcher.write('\n')
	    self.watcher.write(str)
            self.watcher.write('\n')

            errstr = '\n \n' + tmp + '\n'
	    raise ValueError(errstr)
        
        
    def str(self,macfile,column=15):
	"provides formatted info about postprocessor runtime specs for a given macro file"

	info  = '\n'
	tmp   = 'Macro file(s)     : ' + macfile
        info += self.rjustln(column, tmp)
	tmp   = 'Output file       : ' + self.outfile
        info += self.rjustln(column, tmp)
	tmp   = 'Current dir       : ' + self.dir
        info += self.rjustln(column, tmp)
	tmp   = 'tools dir         : ' + self.toolsdir
        info += self.rjustln(column, tmp)
	tmp   = 'PostProcessor exe : ' + self.progname
        info += self.rjustln(column, tmp)
	tmp   = 'Log file          : ' + self.logfile
        info += self.rjustln(column, tmp)
        info += '\n'
	tmp   = 'Previous RunSpecs : ' + self.runspecs.str()
        info += self.rjustln(column, tmp)
	info += '\n'

	return info

if __name__=='__main__':
    currdir     = os.getcwd()
    toolsdir    = thisdir + '/../../../'
    toolsdir    = os.path.abspath(toolsdir)

    TestRunner  = unittest.TextTestRunner

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
    #determine ru        ntime directories
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
            
                    # run PostProcessor tests
                    suite      = unittest.TestSuite()
                    macrofiles = ['ensight.mac','gmv.mac','gmv_ascii.mac','query.mac','restart.mac','stat.mac','tecplot.mac','vtk.mac']
                    for macro_file in macrofiles:
                        T  = TruchasPostProcessorRunTest(currdir,toolsdir,macro_file,LastRunSpecs) 
                        suite.addTest(T)
                    print T.str(str(macrofiles))
        
                    runner = TestRunner(verbosity=2)
                    result = runner.run( suite )
            
    except:

	print "---> Test failed in some aspect <---"
	print "\nNature-Of-Error:",sys.exc_info()[0],"\n"
	if opt.d: raise
        
        
         
