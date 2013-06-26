#!/usr/bin/env python
"""
RunBaseTestSuite

-----------------------------------------------------------------------------
   Purpose:
  
      Provides basic functions to run the following TestSuites:
      -Capability
      -PostProcessor
      -Hierarchical
      -List
  
   Public Interface:
  
      T = RunBaseTestSuite()
      T.defineUsage(self)
      T.getOptions(self)
      T.parseOptions(self)
      T.setLogDir(self,outputdir)
      T.checkOptions(self)
      T.loadAndRun(self,loader,currdir,executable,fpsummary)
      T.cleanUp(self,currdir)
      T.dirname(self,currdir,logdir)
      T.deconstructTruchasExe(self,executable)
      T.printRunTimeSpecs(self,executable,logdir)
      T.printSummary(self,thedir,success)
  
   Contains:
   
      class RunBaseTestSuite
            defineUsage(self)
            getOptions(self)
            parseOptions(self)
            setLogDir(self,outputdir)
            checkOptions(self)
            loadAndRun(self,loader,currdir,executable,fpsummary)
            cleanUp(self,currdir)
            dirname(self,currdir,logdir)
            __usageExit(self,msg)
            __checkChoices(self)
            deconstructTruchasExe(self,executable)
            __rjustln(self,column,s)
            __printChoices(self)
            __getValidExecutables(bindir)
            __printExecutables(self)
            printRunTimeSpecs(self,executable,logdir)
            
   Author: Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""

import os, sys, string, unittest, sets, commands, fnmatch, shutil

if __name__=='__main__':
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    testingdir  = thisdir + '/../'
    sys.path.append(testingdir)
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

from optparse    import OptionParser 
from Loggers     import TestLoggers
from PYTHONutils import getAbortString

import mytestrunner
TestRunner = mytestrunner.TextTestRunner
#if self.debug:
#TestRunner = unittest.TextTestRunner
        
class getArch:

    def __init__(self):

        if 'linux' in sys.platform:
            self.platform  = 'linux'
            self.arch      = commands.getstatusoutput('uname -i')[1]
            if 'unknown' in self.arch:
                self.arch  = commands.getstatusoutput('uname -m')[1]
            # hack for ubuntu - doesn't know uname -i or -p, -m returns i686
            if self.arch == 'i686':
                self.arch = 'i386'

        elif 'darwin' in sys.platform:
            self.platform  ='darwin'
            self.arch      = commands.getstatusoutput('uname -p')[1]

        elif 'osf1' in sys.platform:
            self.platform  ='osf1'
            tmp            = commands.getstatusoutput('uname -a')[1]
            self.arch      = string.split(tmp,' ')[4]

        else:
            self.arch      = None
            self.platform  = None


class RunBaseTestSuite:
    """
    Base class providing base methods for running the following TestSuites
      -Capability
      -PostProcessor
      -Hierarchical
      -List    
    """

    def defineUsage(self):

        self.baseusage = """\
        python RunTestSuite.py [options] 

        Examples:
        python RunTestSuite.py
             - runs TestSuite using all existing executables in 
             the binary directory truchas/bin and all TestCases
             that exist in subdirectories underneath the current working dir.
             Here the location of the truchas directory must be specified in the RunTestSuite.py file.

        python RunTestSuite.py -d dir
             - runs TestSuite using all existing executables in 
             the binary directory truchas/bin. The resulting TestSuite contains all TestCases
             that exist in subdirectories underneath the directory 'dir'.
                                             
        python RunTestSuite.py -b bindir
             - runs TestSuite using all existing executables in 
             the directory 'bindir'
                                             
        python RunTestSuite.py -b bindir -e t-linux.intel.lahey.serial.dbg-2.1.dev
             - runs TestSuite using the t-linux.intel.lahey.serial.dbg-2.1.dev
             executable in the directory 'bindir'
                                             
        python RunTestSuite.py -b bindir -m serial-dbg
             - runs TestSuite using the serial-debug executable in the
             directory 'bindir'
                                             
        python RunTestSuite.py -b bindir -c nag
             - runs TestSuite using executable(s) created from the nag 
             compiler that exists in the directory 'bindir'
                                             
        python RunTestSuite.py --nr
             - don't run any BasicRun, RestartRun or PostProcessorRun TestCases,
             only run the Capability TestCases. This option provides a
             valid result only when an outputs.* directory exists in directory.
        """

    def getOptions(self):

        thisdir       = os.path.abspath(os.path.dirname(__file__))
        testingdir    = thisdir + '/../'
        #truchas_bin   = testingdir + '/../../../bin'
        #truchas_bin   = os.path.abspath(truchas_bin)
        #self.srcdir   = os.path.abspath(thisdir + '/../../../../regressiontests')

        #parse command line arguments
        self.mycommands    = OptionParser(usage=self.usage)
        self.mycommands.add_option("-d","--dir", action="store",type="string",
                              help="DIR = directory containing all TestCases to be run")
        self.mycommands.add_option("-c","--compiler", action="store",type="string",
                              help="Choice of compiler. For linux COMPILER = nag or COMPILER = lahey. For osf1 COMPILER = compaq.")
        self.mycommands.add_option("-m","--mode", action="store",type="string",
                              help="Compilation mode. MODE = serial-dbg, serial-opt, parallel-dbg, parallel-opt.")
        self.mycommands.add_option("-b","--binarydir", action="store",type="string",
                              help="BINARYDIR = directory containing the Truchas executable(s)")
        self.mycommands.add_option("-e","--executable", action="store",type="string",
                              help="EXECUTABLE = full name of Truchas executable that lives in BINARYDIR specified above")
        self.mycommands.add_option("-n","--np","--processors", action="store",type="int",
                              help="NP = number of processors used in a parallel run")
        self.mycommands.add_option("--nr","--norerun", action="store_true",
                              help="BasicRun, RestartRun and PostProcessorRun TestCases will not be run,only the Capability TestCase(s) will only be run.")

    def parseOptions(self):
        
        (self.options, args) = self.mycommands.parse_args()


    def setLogDir(self,outputdir):

        if not os.path.isdir(outputdir):
            os.mkdir(outputdir)

        runnername = 'TestSuiteRunner'
        if self.debug:
            self.debugdir = 'for_debugging'
            if not os.path.isdir(os.path.join(outputdir,self.debugdir)):
                os.mkdir(os.path.join(outputdir,self.debugdir))
            
                debugdir      = os.path.abspath(os.path.join(outputdir,self.debugdir))
                runnername    = '%s/TestSuiteRunner' %(debugdir)
        
        self.logger   = TestLoggers(runnername,debugging=self.debug).getThem()

    def checkOptions(self):
       
        print 'self.bindir=', self.bindir 
        self.__getValidExecutables(self.bindir)

        if self.options.binarydir  != None:
            self.bindir = self.options.binarydir
            self.__getValidExecutables(self.bindir)
        if self.options.mode       != None:
            T = string.split(self.options.mode,'-')
            errstring     = 'Error!! Compilation mode must be of the form *-*' 
            assert len(T) == 2, errstring
            self.parallel_env = T[0]
            self.compile_mode = T[1]
        if self.options.executable != None:
            self.exevalue     = string.split(self.options.executable,'/')[-1]
            self.executables  = []
            self.executables.append(self.exevalue)            
        if self.options.compiler != None:
            self.compiler     = self.options.compiler    
        if self.options.np  != None:
            self.np           = self.options.np
        if self.options.nr  != None:
            self.no_rerun     = 1
        if self.options.dir != None:
            self.currdir      = os.path.abspath(self.options.dir)

        self.logger.debug(self.__printChoices())

        self.__checkChoices()

        self.logger.debug(self.__printExecutables())

    def loadAndRun(self,loader,executable,fpsummary=sys.stdout):

        if not self.no_rerun:
            #load and run TestSuite given executable resulting from user input choices 
            self.deconstructTruchasExe(executable)

        X = loader(self.no_rerun,
                   self.srcdir,
                   self.currdir,
                   self.debug,
                   self.parallel_env,
                   self.np,
                   self.compiler,
                   self.compile_mode,
		   self.version,
                   self.bindir,
                   executable)
                
        X.getSuite()

        if not os.path.isdir(os.path.abspath(os.path.join(self.currdir,X.logdir))):
            os.mkdir(os.path.abspath(os.path.join(self.currdir,X.logdir)))

        thedir  = self.dirname(self.currdir,X.logdir)
        info    = self.printRunTimeSpecs(executable,thedir)

        self.logger.debug(info)
        self.logger.debug(X.suite)

        fpverb = sys.stdout
        fpverb.write(info)

        runner  = TestRunner(stream=fpverb,verbosity=2)
        result  = runner.run( X.suite )

        status_file = os.path.abspath(os.path.join(self.currdir,X.logdir,'STATUS'))
        fpstatus    = open(status_file,'w')

        if result.wasSuccessful():
            fpstatus.write('PASS')
        else:
            fpstatus.write('FAIL')

        try:
            t         = '%5.3f' %(result.timeTaken)
            timeTaken = '(' + t + 's)'
        except:
            timeTaken = 0
            
        info = self.printSummary(thedir,result.wasSuccessful(),timeTaken)
        fpsummary.write(info)
                    
        fpstatus.close()


    def cleanUp(self):

        filepatterns  = ['outputs.t-*', 'for_debugging']

        files = os.listdir(self.currdir)

        for file in files:
            for pattern in filepatterns:
                if fnmatch.fnmatch(file,pattern):

                    if file in files and os.path.isdir(file):

                        try:
                            shutil.rmtree(os.path.join(self.currdir,file))
                        except:
                            cmd = 'rm -r %s/* > /dev/null  ' %(os.path.join(self.currdir,file))
                            os.system(cmd)
                            cmd = 'rm -rf  %s/ > /dev/null  ' %(os.path.join(self.currdir,file))
                            os.system(cmd)
                            self.logger.debug(self.currdir)
                            self.logger.debug(cmd)
                    if file in files and os.path.isfile(file):
                        try:
                            os.remove(os.path.join(self.currdir,file))
                        except:
                            cmd = 'rm %s > /dev/null' %(os.path.join(self.currdir,file))
                            os.system(cmd)                            
                            
                    files = os.listdir(self.currdir)


    def dirname(self,currdir,logdir):
        "for outputting runtime specifications"
        
        return logdir


    def __usageExit(self, msg=None):
        if msg: print msg
        print self.usage % self.__dict__
        sys.exit(2)

    def __checkChoices(self):
        "checks for valid run choices creates a list of valid executables"
        
        abortstring   = getAbortString()

# this is silly.  __getValidExecutables should have checked itself!
#        #check to ensure there are executables in bindir
#        errstring     = '\n \n    No executables exist in %s     \n' %(self.bindir)
#        errstring    += abortstring 
#        assert len(self.executables) > 0, errstring

        #if executable specified...check validity
        if self.exevalue != None:
            errstring  = '\n \n    The executable %s is not in %s  \n' %(self.exevalue,self.bindir)
            errstring += abortstring
            valid_exe  = 0
            exelist    = self.executables
	    self.executables = []
            if self.exevalue in exelist:
                valid_exe = 1
		self.executables.append(self.exevalue)
            assert valid_exe == 1, errstring

            errstring  = '\n \n Specify the full name of the executable i.e "-e t-linux.i386.nag.serial.opt-2.2.dev" instead of "-e truchas"  \n' 
            errstring += abortstring

            assert 'bin/truchas' not in self.exevalue, errstring

	    return

	compilerexes   = self.executables
	compilemodexes = self.executables
	penvexes       = self.executables
        npexes         = self.executables

        #if compiler specified....check validity
        if self.compiler != None:
            errstring  = '\n \n  An executable compiled with %s is not in %s  \n' %(self.compiler,self.bindir)
            errstring += abortstring
            valid_compiler = 0
	    compilerexes   = []   
            for file in self.executables :
                if self.compiler in file :
                    valid_compiler = 1
		    compilerexes.append(file)
            assert valid_compiler == 1, errstring

	
        #if mode specified....check compile mode validity
        if self.compile_mode != None:
            errstring  = '\n \n  An executable compiled in %s-%s mode is not in %s  \n' %(self.parallel_env,self.compile_mode,self.bindir)
            errstring += abortstring
            valid_compilemode = 0
	    compilemodexes    = []
            for file in self.executables:
                if self.compile_mode in file :
                    valid_compilemode = 1
		    compilemodexes.append(file)
            assert valid_compilemode == 1, errstring
            
        #if mode specified....check parallel environment validity
        if self.parallel_env != None:
            errstring  = '\n \n  An executable compiled in %s is not in %s  \n' %(self.parallel_env,self.bindir)
            errstring += abortstring
	    penvexes       = []
            valid_parallel = 0
            for file in self.executables:
                if self.parallel_env in file :
                    valid_parallel = 1
		    penvexes.append(file)
            assert valid_parallel == 1, errstring

        #if np specified....check for parallel executable

        parallel = 0
        for f in self.executables:
            if 'parallel' in f:
                parallel = 1

        if self.np > 1:         
            npexes = []
        for f in self.executables:
            if 'parallel' in f and self.np > 1:
                npexes.append(f)
            
        if self.np > 1 and not self.no_rerun:
            if not parallel:
                errstring  = '\n \n   WARNING!! Chosen np=%i but no parallel executable in %s  \n' %(self.np,self.bindir)
                errstring += '                   Ignoring np value and running serial executable(s)......\n'
                print errstring
            if self.parallel_env == 'serial':
                errstring  = '\n \n   WARNING!! Chosen np=%i but have chosen %s mode  \n' %(self.np,self.parallel_env)
                errstring += '                 Ignoring np value and running serial executable......\n'
                print errstring


	#find the intersection between the parallel env, compiler, compile mode options
	s1               = sets.Set(self.executables).intersection(sets.Set(compilerexes))
	s2               = sets.Set(s1).intersection(sets.Set(compilemodexes))
	s3               = sets.Set(s2).intersection(sets.Set(penvexes))
        self.executables = s3

	errstring  = '\n \n      Your combination of choices requires an executable  \n'
	errstring += '      that does not exist in %s.                          \n' %(self.bindir)
        errstring += '      Please recheck your choices!                        \n' 
	errstring += abortstring

	assert len(self.executables) > 0, errstring

	return 

    def deconstructTruchasExe(self,executable):
        "given executable..gets parallel_env, compile_mode, compiler"

        tmp               = string.split(executable,'.')
        self.compiler     = tmp[2]
        self.parallel_env = tmp[3]
        tmp2              = string.split(tmp[4],'-')
        self.compile_mode = tmp2[0]
	#BROKENself.version      = tmp2[1] + '.' + tmp[5]
	self.version      = '2' + '.' + '8'

    def __rjustln(self,column,s):

        tmp = s.rjust(column + len(s))
        tmp += '\n'

        return tmp
    
    def __printChoices(self):

	info  = '\n'
        tmp   = '*Preparing to run Truchas Capability TestSuite'
        info += self.__rjustln(self.column, tmp)
        tmp   = 'Developer choices.......'
        info += self.__rjustln(self.column, tmp)
        tmp   = 'Current dir: %s ' %(self.currdir)
        info += self.__rjustln(self.column, tmp)
	tmp   = 'Binary dir : %s ' %(self.bindir)
        info += self.__rjustln(self.column, tmp)
        tmp   = 'Available executables:'
        info += self.__rjustln(self.column, tmp)
        for exe in self.executables:
            tmp   = exe
            info += self.__rjustln(self.column, tmp)
	tmp   = 'Compiler       : ' + str(self.compiler)
        info += self.__rjustln(self.column, tmp)
	tmp   = 'Executable type: ' + str(self.parallel_env) + '-' + str(self.compile_mode)
        info += self.__rjustln(self.column, tmp)
        tmp   = 'np             : ' + str(self.np)
        info += self.__rjustln(self.column, tmp)

        return info

    def __getValidExecutables(self,bindir):

        self.executables = []
        pltfrm           = sys.platform
        A                = getArch()
        arch             = A.arch
        pltfrm           = A.platform
        
        for file in os.listdir(bindir):
            if 't-' in file and pltfrm in file and arch in file: # and compiler in file:
                self.executables.append(file)

        if len(self.executables) == 0:
            print
            print '    can\'t find any suitable executables in %s' % bindir
            print '    looking for filenames containing "t-" and "%s" and "%s"' % (pltfrm, arch)
            print
            print '    ABORT'
            print
            raise AssertionError

    def __printExecutables(self):

        info  = '\n \n'
        tmp   = 'Executable(s) to be run : '
        info += self.__rjustln(self.column1, tmp)
        for exe in self.executables:
            tmp   = exe
            info += self.__rjustln(self.column2, tmp)
            
        return info

   
    def printRunTimeSpecs(self,
                          executable,
                          logdir): 

	separator2 = '-' * 80
        info    = '\n'
	info   += separator2
	info   += '\n'
	tmp    = '{TruchasRunTime}  '
        length = len(tmp)
	tmp   += 'Binary dir     : %s ' %(self.bindir)
	info   += tmp
	info   += '\n'
	tmp    = 'Compiler       : ' + str(self.compiler)
        info  += self.__rjustln(length, tmp)
	tmp    = 'Executable type: ' + str(self.parallel_env) + '-' + str(self.compile_mode)
        info  += self.__rjustln(length, tmp)
        tmp    = 'Executable name: ' + str(executable)
        info  += self.__rjustln(length, tmp)
	tmp    = 'Outputs dir    : ' + str(logdir)
        info  += self.__rjustln(length, tmp)
        if self.parallel_env == 'parallel':
	    tmp   = 'np             : ' + str(self.np)
            info += self.__rjustln(length, tmp)
            
        info += '\n'
	return info


    def printSummary(self,thedir,success,timetaken=0):

        #strip off outputs.t-* directory
        L      = string.split(thedir,'/')
        M      = L[0:-1]
        newdir = string.join(M,'/')

        if success:
            result = '......PASS' 
        else:
            result = '......FAIL' 

        result += '  %s' %(str(timetaken))
        
        tmp     = result.rjust(70-len(str(newdir)))
        #info   = str(newdir) + tmp
        info    = tmp
        info   += '\n'
            
        return info
    
