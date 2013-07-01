#!/usr/bin/env python
"""
 RunHierarchicalTestSuite

-----------------------------------------------------------------------------
   Purpose:
  
      Creates and runs a TestSuite containing:
      multiple Capability TestSuites and multiple PostProcessor
      TestSuites by implementing a TruchasCapabilityLoader
      and a TruchasPostProcessorLoader recursively from each subdirectory
  
   Public Interface:

      T = RunHierarchicalTestSuite(RunBaseTestSuite)
      T.options(self)
      T.setRunTimeEnvironment(self,no_rerun,np,compiler,
                              parallel_env,compile_mode,
                              bindir,executables)
      T.hierarchicalCleanUp(self)
      T.hierarchicalLoadAndRun(self,executable)
      T.dirname(self,currdir,logdir)

  
   Contains:
   
      class RunHierarchicalTestSuite(RunBaseTestSuite):
            __init__(self,argv,currdir,
                     caploader,pproloader,
                     fpsummary = sys.stdout)
            options(self)
            setRunTimeEnvironment(self,no_rerun,np,compiler,
                                  parallel_env,compile_mode,
                                  bindir,executables)
            hierarchicalCleanUp(self)
            hierarchicalLoadAndRun(self,executable)
            __searchSubDirsCleanUp(self,currdir)
            __searchSubDirsLoadAndRun(self,currdir,executable)
            dirname(self,currdir,logdir)

   Unit Test Block
  
   Author: Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""
import os, sys, string, fnmatch

thisdir     = os.path.abspath(os.path.dirname(__file__))
testingdir  = thisdir + '/../'

if __name__=='__main__':
    sys.path.append(testingdir)
    sys.path.append(testingdir+'Loaders')
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

from Loaders import TruchasCapabilityLoader,TruchasPostProcessorLoader

from RunBaseTestSuite import RunBaseTestSuite
import mytestrunner
TestRunner = mytestrunner.TextTestRunner

class RunHierarchicalTestSuite(RunBaseTestSuite):
    """
    Creates and runs a TestSuite by searching recursively
    through subdirectories of 'currdir' and applying the
    TruchasCapabilityLoader and TruchasPostProcesorLoader
    depending on the contents of the given directory.
    """

    def __init__(self,
                 argv,
                 srcdir,
                 currdir,
                 caploader  = TruchasCapabilityLoader,
                 pproloader = TruchasPostProcessorLoader,
                 debug      = 0): 
                 #currdir,
                 #caploader,
                 #pproloader,
                 #fpsummary = sys.stdout,
                 #debug     = 0): 

        self.currdir      = currdir
        self.srcdir       = srcdir
        truchas_bin       = testingdir + '/../../../bin'
        truchas_bin       = os.path.abspath(truchas_bin)
        self.parallel_env = None
        self.np           = 1
        self.compile_mode = None
        self.compiler     = None
	self.version      = '2.3'
        self.executables  = []
        self.exevalue     = None
	self.bindir       = truchas_bin
        self.no_rerun     = 0
        #for formatting output
        self.column       = 5
        self.column1      = 10
        self.column2      = 15

        self.capfilePat   = 'CapabilityTest.py'
        self.pprofilePat  = '*.mac'
	# To Developer:  Make sure each specialty loader has a unique name 
	#                within a given suite
        self.myloaderPat  = 'My*Loader.py'
        self.caploader    = caploader
        self.pproloader   = pproloader
        self.debug        = debug

        self.setLogDir(outputdir = self.currdir)

        self.timeTaken    = 0.0

    def options(self):

        self.defineUsage()
        self.usage        = self.baseusage.replace('RunTestSuite.py','RunHierarchicalTestSuite.py')        
        self.getOptions()
        self.parseOptions()
        self.checkOptions()

    def setRunTimeEnvironment(self,no_rerun,np,compiler,parallel_env,compile_mode,version,bindir,executables):
        "useful when command line options differ to standard options - this provides same info as attained from options"

        self.no_rerun     = no_rerun
        self.np           = np
        self.compiler     = compiler
        self.parallel_env = parallel_env
        self.compile_mode = compile_mode
	self.version      = version
        self.bindir       = bindir
        self.executables  = executables

    def cleanUp(self): 
        
        self.__searchSubDirsCleanUp(self.currdir)


    def loadAndRun(self,
                   loader,
                   executable,
                   fpsummary  = sys.stdout,
                   verbtofile = 0):

        self.__searchSubDirsLoadAndRun(loader,self.srcdir,self.currdir,executable,fpsummary,verbtofile) 

    def dirname(self,currdir,logdir):
        "for outputting runtime specifications"

        J      = string.split(logdir, '/')
        L      = string.split(currdir + '/' + J[-1:][0],'/')
        M      = string.split(self.currdir,'/')
        topdir = M[-1:][0]
        for x in L:
            if topdir in x:
                index = L.index(x)
        thedir = string.join(L[index+0*1:],'/')

        return thedir
    

    def __searchSubDirsCleanUp(self,currdir):

        if os.path.isdir(currdir):
            self.__cleanUp(os.path.abspath(currdir))
            thedirs       = os.listdir(currdir)
            for f in thedirs:
                self.__searchSubDirsCleanUp(os.path.abspath(os.path.join(currdir,f)))   
        else:
            return

    def __searchSubDirsLoadAndRun(self,loader,srcdir,currdir,executable,fpsummary,verbtofile): 

        if os.path.isdir(currdir):
            thedirs          = os.listdir(currdir)
            storecapmatch    = 0
            storeloadermatch = 0
            storeppromatch   = 0

            for f in thedirs:
                loadermatch = fnmatch.fnmatch(f, self.myloaderPat)
                capmatch    = fnmatch.fnmatch(f, self.capfilePat)
                ppromatch   = fnmatch.fnmatch(f, self.pprofilePat)
                if loadermatch:
                    storeloadermatch = 1                    
                if capmatch:
                    storecapmatch    = 1
                if ppromatch:
                    storeppromatch   = 1                

            for f in thedirs:
                capmatch    = fnmatch.fnmatch(f, self.capfilePat)
                loadermatch = fnmatch.fnmatch(f, self.myloaderPat)
                ppromatch   = fnmatch.fnmatch(f, self.pprofilePat)

                if loadermatch:

                    loadername = f
                    #strip the .py extension...
                    if ".py" in f:
                        loadername = f[0:-3] 
                    sys.path.append(currdir)
                    __import__(loadername)
                    module = sys.modules[loadername]
                    loader = vars(module)[loadername]

                    self.__loadAndRun(loader,
                                      srcdir,
                                      currdir,
                                      executable,
                                      fpsummary,
                                      verbtofile)
                    
                if capmatch and not storeloadermatch:
                    self.__loadAndRun(self.caploader,
                                      srcdir,
                                      currdir,
                                      executable,
                                      fpsummary,
                                      verbtofile)
                    
                if ppromatch and not storecapmatch and not storeloadermatch:
                    self.__loadAndRun(self.pproloader,
                                      srcdir,
                                      currdir,
                                      executable,
                                      fpsummary,
                                      verbtofile)

            for f in os.listdir(currdir):
                self.__searchSubDirsLoadAndRun(loader,
                                               os.path.abspath(os.path.join(srcdir,f)),
                                               os.path.abspath(os.path.join(currdir,f)),
                                               executable,
                                               fpsummary,
                                               verbtofile)   
        else:
            return

    def __loadAndRun(self,
                     loader,
                     srcdir,
                     currdir,
                     executable,
                     fpsummary  = sys.stdout,
                     verbtofile = 0):

        if not self.no_rerun:
            #load and run TestSuite given executable resulting from user input choices 
            self.deconstructTruchasExe(executable)

        X = loader(self.no_rerun,
                   srcdir,
                   currdir,
                   self.debug,
                   self.parallel_env,
                   self.np,
                   self.compiler,
                   self.compile_mode,
		   self.version,
                   self.bindir,
                   executable)


        X.getSuite()

        #write out directory name before testcase is run...
        thedir   = self.dirname(currdir,X.logdir)
        L        = string.split(thedir,'/')
        M        = L[0:-1]
        shortdir = string.join(M,'/')
        fpsummary.write(shortdir)

        if not os.path.isdir(os.path.abspath(os.path.join(currdir,X.logdir))):
            os.mkdir(os.path.abspath(os.path.join(currdir,X.logdir)))

        info    = self.printRunTimeSpecs(executable,thedir)

        self.logger.debug(info)
        self.logger.debug(X.suite)

        if verbtofile:
            verb_file = os.path.abspath(os.path.join(currdir,X.logdir,'SUMMARY'))
            fpverb    = open(verb_file,'w')
        else:
            fpverb     = sys.stdout

        fpverb.write(info)

        runner  = TestRunner(stream=fpverb,verbosity=2)
        result  = runner.run( X.suite )

        status_file = os.path.abspath(os.path.join(currdir,X.logdir,'STATUS'))
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

        self.timeTaken += result.timeTaken
        
        info = self.printSummary(thedir,result.wasSuccessful(),timeTaken)
        fpsummary.write(info)
                    
        fpstatus.close()
        if verbtofile:
            fpverb.close()

    def __cleanUp(self,currdir):

        filepatterns = ['for_debugging','outputs.t-*'] #,'for_debugging']
        files        = os.listdir(currdir)

        for file in files:
            for pattern in filepatterns:
                if fnmatch.fnmatch(file,pattern):
                    if file in files and os.path.isdir(os.path.join(currdir,file)):
                        try:
                            shutil.rmtree(os.path.join(currdir,file))
                        except:
                            cmd = 'rm -r %s/* > /dev/null  ' %(os.path.join(currdir,file))
                            os.system(cmd)
                            cmd = 'rm -rf  %s/ > /dev/null  ' %(os.path.join(currdir,file))
                            os.system(cmd)
                            self.logger.debug(currdir)
                            self.logger.debug(cmd)
                    if file in files and os.path.isfile(file):
                        try:
                            os.remove(os.path.join(currdir,file))
                        except:
                            cmd = 'rm %s > /dev/null' %(os.path.join(currdir,file))
                            os.system(cmd)                            
                            
                    files = os.listdir(currdir)


    
if __name__=="__main__":

    argv = sys.argv[1:]
    outfile   = '%s_component_test' %(__file__)
    fpsummary = open(outfile,'w')
    curdir    = os.path.abspath(os.getcwd()+'/../examples/hierarchical_example')
    X    = RunHierarchicalTestSuite(argv,
                                    srcdir  = curdir,
                                    currdir = curdir,
                                    debug   = 0)                                    
    X.options()    
    for executable in X.executables:
        print
        print 'EXECUTABLE : %s' %(str(executable))
        print
        X.loadAndRun(loader     = None, 
                     executable = executable,
                     fpsummary  = fpsummary)

    fpsummary.close()
    
    X.cleanUp()


