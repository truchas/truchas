#!/usr/bin/env python
"""
RunPostProcessorTestSuite

-----------------------------------------------------------------------------
   Purpose:
  
      Wrapper to allow developer to create and run a PostProcessor TestSuite given
      Truchas run-time specifications and a Loader
  
   Public Interface:
  
      RunPostProcessorTestSuite(RunBaseTestSuite)
           options(self)
  
   Contains:
   
      class RunPostProcessorTestSuite(RunBaseTestSuite)
            __init__(self,argv)
            options(self)
            
   Unit Test Block
  
   Author: Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""

import os, sys, string, fnmatch

thisdir     = os.path.abspath(os.path.dirname(__file__))
testingdir  = thisdir + '/../'

if __name__ == "__main__":
    sys.path.append(testingdir)
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

from RunBaseTestSuite import RunBaseTestSuite

class RunPostProcessorTestSuite(RunBaseTestSuite):
    """
    Wrapper that allows the user to load and run a TestSuite containing
    one or more of 
    TruchasBasicRunTest,
    TruchasPostProcessorRunTest,
    TruchasRestartRunTest,
    CheckPPGMV,
    CheckPPTecPlot,
    CheckPPEnsight,
    CheckPPVTK,
    CheckPPDiagnostics.
    """

    def __init__(self,
                 argv,
                 currdir,
                 debug = 0): 

        self.currdir      = currdir
	self.srcdir       = currdir
        truchas_bin       = testingdir + '/../../../bin'
        truchas_bin       = os.path.abspath(truchas_bin)
        self.parallel_env = None
        self.np           = 1
        self.compile_mode = None
        self.compiler     = None
        self.executables  = []
        self.exevalue     = None
	self.bindir       = truchas_bin
        self.no_rerun     = 0
        #for formatting output
        self.column       = 5
        self.column1      = 10
        self.column2      = 15
        self.debug        = debug

        self.setLogDir(outputdir = self.currdir)

    def options(self):
        
        self.defineUsage()
        self.usage        = self.baseusage.replace('RunTestSuite.py','RunPostProcessorTestSuite.py')        
        self.getOptions()
        self.parseOptions()
        self.checkOptions()

if __name__=="__main__":

    # use .inp files from TestCases
    # start by removing any .inp files in the curdir
    files = os.listdir('.')
    for f in files: 
	if fnmatch.fnmatch(f,'*.inp'): os.remove(f)
    # now copy .inp files from ../TestCases
    files = os.listdir('../TestCases')
    for f in files: 
	if fnmatch.fnmatch(f,'static_drop.inp'): os.link('../TestCases/'+f,f)

    from Loaders import TruchasPostProcessorLoader
    argv = sys.argv[1:]
    X    = RunPostProcessorTestSuite(argv,
                                     currdir = os.getcwd())
    X.options()

    outfile   = '%s_component_test' %(__file__)
    fpsummary = open(outfile,'w')
    
    for executable in X.executables:

        print
        print 'Executable : %s' %(str(executable))
        print        
        
        X.loadAndRun(loader     = TruchasPostProcessorLoader,
                     executable = executable,
                     fpsummary  = fpsummary)

    fpsummary.close()
    
    X.cleanUp()
    
