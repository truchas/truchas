#!/usr/bin/env python
"""
RunThisCapability

-----------------------------------------------------------------------------
   Purpose:
  
      Wrapper to allow developer to create and run a Capability TestSuite
      and clean the results of the Capability TestSuite.
  
   Public Interface:
  
      RunThisCapability(testdir,loadername,clean)
  
   Contains:
   
      class RunThisCapability
            __init__(testdir,loadername,clean) 
            
   Unit Test Block
  
   Author: Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""

import os, sys

class RunThisCapability:

    def __init__(self,testdir,
                 runnername = 'RunThisTestSuite', 
                 loadername = 'TruchasCapabilityLoader',
                 clean      = 0,
                 debug      = 0):

        self.testdir = testdir
        self.clean   = clean
        self.debug   = debug

        parserdir   = self.testdir + '/../TBrookParser'
        sys.path.append(parserdir)

        if clean:
            from Clean import Clean
            Clean()
        else:
            try:    
                loader     = 'Loaders.' + loadername
                __import__(loader)
            except ImportError:
                loader     = loadername
                __import__(loader)
            except:
                raise sys.exit(), 'Unable to import loader %s' %(loader)
            module         = sys.modules[loader]
            thisloader     = vars(module)[loadername]

            try:    
                runner     = 'Runners.' + runnername 
                __import__(runner)
            except ImportError:
                runner     =  runnername
                __import__(runner)
            except:
                raise sys.exit(), 'Unable to import runner %s' %(runner)
            module     = sys.modules[runner]
            thisrunner = vars(module)[runnername]
            
            argv  = sys.argv[1:]

            X     = thisrunner(argv,
                               srcdir  = os.getcwd(), #os.path.abspath(self.testdir + '/../../../regressiontests'),
                               currdir = os.getcwd(),
                               debug   = self.debug)
            X.options()
            
            fpsummary = open('/dev/null','w')
            
            if X.no_rerun:

                X.loadAndRun(loader     = thisloader,
                             executable = X.exevalue,
                             fpsummary  = fpsummary)
                
            else:
                
                for executable in X.executables:
                
                    X.loadAndRun(loader     = thisloader, 
                                 executable = executable,
                                 fpsummary  = fpsummary)
            
            fpsummary.close()
        
