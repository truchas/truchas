#!/usr/bin/env python
"""
RunTimeSpecs

-----------------------------------------------------------------------------
   Purpose:

      Class to hold and print run-time specifications employed in 
      a Truchas BasicRun or Truchas RestartRun
  
   Public Interface:
  
      RunTimeSpecs
           str()
  
   Contains:
   
      class RunTimeSpecs
            __init__(self,logdir,
                     logfile,
		     parallel_env,
		     np,
		     compiler,
		     compile_mode,
                     success,
                     version)
            str(self)
            
   Author: Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""

import os, string

class RunTimeSpecs:

    def __init__(self,logdir,
                 logfile,
		 parallel_env = 'serial',
		 np           = 1,
		 compiler     = 'lahey',
		 compile_mode = 'opt',
                 version      = '2.3',
                 success      = 'truchas terminated normally'):

	self.logdir       = logdir
        self.logfile      = logfile
	self.parallel_env = parallel_env
	self.np           = np
	self.compiler     = compiler
	self.compile_mode = compile_mode
        self.version      = version
        self.basename     = None
        self.success      = success

        #extract problem basename from logfile...

        if self.logdir in os.getcwd():
            for file in os.listdir(self.logdir):
                if '_logs' in file:
                    L2            = string.split(file,'_logs')
                    self.basename = L2[0]

    def str(self,column=10):

	tmp   = 'Parallel env  : ' + self.parallel_env
	strng = tmp.rjust(column + len(tmp))
	info  = '\n'
	info += strng

	tmp   = 'np            : ' + str(self.np)
	strng = tmp.rjust(column + len(tmp))
	info += '\n'
	info += strng

	tmp   = 'Compiler      : ' + self.compiler
	strng = tmp.rjust(column + len(tmp))
	info += '\n'
	info += strng

	tmp   = 'Compile mode  : ' + self.compile_mode
	strng = tmp.rjust(column + len(tmp))
	info += '\n'
	info += strng

	tmp   = 'Version       : ' + self.version
	strng = tmp.rjust(column + len(tmp))
	info += '\n'
	info += strng

        L  = string.split(self.logfile,'/')
        i  = 0
        for word in L:
            if 'outputs.t' in word: i=L.index(word)
        thelogfile = string.join(L[i:],'/')
        
	tmp   = 'Log file      : ' + thelogfile
	strng = tmp.rjust(column + len(tmp))
	info += '\n'
	info += strng
	
	return info    
        
