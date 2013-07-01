#!/usr/bin/env python
"""
 Clean

-----------------------------------------------------------------------------
   Purpose:
  
      Clean up after component tests
  
   Public Interface:
  
      Clean
  
   Contains: Clean

   Unit Test Block
  
   Author: Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""

import os,fnmatch

class Clean:

    def __init__(self):
        
        filepatterns = ['*.out', '*_output', '*.ppout', 'gmv.bin', 'outputs.*','for_debugging','debug.*','*.txt','*.ascii','stderr_*','*.log']

        print '\n\n        Cleaning up all outputs from any previous runs...  \n'

        files = os.listdir(os.getcwd())
        for file in files:
            for pattern in filepatterns:
                if fnmatch.fnmatch(file,pattern):
                    if file in files:
                        try:
                            os.remove(file)
                        except:
                            cmd = 'rm -r %s ' %(file)
                            os.system(cmd)
                    files = os.listdir(os.getcwd())
 
        print '\n                     Done cleaning......                        \n'

if __name__ == "__main__":
    Clean()
