#!/usr/bin/env python
"""
 TestingInfrastructure

-----------------------------------------------------------------------------
   Purpose:
  
      provides helper Python modules for developers writing test scripts
  
   Public Interface:
  
      import TestingInfrastructure
  
   Contains: Loaders      (loads a set of TestCases to create a TestSuite)
             TestCases    (a type of Truchas test; is one of
                           Basic Run, Restart Run, PostProcessor Run, Capability) 
             Runners      (runs a TestSuite)
             DataHandlers (creates, manipulates data storages generated from
                           postprocessing Truchas .TBrook.xml files)
             Loggers      (Python logging module, used for debugging)
  
   Author : Sharen Cummins (scummins@lanl.gov)
 -----------------------------------------------------------------------------
"""

import os, sys
thisdir     = os.path.abspath(os.path.dirname(__file__))
sys.path.append(thisdir)

import Loggers
import Runners
import TestCases
import Loaders
import DataHandlers
import Clean
