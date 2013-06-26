"""
Purpose:
  Template Parameters file for input to test_mapping.py, which is
  and external unit-test driver for MAPVariable.py
"""

"---------------------------------------------------------------------"
" !!! Do not change this directory magic !!! "
import os
if ('TRUCHAS' in os.environ.keys()):
    "pick up TRUCHAS home directory from environment by default"
    __Tdir = os.environ['TRUCHAS']
else:
    "otherwise, <assume> a relative location to this file."
    __cwd  = os.path.dirname(__file__)
    __Tdir = os.path.abspath(__cwd + '../../../../../../')
"---------------------------------------------------------------------"

"""
USER: Modify the values below to define
      your mesh selections, source location and outputs
"""
if (True):
    " Provide a directory from which to obtain mesh and XML files."
    InputsDir    = __Tdir + 'tools/scripts/test_TBrookParse/samples/'

    " Define DICTs of A & B EXO files and mesh types (HEX/TET)"
    mesha_files  = {}
    
    meshb_files  = {}
    
    " Provide an output filename (written to current working directory)"
    filename     = ''
    
    " Provide a graphics filename prefix"
    graphicsname = ''

