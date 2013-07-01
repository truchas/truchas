"""
Purpose:
  Example and default parameters file for input to MAPVariable.py
  internal unit testing.
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
USER:
 Modify the values below to define
 your mesh selections, source location and outputs
"""
if (True):
    " Provide a directory from which to obtain mesh and XML files."
    InputsDir = __Tdir + '/tools/scripts/test_TBrookParse/samples/'

    " Provide a TBrook.xml file prefix"
    TBprefix  = 'map'
    TBfile    = '%s_output/%s.TBrook.xml' %(TBprefix,TBprefix)
    
    " Define DICTs of A & B EXO files and mesh types (HEX/TET)"
    " An empty DICT for either file implies use mesh from XML file"
    mesha     = {}  
    meshb     = {'mesh-a_f.exo' : 'HEX'}

    "Provide a coordinate scale factor for use with mesh 'a'"
    CSF       = 1.0

    "Provide a mapping rule: 'constants_preserving' or 'conservative'"
    mapRule   = 'conservative'

