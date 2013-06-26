"""
Purpose:
  Parameters file two extract meshes from 2 Truchas simulations
  and output them in ascii form. Required to test the map_wrapper.F90
  file. 
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
    "Provide a directory from which to obtain mesh and XML files."
    InputsDir    = __Tdir + '/tools/scripts/test_TBrookParse/samples/'

    "Provide two Truchas input file prefixs"
    TBprfx       = ['map', 'map-b']
    
    "Conversion to filenames relative to InputsDir"
    TBfile       = '%s_output/%s.TBrook.xml' %(TBprfx[0],TBprfx[0])
    TBfile2      = '%s_output/%s.TBrook.xml' %(TBprfx[1],TBprfx[1])

    mesha        = {}
    meshb        = {}

    "Coordinate scale factor"
    CSF          = None
    
    "Provide a mapping rule: 'constants_preserving' or 'conservative'"
    mapRule      = 'conservative'
