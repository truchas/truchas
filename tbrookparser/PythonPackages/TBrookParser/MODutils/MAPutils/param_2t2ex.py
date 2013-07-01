"""
Purpose:
  Example and default parameters file for input to test_2t_to_ex.py,
  which is and external unit-test driver for MAPVariable.py
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

    "Provide a directory name for output files"
    OutputDir    = 'mapa_and_mapb_to_meshab' 
    
    "Provide two Truchas input file prefixs"
    TBprfx       = ['map', 'map-b']
    
    "Conversion to filenames relative to InputsDir"
    TBfile       = '%s_output/%s.TBrook.xml' %(TBprfx[0],TBprfx[0])
    TBfile2      = '%s_output/%s.TBrook.xml' %(TBprfx[1],TBprfx[1])
    T_Sim        = '%s_output/%s.TBrook.xml' %(TBprfx[0],TBprfx[0])
    T_Sim2       = '%s_output/%s.TBrook.xml' %(TBprfx[1],TBprfx[1])

    "Coordinate scale factor"
    CSF          = None
    
    "Cycles for each TBrook file"
    CYCLE        = (-1,-1) # use -1 to get the last cycle number

    "Define DICTs of A & B EXO files and mesh types (HEX/TET)"
    EXOfile      = {'mesh-ab_f.exo' : 'HEX'}

    "Provide a mapping rule: 'constants_preserving' or 'conservative'"
    mapRule   = 'conservative'
