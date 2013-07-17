#!/usr/bin/env python
# TBrookParse.py

# this is the main driver script for the truchas output parser/translator

import os, sys, getopt

if (sys.version[0:3] == '2.6' or sys.version[0:3] == '2.7'):
    pass
else:
    print 'FATAL: TBrookParse.py requires Python 2.6 or 2.7'
    sys.exit(1)

# check if we have numpy.oldnumeric available - this is what we really want
try:
    import numpy.oldnumeric as Numeric
except ImportError:
    print
    print 'warning: this python installation doesn\'t have numpy installed'
    print '    attempting to use Numeric instead, although this is deprecated'
    print '    suggest installing numpy'

# no numpy.oldnumeric, check for real Numeric - this is deprecated - code should be converted to numpy
    try:
        import Numeric
    except ImportError:
        print
        print 'fatal: this python installation has neither numpy nor Numeric'
        print '    can\'t continue'
        print '    suggest installing numpy'
        print
        sys.exit(1)

# add directories to sys.path - directories are searched for importing modules
#     ${TRUCHAS_ROOT}/bin      location of compiled modules (shared objects)
#     ${TRUCHAS_ROOT}/tools    location of various python modules

scriptdir = os.path.abspath(os.path.dirname(sys.argv[0]))
bindir    = '@Truchas_BIN_INSTALL_DIR@'
tooldir   = '@TBrookParser_SOURCE_DIR@'
sys.path.append(bindir)
sys.path.append(tooldir)
# print 'sys.path:', sys.path

# check if we can import Truchas python modules - failure would be an internal error
# this is crude, but there may be sriramic effects involved with changing it and
# the * may walk all over our current namespace
try:
    from PythonPackages import *
except:
    print
    print 'fatal: unable to import Truchas Python modules'
    print '    internal error'
    print
    raise

# now that we knmow we can get at the PythonPackages, import some we need right now
from getPostProcessorObject import getPostProcessorObject
from processInput import processInput

# this is really ugly - why are we setting defaults here?  why do we need default files at all?
def_file    = os.path.abspath(os.path.dirname(__file__)) + '/test_TBrookParse/samples/phasechange_mixed_output/phasechange_mixed.TBrook.xml'
def_exofile = os.path.abspath(os.path.dirname(__file__)) + '/test_TBrookParse/samples/hytec_quarter_coarse_side_1.exo'

# instantiate a post processor
pp = getPostProcessorObject(def_file,def_exofile,input,debug=0)

description = """
          help:            display this message
          load:            load a TBrook file
          quit:            exit the postprocessor
          list:            list either loaded files, regions or variables
          region:          define a spatial region
          query:           query variables
          stat:            list statistics for a single variable
          write:           write out a visualization file
          restart:         write a restart file
          probe:           interrogate Truchas probe variables
          timeseries:      create a time series of a Truchas variable
          deleteregions:   delete all regions defined specifically by the user
          legacy:          provide legacy Truchas output files
          defaults:        set postprocessor defaults (precision)
"""            
    
commands    = {'load'            :pp.load,
               'help'            :pp.help,
               'h'               :pp.help,
               'quit'            :quit,
               'q'               :quit,
               'end'             :quit,
               'e'               :quit,
               'list'            :list,
               'define'          :pp.createregion,
               'region'          :pp.createregion,
               'query'           :pp.query,
               'stat'            :pp.stat,
               'statistics'      :pp.stat,
               'write'           :pp.write,
               'restart'         :pp.restart,
               'probe'           :pp.probe,
               'timeseries'      :pp.timeseries,
               'deleteAllRegions':pp.deleteregions,
               'deleteregions'   :pp.deleteregions,
               'defaults'        :pp.defaults,
               'description'     :description
               }

# process command line
inBuffer  = ''
opts, files = getopt.getopt(sys.argv[1:], 'hf:e:')

for o,a in opts:
    if (o == '-h'):
        pp.help(commands)
        sys.exit()
    elif (o == '-f'):
        inBuffer = '@' + a.rstrip().lstrip() + ' q'
    elif (o == '-e'):
        inBuffer = a.rstrip().lstrip() + ' ' + ' q'
        if len(files):
            inBuffer = a.rstrip().lstrip() + ' ' + files[0] + ' q'
            files = []

# run commands according to user input
#     seems like processInput should be a post processor method: pp.processInput()
T = processInput(inBuffer, sys.stdout, commands, files, def_file, def_exofile) 
