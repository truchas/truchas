"""
SimInfo
----------------------------------------------------------------------------
  Purpose:
    Stand-alone script to gather information about a Truchas simulation run.
    This script uses the underlying post-processor utilities but not the 
    parser intferace.

    Information extracted:
        Name
	Cycles/Times
	Variables
	Blocks
	Probes
	Meshes
	Scale

  Author: Erin Iesulauro Barker (eibarker@lanl.gov)
----------------------------------------------------------------------------
"""

import os, sys, string, getopt

currdir   = os.path.abspath(os.path.dirname(sys.argv[0]))
bindir    = os.path.normpath(currdir + '/../../bin')
tooldir   = os.path.normpath(currdir + '/../tools')
tbrookdir = os.path.normpath(tooldir + '/../PythonPackages/TBrookParser')

sys.path.append(bindir)
sys.path.append(tooldir)
sys.path.append(tbrookdir)

from getPostProcessorObject import getPostProcessorObject
from XMLgetStorageObject    import getStorageObject
from POSTPROCESSORutils     import usubs,getFileObject
from XMLutils               import getSimSpecs
from PYTHONutils            import unique

def SimInfo(File, debug=False):

    if debug: print "Debug - welcome to SimInfo"
    # set variables
    def_exofile = os.path.abspath(os.path.dirname(__file__)) + 'test_TBrookParse/samples/hytec_quarter_coarse_side_1.exo'
    input = usubs.input('')
    options = {'inFiles':[File],
	       'outFile':None,
	       'exoFile':def_exofile,
	       'mapping':0,
	       'cycle'  :{},
	       'meshes' :[],
	       'scale'  :None,
	       'outFmt' :'binary',
	       'vars'   :{},
	       'blocks' :{} }
    
    # load provided xml file
    ppobject = getPostProcessorObject(File,def_exofile,input,debug=debug)
    X        = ppobject.load(File)
    options['cycle'][File] = X.storage.tlist[-1].cycle


    # extract avaiable timesteps
    if debug: print "Debug - extract available timesteps"
    cycles = []
    times  = []
    ids    = {}
    for step in X.storage.tlist:
	cycles += [step.cycle]
	times  += [step.time]
    if options['cycle'][File] in cycles:
	ids[0] = [cycles.index(options['cycle'][File])]

    # extract available mesh names
    if debug: print "Debug - extract mesh names"
    meshes = []
    for mesh in X.storage.mlist:
	meshes += [str(mesh.name)]

    # extract available features
    if debug: print "Debug - extract features"
    these_files = [X.file]
    files_storages = {}
    files_storages[X.file] = [X.storage]
    mapcount = 0
    stepids = []
    for i in range(len(these_files)): stepids.append(ids[i][0])

    if debug: print "Debug - call getFileObject"
    Y = getFileObject(input, options['inFiles'][0], these_files, 
	              files_storages, stepids, mapcount, sys.stdout, debug, 
		      options)

    if debug: print "Debug - extract regions"
    regions = []
    regname = []
    for reg in Y.regions: 
	regions.append(reg)
	regname.append(reg.name)

    print " Available Simulation Information:"
    print "    Meshes    : ",meshes
    print "    Cycles    : ",cycles
    print "    Times     : ",times
    print "    Features  : ",X.storage.specs.feats
    print "    Phases    : ",X.storage.specs.nphases
    print "    Components: ",X.storage.specs.ncomps
    print "    Scale     : ",X.storage.specs.csf
    print "    Regions   : ",regname
    print 

    for mesh in X.storage.mlist: 
        print "    Mesh:",mesh.name
	print "            ",mesh.type
	for sp in mesh.mslist:
	    print "            ",sp.name,":",sp.size
	    for var in sp.vlist: 
		#print var.str()
	        if var.name == 'SIDESETS':
	            print "             NUMSIDESETS:",var.rank
        if mesh.cells.has_key('BLOCKID'):
	    print "             BLOCKID:",unique(mesh.cells['BLOCKID'])
	    

    return None

def usage():

    print
    print "Current usage:"
    print "  SimInfo -d -h TBrookFile"
    print
    print " available options:"
    print "    -d:        sets debug flag to True"
    print "    -h/--help: prints usage message and exits"
    print

if __name__=='__main__':

    debug = False
    def_file    = os.path.abspath(os.path.dirname(__file__)) + '/test_TBrookParse/samples/phasechange_mixed_output/phasechange_mixed.TBrook.xml'

    try:
	opts, args = getopt.getopt(sys.argv[1:], "dh",["help"])

	for o, a in opts:
            if o == "-d": 
		debug = True
	    if o in ("-h","-help","--help"):
    	        usage()
	        sys.exit()
	
	if debug: 
	    print "Debug - opts: ",opts
	    print "Debug - args: ",args

	if len(args) == 0: 
	    filename = def_file
	else:              
	    #filename = os.path.abspath(os.path.dirname(__file__)) + args[0]
	    filename =  args[0]

	if debug: print "Debug - filename: ",filename

        SimInfo(filename, debug )

    except getopt.GetoptError,err:
	print " Sorry, encountered the following error"
	print str(err)
	sys.exit(2)
