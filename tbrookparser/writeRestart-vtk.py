"""
writeRestart
----------------------------------------------------------------------------
  Purpose:
    Stand-alone script to generate a Truchas restart file.  This script uses 
    the underlying post-processor utilities but not the parser intferace.

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

		
def restartNoIO(Xs, options, input, debug=False):

	if debug: print "Debug > welcome to restartNoIO"
        meshchs      = []
        ignoredfiles = []
	# mmove verification to main or seperate functions
	if debug: print "Debug > verifying ",options['meshes']," exists"
        for X in Xs:
            m      = []   
            meshes = X.storage.mlist
            if (len(meshes)):
                meshlist = ''
                for mesh in meshes:
                    m += [str(mesh.name)]
                for mesh in options['meshes']:
		    if mesh in m: meshchs.append(mesh)
            else:
                print 
                print '     No meshes in the simulation: %s' %(X.file)
                print '     Will not create a restart for this simulation'
                print 
                ignoredfiles.append(Xs.index(X))
	if debug: print "Debug > meshchs:",meshchs

        T  = []
        if len(ignoredfiles)> 0 and len(Xs) > 0:
            for ind in ignoredfiles: T.append(Xs[ind])
            if len(T) == len(Xs): return
            for val in T:
                if val in Xs: del Xs[Xs.index(val)]

	ids = {}
	cycles = {}
	cnt = 0
	ignoredfiles = []

	for X in Xs:
	    ignorfile = 0
	    #ids[cnt],ignorfile = chooseRestartSteps(X,type)
	    steps = X.storage.tlist
	    c = []
	    if (len(steps)):
		for step in steps:
		    c += [step.cycle]
	    if options['cycle'][X.file] in c:
                ids[cnt] = [c.index(options['cycle'][X.file])]
		if debug: print "Degbu > confirming cycle ",ids[cnt]
	    else: 
	        print "Error - cycle",options['cycle'],"not found in ",X.file
		print "Please try again with a cycle from:"
		print c
		sys.exit(2)
            if len(c) == 0:
		ignoredfiles.append(Xs.index(X))
	    else:
		cnt += 1
        T  = []
        if len(ignoredfiles)> 0 and len(Xs) > 0:
            for ind in ignoredfiles: T.append(Xs[ind])
            if len(T) == len(Xs): return
            for val in T:
                if val in Xs: del Xs[Xs.index(val)]
	


        these_storages = []
        these_files    = []
        these_meshes   = []
	files_storages = {}
        
        for X in Xs:
            these_storages.append(X.storage)
            these_files.append(X.file)
	    files_storages[X.file] = [X.storage]
            for mesh in X.storage.mlist:
                for meshch in meshchs:
                    if meshch == mesh.name:
                        these_meshes.append(mesh)
                            
        if options['mapping']:
	    if debug: print "Debug > doing mapping"
            try:
                #check for any existing exodus files
                exof    =  options['exoFile']
                for thisfile in os.listdir(os.getcwd()):
                    if thisfile[-4:] == '.exo':
                        exof = thisfile

                stepids = []
                for i in range(len(these_files)):
                    stepids.append(ids[i][0])

		mapcount = 1
	        if debug: print "Debug > calling getFileObject with ",exof
		  
                Y  = getFileObject(input,
				   options['exoFile'],
                                   these_files, 
				   files_storages,
                                   stepids, 
				   mapcount,
                                   debug = debug,
				   options = options)

		if debug: print "Debug > removing features: ",Y.storage.specs.feats
                """
                Since we have chosen to remap we must now alter the storage 
                 restart features so that the features:
                      1.temperature
                      2.solid_mechanics
                      3.fluid_flow
                      4.joule_heat
                 are ignored in the restart
                (i.e will not be read in in the newly mapped restart file) 
                """
                feats  = Y.storage.specs.feats
		if debug: print feats

                for feat in feats:
                    if 'temperature' in feat:
			if debug: print 'Debug - temperature features: ',feats[feats.index(feat)]
                        del feats[feats.index(feat)]
                for feat in feats:
                    if 'fluid_flow' in feat:
			if debug: print 'Debug - fluid_flow features: ',feats[feats.index(feat)]
                        del feats[feats.index(feat)]
                for feat in feats:
                    if 'joule_heat' in feat:
			if debug: print 'Debug - joule_heat features: ',feats[feats.index(feat)]
                        del feats[feats.index(feat)]
                for feat in feats:
                    if 'solid_mechanics' in feat:
			if debug: print 'Debug - solid_mechanices features: ',feats[feats.index(feat)]
                        del feats[feats.index(feat)]
                        
                Y.storage.specs.feats  = feats
                Y.storage.specs.nfeats = len(feats)

		if debug: print "Debug > filling these_*"
                these_files       = []
                these_storages    = []
                these_meshes      = []
                these_files.append(Y.file)
                these_storages.append(Y.storage)
                these_meshes.append(Y.storage.mlist[0])
                # if we have chosen to map to a new mesh,
                # adjust the cycles list and id list
		if debug: print "Debug > doing some cycle stuff"
                cnt          = 0
                cycles[cnt]  = Y.storage.timesteps.cyclelist

                ids[cnt]     = []
                ids[cnt].append(Y.storage.timesteps.idlist[-1])
                if len(cycles[cnt]) > 1 :
                    print 
                    print '     In this postprocessing event you have performed'
                    print '     %i mappings (i.e cycles %s) with these mesh',\
			  'choices.' %(len(cycles[cnt]), str(cycles[cnt]))
                    print 
                    cycles[cnt][0] = cycles[cnt][-1] 
                    ids[cnt]       = [cycles[cnt][0]]
                    
                #now save this file object in the post-processor's 
		#storage structure
		if debug: print "Debug > doing some file stuff"
                files_storages[Y.file] = [Y.storage]

            except:
                print 
                print 'Error in obtaining new restart file from %s '%(exofile)
                print 'Aborting restart writer.'
                print 
                print 'Please ensure your original file:'
                print '  1.is an XML file'
                print 'Please ensure your new mesh: '
                print '  1. Is a non-degenerate HEX Mesh'
                print '  2. Overlaps with the Truchas source mesh'
                print '  3. Is in Exodus II format'
                print 


        # generate restart files now
	if debug: print "Debug > Creating all restart files now.."
        
        for this_file in these_files:
            cnt     = these_files.index(this_file)
            steps   = these_storages[cnt].tlist            
            iend    = string.find(this_file,".TBrook.xml")
            thefile = this_file[0:iend]
            
	    if debug: print "Debug > looping over ids:",ids[cnt]
            for idx in ids[cnt]:
		if debug: print "  Debug > on idx = ",idx
                # get the values for this timestep
                these_storages[cnt].getValues(steps[idx].vlist,options)
                from RESTARTwriteStorageObject import RESTARTwriteStorageObject
		if debug: print "Debug > calling RESTARTwriteStorageObject"
                RESTARTwriteStorageObject(options['outFile'], 
				          options['outFmt'],
				          these_storages[cnt],
                                          these_meshes[cnt], 
					  idx,
                                          steps[idx].time, 
					  steps[idx].vlist,
                                          dformat = '%20.5E',
                                          debug   = debug)
		mesh = these_storages[cnt].mlist[0]
		mesh.fillMesh()

		if debug: print "Debug > calling GMVwriteSO"

                from VTKwriteStorageObject import VTKwriteStorageObject
                VTKwriteStorageObject(options['vtkFile'],
				          options['outFmt'],
					  these_meshes[cnt],
					  these_storages[cnt].tlist[idx].time,
					  idx,
					  these_storages[cnt].tlist[idx].vlist,
					  fpwatch = sys.stdout,
					  debug = debug)







def usage():

    print
    print "Current usage:"
    print "  writeRestart.py -i <inFile1,inFile2,...> -o <outFile> "
    print 
    print " available options:"
    print "    -i: input file(s) [phasechange_mixed.TBrook.xml]"
    print "    -o: output file [restart.bin]"
    print "    -S: indicates a standard restart [default]"
    print "    -M: indicates mapped restart"
    print "    -e: destination Exodus II mesh file *ONLY FOR MAPPING*"
    print "        [hytec_quarter_coarse_side_1.exo]"
    print "    -c: cycle to create restart from [last cycle]"
    print "    -m: which mesh to use [DefaultMesh]"
    print "    -s: mesh scale factor [defaulted to exo meshfile value]"
    print "    -f: file format ascii/binary [binary]"
    print "    -d: sets debug flag to True "
    print "    -v: provide a list of variables and default values"
    print "        This options overrides the default values assigned"
    print "        to variables during the mapping process."
    print "        <varname1:val,varname2:val,...>"
    print "        NOTE: follow format carefully for this options"
    print "    -b: provide block info for mapping"
    print "        <0:blockid1:blockid2:...,1:blockid3:...>"
    print "         0,1,... coorspond to listing of input files"
    print "    -x: provide override info for mapping"
    print "        <varname:blockid:value,varname2:blockid:value,...>"
    print
    print " NOTE: Follow formating for -v and -b options carefully"
    print " NOTE: Option flags are case sensitive"
    print " NOTE: If overriding VOF* values, be sure to override all "
    print "       desired VOFs.  "
    print " NOTE: By default all VOFs are set to 0.0 in mesh blocks not"
    print "       in the source mesh(s).  The user must set VOF values"
    print "       appropriately.  Override the destination mesh values"
    print "       using the -x option."
    print " NOTE: This routine does not check for valid variable values."
    print "       That is up to the user."
    print 
    print " Default variable mapping is listed in the file "
    print "    truchas/tools/PythonPackages/TBrookParser/MODutils/MAPutils/defaultmaps.txt"
    print
    return None


if __name__=='__main__':


    # set up some defaults to get started
    def_file    = os.path.abspath(os.path.dirname(__file__)) + '/test_TBrookParse/samples/phasechange_mixed_output/phasechange_mixed.TBrook.xml'
    def_exofile = os.path.abspath(os.path.dirname(__file__)) + 'test_TBrookParse/samples/hytec_quarter_coarse_side_1.exo'
    def_outfile = 'restart.bin'
    def_vtkfile = 'restart.vtk'
    inbuffer = ''
    input = usubs.input(inbuffer)
    debug = False
    input = None

    # process cmd line options
    try:
	opts, args = getopt.getopt(sys.argv[1:], "i:e:o:c:f:m:s:v:b:x:SMdh",
			           ["blocks","help"])
    except getopt.GetoptError,err:
	print " Sorry, encountered the following error"
	print str(err)
	print 
	#usage()
	sys.exit(2)

    for o, a in opts:
	if o == "-d":
	    debug = True
	if o == "-f":
	    if a.lower() in ('ascii'): 
		def_outfile = 'restart.ascii'

    if debug: print "Debug > setting up defaults"
    Xs      = []
    options = {'inFiles' :[],
	       'outFile' :def_outfile,
	       'exoFile' :def_exofile, 
	       'vtkFile' :def_vtkfile,
	       'mapping' :0, 
	       'cycle'   :{}, 
	       'meshes'  :['DefaultMesh'], 
	       'scale'   :None,
	       'outFmt'  :'binary',
	       'vars'    :{},
	       'blocks'  :{},
	       'override':{}} 
    if debug: 
	print "Debug > set options based on cmd line input"
	print "Debug > ",opts
    for o, a in opts:
	if o == "-i":
	    files = string.splitfields(a,',')
	    options['inFiles'] = files
	if o == "-e":
	    options['exoFile'] = a
	if o == "-o":
	    options['outFile'] = a
	    fname = a
	    if fname[-4:] == '.bin': fname = fname[:-4]
	    options['vtkFile'] = fname + '.vtk'
	if o == "-S":
	    options['mapping'] = 0
	if o == "-M":
	    options['mapping'] = 1
	if o == "-m":
	    options['meshes'] = string.splitfields(a,',')
	if o == "-c":
  	    cycles = string.splitfields(a,',')
            if len(cycles) > 1: 
                cnt = 0
	        for c in cycles:
	            options['cycle'][options['inFiles'][cnt]] = int(c)
	            cnt += 1
            else:
	        options['cycle'][options['inFiles'][0]] = int(cycles[0])
	if o == "-s":
	    options['scale'] = float(a)
	if o == "-f":
	    options['outFmt'] = a
	    if a.lower() not in ('binary','ascii'):
		print "  Invalid output file format of ",a
		print "  Please choose binary or ascii and try again"
		sys.exit()
	if o in ("-v"):
	    sets = string.splitfields(a,',')
	    for item in sets:
	        name,val = string.splitfields(item,':')
	        options['vars'][name] = float(val)
	if o in ("-b"):
	    # parse the input string to build blocks dictionary
	    inp = string.splitfields(a,',')
	    for i in inp:
	        set = string.splitfields(i,':')
		id = int(set[0])
		options['blocks'][id] = []
		for b in set[1:]: 
		    options['blocks'][id].append(int(b))
	    # build a reverse dictionary, remove any entires with < 2 entries
	    if len(options['blocks']) > 1:
		multi = {}
		for inp in options['blocks']:
		    for bid in options['blocks'][inp]:
			if bid not in multi: multi[bid] = [inp]
			else: multi[bid].append(inp)
		for id in multi.keys():
		    if len(multi[id]) < 2: del multi[id]
		if len(multi) > 1: options['multi'] = multi
	if o in ("-x"):
	    vars = string.splitfields(a,',')
	    for v in vars:
		items = string.splitfields(v,':')
		try:
		    varname = items[0]
		    blockid = int(items[1])
		    value   = float(items[2])
		    if varname not in options['override'].keys():
		        options['override'][varname] = [(blockid,value)]
		    else:
		        options['override'][varname].append((blockid,value))
		except:
		    print "    Problem reading override option"
		    usage()
		    sys.exit()

	if o in ("-h","-help","--help"):
	    usage()
	    sys.exit()

    # set default inFile if none provided by user
    if len(options['inFiles']) < 1: 
	options['inFiles'].append(def_file)

    if debug: 
	print "Debug > options:"
	for key in options.keys(): print "    ",key,":",options[key]
	    
    try:
	if debug: print "Debug > init getPPO with input =",input
	ppobject = getPostProcessorObject(options['inFiles'][0],
			                  options['exoFile'],
		                          input,fp=sys.stderr,debug=debug)
	ppobject.noIO = 0
	ppobject.options = options
	for file in options['inFiles']:
            if file not in ppobject.files:
	        if debug: print "Debug > load fileObj ",file
		X = ppobject.load(file)
	        if debug: print "Debug > done load fileObj ",file
		# if no cycle number input, default to last cycle
		if (file not in options['cycle'].keys()) or \
		   (options['cycle'][file]==None):
		    options['cycle'][file] = X.storage.tlist[-1].cycle
		    if debug: 
			print "Debug > just defaulted cycle to ",\
			       options['cycle'][file]
		elif (debug): 
		    print "Debug > using user specified cycle",
		    print options['cycle'][file]
		# add section to gather a cycle list and check if cycle exists
		#if options['cycle'][file] not in X.storage.tlist:
		#    print "Error - Cycle ",options['cycle'][file],\
		#          "not found in file",file
		#    sys.exit(2)
		Xs.append(X)
	if debug: print "Debug > calling restart"
	restartNoIO(Xs, options, input, debug)
    except:
	print "--> Test failed <--"
	print "Nature-Of-Error:", sys.exc_info()[0]

