

"""

 MODStorageObject

 -----------------------------------------------------------------------------
  Purpose:

    Instantiates and loads the input data from the mesh and simulation files.

  Public Interface(s):

    mSO = modStorageObject(TSObjects,TMeshs)
    mSO.modThisVar(var)
    mSO.getValues(vars)
    mSO.str()

  Contains:
    class modStorageObject
        __init__(self,TSObjects,TMeshs,CurrSObject=None,InpMesh='',
                 meshscale_factor=1.0,stepids=[0],debug=0):
        modThisVar(self,var)
        getValues(self,vars)
        str(self)

    Unit Test Block

  Version:
    $ID$

  Author(s): Sharen Cummins (scummins@lanl.gov)
 -----------------------------------------------------------------------------
"""

if __name__=='__main__':
    " Set sys.path for component testing mode "
    import os, sys
    thisdir   = os.path.dirname(__file__)
    # modify dot-dots appropriately
    parsedir  = os.path.abspath(thisdir+'../') 

    sys.path.append(parsedir)

from MODutils    import modSimSpecs,modMesh,modTimeSteps,\
                        modFieldVariable,createVariable
from PYTHONutils import getdir
import sys

class modStorageObject:

    def __init__(self,TSObjects,TMeshs,CurrSObject=None,InpMesh='',
                 meshscale_factor=1.0,stepids=[0],fp=sys.stdout,debug=0):

        #initialise the storage object

        self.specs          = None    #the simulation specifications
        self.mlist          = []      #a set of meshes
        self.tlist          = []      #a set of timesteps
        self.varlist        = []      #a set of user specified variables
        self.plist          = []      #a set of probes (yet to be implemented)
        if CurrSObject == None:
            self.timesteps  = None    #a timesteps object
                                      #(contains id lists,
                                      #cycle lists, times list
                                      #for the whole simulation)
        else:
            self.timesteps  = CurrSObject.timesteps

        self.spacelu        = TSObjects[0].spacelu
        self.namelu         = TSObjects[0].namelu

        #now define any extra attributes required for this particular class

        self.inpmesh        = InpMesh #the file name of the inputted mesh
        self.tstorages      = TSObjects
        self.tmeshs         = TMeshs
        self.stepids        = stepids
        self.debug          = debug
        self.fp             = fp   #file pointer to direct screen output to

        #ensure input mesh file specified is an exodus file
        import string
        try:
            correctfile =    string.find(self.inpmesh,".exo") > -1 \
                          or string.find(self.inpmesh,".EXO") > -1 \
                          or string.find(self.inpmesh,".gen") > -1 \
                          or string.find(self.inpmesh,".GEN") > -1
            assert correctfile
        except:
            print >> self.fp
            print >> self.fp, 'ERROR: Target mesh files must be end in either',\
                  '.exo, .EXO, .gen, .GEN extensions'
            print >> self.fp, 'Will not continue to map onto the target mesh',\
                  'as it is not in Exodus II format.'
            print >> self.fp
            sys.exit(1)

        """
        ensure numbers of inputted T storage objects,
        T meshes and stepids are consistent
        """
        try:
            assert len(stepids)   == len(TSObjects)
        except:
            print >> self.fp
            print >> self.fp, 'Number of ids inputted (%i) does not ',\
			      'correspond to number of Truchas simulations',\
			      'inputted (%i)' %(len(stepids),len(TSObjects))
            sys.exit(1)

        try:
            assert len(TSObjects) == len(TMeshs)
        except:
            print >> self.fp
            print >> self.fp, 'Number of meshes inputted (%i) does not ', \
			      'correspond to number of Truchas simulations ',\
			      'inputted (%i)' %(len(TMeshs),len(TSObjects))
            sys.exit(1)

        """
        get the simulation specifications
        - these come directly from the TSObjects
        - for multiple inputted T simulations we'll assume that the
          simulations specifications come from the first T simulation

        LJCox Note: this next line appears to be essentially identical to
        > self.specs    = self.tstorages[0].specs
        as modSimSpecs just copies the data from its argument.
        """
        self.specs    = modSimSpecs(self.tstorages[0].specs)

        if (self.debug > 0):
            print >> self.fp
            print >> self.fp, 'simulation specs'
            print >> self.fp, self.specs.str()
            print >> self.fp

        #obtain relevant data on all tmeshes
        msnames         = []
        for space in self.tmeshs[0].mslist:
            msnames.append(space.name)

        cnt             = 0
        for mesh in self.tmeshs:
            for space in mesh.mslist:
                self.tstorages[cnt].getValues(space.vlist)
            cnt         = cnt + 1

        #now form the new Exodus mesh where the data will be mapped too

        wdir        = getdir(__file__)
        meshgendir  = wdir + '/MODutils/EXOutils'
        mesh        = modMesh(self.tmeshs[0].type,msnames,self.inpmesh,
                              meshgendir,scale_factor=meshscale_factor,
                              clean='clean',fp=self.fp,debug=self.debug)
        self.mlist.append(mesh)

        #get all relevant timesteps

        self.timesteps      = modTimeSteps(self.tstorages,
                                           self.timesteps,stepids)
        self.tlist          = self.timesteps.tlist

    def modThisVar(self,var):

        #check to decide if we should modify variable x or not

        tmp  =     (not 'MU'            in var.name) \
               and (not 'SIGMA'         in var.name) \
               and (not 'COIL'          in var.name) \
               and (not 'sens_function' in var.name)

        return tmp

    def getValues(self,vars,options=None):
        """
        Provides user data values for a list of variables [vars]
        - we seperate this from obtaining other variable attributes so
          that the data values are obtained only when the user requests
          them

        Then mesh has been specified
        - use hex-to-hex mapping where appropriate
        - two ways of generating data depending on what variables
          are passed through :
           1. mesh variables  - already obtained from reading in mesh file
           2. field variables - obtained here from hex to hex mapping
        """

        try:
           import numpy.oldnumeric as Numeric
        except ImportError:
           import Numeric
        except:
           raise

        # first get variable data from all inputted Truchas simulations
        tmpstor = {}
        varstor = {}
        cnt     = 0

        for storage in self.tstorages:
            thisid       = self.stepids[cnt]
            tmplist      = storage.tlist[thisid].vlist
            tmpstor[cnt] = []
            for var in vars:
                for tmpvar in tmplist:
                    if var.name == tmpvar.name:
                        tmpstor[cnt].append(tmpvar)
            cnt = cnt + 1

        newvars = []

        cnt     = 0
        tmpdict = {0:'1st',1:'2nd'}

        """
        SJC Note: For now provide default value for all target variables
                  on regions outside the source mesh domain
        """
        defaultval = 1.0
	

	"""
	EIB Note: The following allows for the user to override the default
	          values being assigned to specific variables.  This follows
		  from the one-off version of the post-processor produced
		  for the Basi Hemi Problem (a.k.a. the opencalc) created by
		  Sharen Cummins.
	"""
	varvaldict = {}
	if options :
	    if options['vars']:
	        for key in options['vars'].keys():
	            varvaldict[key] = options['vars'][key]

        for storage in self.tstorages:
            varstor[cnt] = []
            for var in vars:
                thisvar = var
		if self.debug: print 'Debug - check ',var.name,' for modify'
                if var.data == None and var.rank > 0 and self.modThisVar(var):

                    print >> self.fp, '\n','*'*78
                    print >> self.fp, '    About to modify field variable',\
                          '%s from the %s Truchas simulation :\n' \
                          %(var.name, tmpdict[cnt])

		    """
		    EIB Note: Again, addition for default value override.
		    """
		    if 'VOF' in var.name: defaultval = 0.0
		    if varvaldict.has_key(var.name):
			defaultval = varvaldict[var.name]

                    for tvar in tmpstor[cnt]:
                        if var.name == tvar.name:
                            tvarl       = []
                            tvarl.append(tvar)
                            tvarl       = storage.getValues(tvarl)
                            #modify this Truchas field variable
                            wdir       = getdir(__file__)
                            mapvardir  = wdir + '/MODutils/MAPutils'
                            thisvar    = modFieldVariable(tvarl[0], 
					                  self.tmeshs[cnt],
                                                          self.mlist[0], 
							  mapvardir,
                                                          defaultval,
							  fp=self.fp,
							  debug=self.debug)
                    print >> self.fp, '\n    Finished modifying field ',\
				      'variable %s' %(var.name)
                    print >> self.fp, '*'*78,'\n'

                varstor[cnt].append(thisvar)

            cnt = cnt + 1

        """
        Now add data from the mapping(s) together to form a new variable
        For multiple mappings, subtract defaultval from the
        entire target variable data
        """
        for var in vars:

            if var.data == None:

                cnt       = 0
                thisindex = vars.index(var)
                thisvar   = varstor[0][thisindex]
		data      = None

                for storage in self.tstorages:
		    thedata = varstor[cnt][thisindex].data
		    if thedata != None: 
		        blockids = self.mlist[0].cells['BLOCKID']
		        try:
		    	    bids = options['blocks'][cnt]
			    cond = blockids == bids[0]
			    if len(bids) > 1:
			        for i in xrange(1,len(bids)):
			            cond = cond + (blockids == bids[i])
			    d = Numeric.where(cond,thedata,defaultval)
			    try:
			        data = d + data
			    except:
			        data = d
		        except:
			    try:
			        data = data + thedata
			    except:
			        data = thedata

                    cnt = cnt + 1

                vars[thisindex] = createVariable(thisvar.name, thisvar.nickname,
                                                 thisvar.rank, thisvar.shape,
                                                 thisvar.type, data,
                                                 self.mlist[0].name, 
						 thisvar.meshspace)

	# override values in the final arrays, specified by user
	if options != None:
	    if options['override'] != None:
		# get list of blockids!!
		for var in vars:
		    if var.name in options['override']:
			for block in options['override'][var.name]:
  			    olddata = var.data
		            blockids = self.mlist[0].cells['BLOCKID']
			    cond = (blockids == block[0])
			    value = block[1]
			    newdata = Numeric.where(cond,value,olddata)
			    var.data = newdata

        return vars


    def __comboFields(self, defaultval, varstor, data):
	"""
        EIB Note: This is the old default action of combining fields from
	          overlapping source results
		  verstor -> varstor[cnt][thisindex]
	"""
        try:
            tmpv = defaultval * int(varstor.mustmap)
            tmpd = Numeric.array(tmpv, varstor.type)
            tmpd = Numeric.resize(tmpd, varstor.shape)
            data = data + varstor.data - tmpd
        except:
            data  = varstor.data


    def str(self):

        print >> self.fp, 'simulation specs:'
        self.specs.str()

        print >> self.fp, 'meshes:'
        for mesh in self.mlist:
            mesh.fillMesh()
            #mesh.str()

        print >> self.fp, 'timesteps:'
        for timestep in self.timesteps.tlist:
            timestep.vlist = self.getValues(timestep.vlist)
        #self.timesteps.str()

if __name__=='__main__':

    "Test modifying an XML storage object"

    from PYTHONutils import uTestOpts

    """
    Command-line options are for two XML files, a meshfile, 
    and output file and debug.
    The XML file location is hardwired in dfltdir
    """
    dfltdir    = '../../scripts/test_TBrookParse/samples/'
    T_Sim      = 'map_output/map.TBrook.xml'
    T_Sim2     = 'map-b_output/map-b.TBrook.xml'
    meshfile   = 'mesh-ab_f.exo'
    fileout    = 'mod_TSim1_and_2'

    opts = uTestOpts('abmod',
                     defaults={'a' : dfltdir + T_Sim,  # first XML file
                               'b' : dfltdir + T_Sim2, # second XML file
                               'm' : meshfile,         # mesh file
                               'o' : fileout,          # GMV output file
                               'd' : False },          # debug flag
                     dir=dfltdir)
    (opt,args) = opts.parse_args()
    opts.header(__file__)

    fpwatch    = open('crap.dat','w') #sys.stdout
    try:
        from XMLgetStorageObject   import getStorageObject
        from GMVwriteStorageObject import GMVwriteStorageObject

        tstorages     = []
        tmeshs        = []

        "Load First simulation"
        tstorage      = getStorageObject(opt.a,fp=fpwatch,debug=0)
        itimeout      = 0
        t             = tstorage.tlist[itimeout].time
        vars          = tstorage.tlist[itimeout].vlist

        tstorages.append(tstorage)
        tmeshs.append(tstorage.mlist[0])

        "Load Second simulation"
        tstorage      = getStorageObject(opt.b,fp=fpwatch,debug=0)
        itimeout      = 1
        t             = tstorage.tlist[itimeout].time
        vars          = tstorage.tlist[itimeout].vlist

        tstorages.append(tstorage)
        tmeshs.append(tstorage.mlist[0])
        csf           = tstorage.specs.csf

        """
        Modify the tstorage object by reading in a different mesh
        from a meshfile
        """
        # copy in the meshfile
        wdir          = os.getcwd()
        cmd           = 'cp ' + dfltdir + '%s .' %(opt.m)
        os.system(cmd)

        # create the modified storage object using the meshfile
        modstorage    = None
        modstorage    = modStorageObject(tstorages,tmeshs,modstorage,opt.m,
                                         meshscale_factor=csf,
                                         stepids=[0,1],fp=fpwatch,debug=opt.d)
        modstorage.str()
        itimeout      = 0
        t             = modstorage.tlist[itimeout].time
        vars          = modstorage.tlist[itimeout].vlist
        writer        = GMVwriteStorageObject(opt.o, 'ascii',
                                              modstorage.mlist[0],
                                              t,itimeout,vars,fpwatch=fpwatch,
					      debug=opt.d)

        # remove the meshfile
        cmd = 'rm %s ' %(meshfile)
        os.system(cmd)

    except:
        print >> fpwatch, "---> Test failed in some aspect <---"
        print >> fpwatch, "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise
