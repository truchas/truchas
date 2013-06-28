"""

 MAPVariable

 -----------------------------------------------------------------------------
  Purpose:

     Provides functionality to map Truchas variables from one mesh to another.
     
  Public Interface(s):

    v = mapVariable(oldvar,oldmesh,newmesh,wdir)
    v.getData() 
    v.str()     

  Contains:
    class mapVariable
        __init__(self,oldvar,oldmesh,newmesh,wdir,
                 mustmap=0,defval=0.0,maprule='constants_preserving',
                 clean='clean',debug=0):
        getData(self)
        metastr(self)
        str(self)

    Unit Test Block
     type 'python MAPVariable.py -H' for guidance on testing

  Version:
    $ID$

  Author(s): Sharen Cummins (scummins@lanl.gov)
             Larry Cox      (ljcox@lanl.gov) (unit testing)
 -----------------------------------------------------------------------------
"""
import sys, os
thisdir = os.path.abspath(os.path.dirname(__file__))

if __name__=='__main__':
    print "\nFor component test in %s \n" %(__file__)
    parserdir   = thisdir + '/../../../TBrookParser'
    sys.path.append(parserdir)

class mapVariable:

    def __init__(self,oldvar,oldmesh,newmesh,
                 mustmap = 0,
                 defval  = 0.0,
                 maprule = 'constants_preserving',
                 fp      = sys.stdout,
                 debug   = 0):

        self.name      = oldvar.name
        self.nickname  = oldvar.nickname
        self.rank      = oldvar.rank
        self.shape     = []
        for i in range(len(oldvar.shape)):
            self.shape.append(oldvar.shape[i])
        self.type      = oldvar.type
        self.oldata    = oldvar.data
        self.meshspace = oldvar.meshspace
        self.mesh      = newmesh.name
        self.data      = None

        self.oldmesh   = oldmesh
        self.newmesh   = newmesh
        self.mustmap   = mustmap
        self.maprule   = maprule
        self.oldintgrl    = 0.0    #diagnostics to check conservative mapping
        self.newintgrl    = 0.0
        self.roundoff     = 0.001  # provide round off required in converting to integer format
        self.deftargetval = 0.0    # default value for all variables that cannot (or will not)
                                   #  be mapped onto the new mesh
        self.default      = defval # default value for target variable in regions not contained
                                   #  within the source mesh domain
        self.maxwarn      = 1000   # number of warning messages for searches performed by 
                                   #  mapper for cells outside of the domain of the source mesh 
        self.debug        = debug
        self.fp           = fp

        #input assertions for this class
        
        for meshspace in oldmesh.mslist:
            if meshspace.name == self.meshspace:
                if oldvar.name !='RHS':
                    try:
                        assert len(self.oldata) == meshspace.size
                    except:
                       print >> self.fp, 'Variable %s size (%i) does not equal the %s meshspace size (%i)' \
                             %(oldvar.name,len(self.oldata),self.meshspace,meshspace.size)
                       print >> self.fp, 'We cannot modify this variable..please consult postprocessor',\
                             'support with this error'
                       sys.exit(1)
                if oldvar.name == 'RHS':
                    #special care for the solid mechanics RHS variable
                    try:
                        assert len(self.oldata) == 3*meshspace.size
                    except:
                        print >> self.fp, 'RHS variable size (%i) does not' %(len(self.oldata)) ,\
                              'correspond to ndim x %s meshspace size (%i)' \
                              %(self.meshspace,meshspace.size)
                        print >> self.fp, 'We cannot modify this variable..please consult postprocessor',\
                              'support with this error'
                        sys.exit(1)

        try:
            assert self.meshspace != 'EDGE'
        except:
            print >> self.fp, 'Currently we cannot modify variables that live on the %s meshspace' \
                  %(self.meshspace)
            print >> self.fp, 'This feature will be made available in future releases'
            sys.exit(1)
    
    def getData(self):

        try:
            import numpy.oldnumeric as Numeric
        except ImportError:
            import Numeric
        except:
            raise

        try:
            import gridmap
        except ImportError:
            print >> self.fp, 'error importing gridmap'
            raise

        global old_ncells
        global old_nvc
        global old_connectivity
        global old_cblocks
        global old_nnodes
        global old_coords
        global old_nfaces

        global new_ncells
        global new_nvc
        global new_connectivity
        global new_cblocks
        global new_nnodes
        global new_coords
        global new_nfaces

        try:    
            """
            old mesh values already exist
            so just use the global old_* variables
            """
            X                          = old_ncells
            assert len(self.old_data) == old_ncells
        except:
            "get old mesh values"
            for meshspace in self.oldmesh.mslist:
                if meshspace.name == 'CELL':
                    self.oldmesh.fillCells(meshspace)
                    old_ncells            = self.oldmesh.cells['N']
                    old_connectivity      = Numeric.array([1],'i')
                    old_connectivity      = Numeric.resize(old_connectivity,
                                                           self.oldmesh.cells['VERTICES'].shape)
                    old_connectivity[:,:] = self.oldmesh.cells['VERTICES'][:,:]
                    old_connectivity      = Numeric.transpose(old_connectivity)                    
                    if self.oldmesh.cells.has_key('NVC'):
                        old_nvc      = self.oldmesh.cells['NVC']
                    else:
                        old_nvc      = len(old_connectivity[:,0])
                    if self.oldmesh.cells.has_key('BLOCKID'):
                        old_cblocks  = self.oldmesh.cells['BLOCKID']
                    else:
                        #for now create dummy array for the mesh cell blocks
                        old_cblocks  = Numeric.array([1],'i')
                        old_cblocks  = Numeric.resize(old_cblocks,[old_ncells])
                if meshspace.name == 'VERTEX':
                    self.oldmesh.fillVertices(meshspace)
                    old_nnodes       = self.oldmesh.vertices['N']
                    old_coords       = Numeric.array([0.0],'d') #type(self.oldmesh.vertices['COORDS'][0,0])) #'d') 
                    old_coords       = Numeric.resize(old_coords,self.oldmesh.vertices['COORDS'].shape)
                    old_coords[:,:]  = self.oldmesh.vertices['COORDS'][:,:]
                    old_coords       = Numeric.transpose(old_coords)

        try:
            """
            new mesh values already exist
            so just use the global new_* variables 
            """
            X                          = new_ncells
            assert len(self.new_data) == new_ncells
        except:
            "get new mesh values"
            for meshspace in self.newmesh.mslist:
                if meshspace.name == 'CELL':
                    self.newmesh.fillCells(meshspace)
                    new_ncells            = self.newmesh.cells['N']
                    new_connectivity      = Numeric.array([1],'i')
                    new_connectivity      = Numeric.resize(new_connectivity,
                                                           self.newmesh.cells['VERTICES'].shape)
                    new_connectivity[:,:] = self.newmesh.cells['VERTICES'][:,:]
                    new_connectivity      = Numeric.transpose(new_connectivity)
                    if self.newmesh.cells.has_key('NVC'):
                        new_nvc      = self.newmesh.cells['NVC']
                    else:
                        new_nvc      = len(new_connectivity[:,0])
                    if self.newmesh.cells.has_key('BLOCKID'):
                        new_cblocks  = self.newmesh.cells['BLOCKID']
                    else:
                        #for now create dummy array for the mesh cell blocks
                        new_cblocks  = Numeric.array([1],'i')
                        new_cblocks  = Numeric.resize(new_cblocks,[new_ncells])
                if meshspace.name == 'VERTEX':
                    self.newmesh.fillVertices(meshspace)
                    new_nnodes       = self.newmesh.vertices['N']
                    new_coords       = Numeric.array([0.0],'d') #type(self.oldmesh.vertices['COORDS'][0,0])) #'d') 
                    new_coords       = Numeric.resize(new_coords,self.newmesh.vertices['COORDS'].shape)
                    new_coords[:,:]  = self.newmesh.vertices['COORDS'][:,:]
                    new_coords       = Numeric.transpose(new_coords) 
                if meshspace.name == 'FACE':
                    self.newmesh.fillFaces(meshspace)
                    new_nfaces       = self.newmesh.faces['N']

        "alter shape of variable depending on the meshspace the variable lives on"
        if self.meshspace == 'CELL':
            self.shape[0] = new_ncells
        if self.meshspace == 'VERTEX':
            self.shape[0] = new_nnodes
        if self.meshspace == 'FACE':
            self.shape[0] = new_nfaces

        ndim       = len(new_coords[:,0])

        #special treatment for solid mechanics RHS variable

        if self.name == 'RHS':
            self.shape[0] = ndim*new_nnodes

        """
        must resize all variables to the new mesh size regardless
        if they are to be mapped or not
        """
        self.data  = Numeric.array([self.deftargetval],self.type)
        self.data  = Numeric.resize(self.data,self.shape)

        assert self.newmesh.type == 'HEX' or self.newmesh.type == 'TET', \
               'New mesh is not a HEX or TET mesh so will not map variable %s'\
               %(self.name)

        if self.mustmap:

            "map restart cell fields using the hex to hex mapping capability"

            data       = None
            data       = Numeric.array([self.default],self.type)
            data       = Numeric.resize(data,[self.shape[0]])

            print >> self.fp
            print >> self.fp, 'Mapping %s variable %s'%(self.meshspace,self.name),\
                  'using hex-to-hex mapping with map rule %s'%(self.maprule) 
            print >> self.fp

            ierr  = 0

            if self.rank == 1:

                self.oldintgrl, self.newintgrl, ierr = gridmap.getfield(
                    self.maprule,
                    self.maxwarn,
                    self.default,
                    ndim,
                    old_ncells,
                    old_nvc,
                    old_nnodes,
                    self.oldata,
                    old_connectivity,
                    old_coords,
                    old_cblocks,
                    new_ncells,
                    new_nvc,
                    new_nnodes,
                    data,
                    new_connectivity,
                    new_coords,
                    new_cblocks)

                if str(self.type) == 'i':
                    data      = (1.0+self.roundoff)*data
                
                self.data = data.astype(str(self.type))

                if self.oldintgrl > 0 and self.newintgrl < self.roundoff*self.oldintgrl:
                    print >> self.fp
                    print >> self.fp, 'Problems mapping data for variable %s ' %(self.name)
                    print >> self.fp, 'Please check each mesh for internal consistencies'
                    print >> self.fp, 'and for consistencies with each other'
                    return
        
            if self.rank == 2:

                for j in range(self.shape[1]):

                    self.oldintgrl, self.newintgrl, ierr = gridmap.getfield(
                        self.maprule,
                        self.maxwarn,
                        self.default,
                        ndim,
                        old_ncells,
                        old_nvc,
                        old_nnodes,
                        self.oldata,
                        old_connectivity,
                        old_coords,
                        old_cblocks,
                        new_ncells,
                        new_nvc,
                        new_nnodes,
                        data,
                        new_connectivity,
                        new_coords,
                        new_cblocks)
                    
                    if str(self.type) == 'i':
                        data       = (1.0+self.roundoff)*data

                    self.data[:,j] = data.astype(str(self.type))

            if self.rank == 3:

                for j in range(self.shape[1]):

                    for k in range(self.shape[2]):

                        self.oldintgrl, self.newintgrl, ierr = gridmap.getfield(
                            self.maprule,
                            self.maxwarn,
                            self.default,
                            ndim,
                            old_ncells,
                            old_nvc,
                            old_nnodes,
                            self.oldata,
                            old_connectivity,
                            old_coords,
                            old_cblocks,
                            new_ncells,
                            new_nvc,
                            new_nnodes,
                            data,
                            new_connectivity,
                            new_coords,
                            new_cblocks)
                        
                        if str(self.type) == 'i':
                            data         = (1.0+self.roundoff)*data
                    
                        self.data[:,j,k] = data.astype(str(self.type))

            if ierr < 0:
                print >> self.fp
                print >> self.fp, 'Meshes are inconsistent so cannot map %s variable ' %(self.name)
                print >> self.fp, 'Please check the geometry of both meshes before \n',\
                      ' proceeding further with this map'
                print >> self.fp, 'In particular check that the coordinate scale factor \n',\
                      ' used in the original Truchas simulation is consistent \n',\
                      ' with what you chose in the postprocessor.'
                print >> self.fp
                sys.exit(1)

        else:

            #fill in all other variables with default values

            print >> self.fp
            print >> self.fp, 'We do not map %s variable %s, so this array set to default values ' \
                  %(self.meshspace,self.name) 
            print >> self.fp
            

    def metastr(self):
	print >> self.fp, 'name:      ',self.name
	print >> self.fp, 'rank:      ',self.rank
	print >> self.fp, 'shape:     ',self.shape
	print >> self.fp, 'type:      ',self.type,
	print >> self.fp, 'meshspace: ',self.meshspace
	print >> self.fp, 'maprule:   ',self.maprule

    def str(self):
	self.metastr()
	print >> self.fp, 'Data:',
        print >> self.fp, self.data
	print >> self.fp, 'old integral: ',self.oldintgrl,' vs. new integral: ',self.newintgrl


if __name__== '__main__':
    
    docstr = \
    """
    Built-in unit tests for MAPVariable.py
    Usage: python MAPvariable.py [-H] [-d] [-i] [-f <params>]
           -H          : print this message and exit
           -d          : set the debug flag to true
           -i          : write out ASCCI intermediate mesh files
           -f <params> : import test parameters from file <params>.py

    Default: use file params_MAPVar.py, debug=False

    See the file README_extex for more information.

    Test parameters are imported from the file given with the -f
    command line option. That file is a python-syntax assignment of
    various parameters.  Provide only the prefix, the suffix must be
    '.py'.  The parameters required are:

      InputsDir   : the source directory for input files
      TBprefix    : the prefix to a a TBrook.xml file
                    (expands to <TBprefix>_output/<TBprefix>.TBrook.xml)
      mesha_file  : a dict of 'a' input mesh file and type
                    an empty DICT means use the TBrook mesh for 'a'
      meshb_file  : a dict of 'b' input mesh file and type
                    an empty DICT means use the TBrook mesh for 'b'
      CSF         : coordinate scale factor for mesh 'a'
      mapRule     : the mapping rule: 'conservative' or 'constants_preserving'
    """

    #-------------------------
    # An output function for this test block
    def writeOutMesh(tstorage,mesh,outfile='tmesh.dat'):
        "write out a truchas mesh in ascii forma"

        #get all data values associated with this mesh
        for meshspace in mesh.mslist:
            tstorage.getValues(meshspace.vlist)
            #place them in a convenient storage 
            if meshspace.name == 'CELL':
                mesh.fillCells(meshspace)
            elif meshspace.name == 'VERTEX':
                mesh.fillVertices(meshspace)
            
        file     = open(outfile,'w')
        fmt      = '%12.5E '
        ifmt     = '%10i '
        nperline = 1
        data     = mesh.cells['VERTICES']
        ncells   = mesh.cells['N']
        cblcks   = mesh.cells['BLOCKID']
        nvc      = 8
        nnodes   = mesh.vertices['N']
        t        = ifmt%ncells + ifmt%nnodes + ifmt%nvc + '\n'
        file.write(t)

        s = ''
        for slice in data:
            for j in slice:
                s = s + ifmt%j
            s = s + '\n'
        file.write(s)

        data   = mesh.cells['BLOCKID']
        s      = ''
        icount = 0
        for j in data:
            s = s + ifmt%j
            if (icount%nperline == 0 and icount < len(data)-1):
                s = s + '\n'
        file.write(s)
    
        data = mesh.vertices['COORDS']
        s = ''
        for slice in data:
            for j in slice:
                s = s + fmt%j
            s = s + '\n'
        file.write(s)
        file.close()
    #-----------------------------
    # Begin Unit test Block
    
    from PYTHONutils import uTestOpts
    "Parse command-line for input parameters file and options."
    o = uTestOpts('fdHi',
                  defaults={'f' : 'param_MAPVar', 'H': False, 'i' : False},
                  actions ={'H' : 'store_true', 'i' : 'store_true'},
                  dir     = thisdir)
    (opt,a) = o.parse_args()
    if (opt.H):
        print docstr
        sys.exit(0)
        
    fpwatch = sys.stdout  
    "Import the parameters file"
    try:
        params = __import__(opt.f) # a bunch of parameters    
    except:
        print >> fpwatch, "\nFailed to import paramaters from file '%s.py'\n" %(opt.f)
        raise

    """
    Step one: Load an XML field variable
    """
    TBfile = params.InputsDir + params.TBfile
    from XMLgetStorageObject import getStorageObject
    tstorage = getStorageObject(TBfile,fp=fpwatch)
    csf      = tstorage.specs.csf
    "Load a cell variable on this truchas mesh"
    for var in tstorage.tlist[0].vlist:
        if var.rank == 1 and var.meshspace == 'CELL' :
            thevar = var
    thevars = []
    thevars.append(thevar)
    tstorage.getValues(thevars)
    thevar.str()
    
    "List of desired meshspace names"
    msnames = ['CELL', 'VERTEX', 'FACE']
    "Establish the directory paths needed."
    wdir          = os.getcwd()
    meshdir       = wdir + '/../EXOutils'
    
    """
    Step Two:
    Load two additional mesh files.
    Mesh 'a' and mesh 'b' can be Exodus or Truchas meshes.  The
    selection is triggered by the presence (or not) of
    params.mesha/meshb dictionary keys.
    """
    from MODutils            import modMesh

    "---Mesh 'a'---"
    if 'mesha' in dir(params) and params.mesha.keys():
        "get mesh 'a' from Exodus file..."
        try:
            meshfile,meshtype = params.mesha.popitem()
            cmd   = 'cp %s%s .' %(params.InputsDir,meshfile)
            os.system(cmd)
            mesha = modMesh(meshtype,msnames,meshfile,meshdir,
                            scale_factor = params.CSF,
                            clean        = 'cleanest',
                            fp           = fpwatch,
                            debug        = opt.d)
            ofile = 'exomesha.ascii'
            cmd   = 'rm %s' %(meshfile)
            os.system(cmd)
        except:
            print >> fpwatch, '\nProblems occurred in generating the',\
                  'Exodus mesh from %s of type %s' %(meshfile,meshtype)
            print >> fpwatch, 'Exiting this component test \n'
            raise

    else:
        "get mesh 'a' from Truchas output"
        mesha = tstorage.mlist[0]
        for meshspace in mesha.mslist:
            tstorage.getValues(meshspace.vlist)
        ofile = 'tmesha.ascii'

    "write out Truchas mesh 'a' to an ascii format for component testing"
    if opt.i: writeOutMesh(tstorage,mesha,outfile=ofile)
    
    "---Mesh 'b'---"
    if 'meshb' in dir(params) and params.meshb.keys():
        "get mesh 'b' from exodus file..."
        try:
            meshfile,meshtype = params.meshb.popitem()
            cmd   = 'cp %s%s .' %(params.InputsDir,meshfile)
            os.system(cmd)
            meshb = modMesh(meshtype,msnames,meshfile,meshdir,
                            scale_factor = tstorage.specs.csf,
                            clean        = 'cleanest',
                            fp           = fpwatch,
                            debug        = opt.d)
            ofile = 'exomeshb.ascii'
            cmd   = 'rm %s' %(meshfile)
            os.system(cmd)
        except:
            print >> fpwatch, '\nProblems occurred in generating the',\
                  'Exodus mesh from %s of type %s' %(meshfile,meshtype)
            print >> fpwatch, 'Exiting this component test \n'
            raise
    else:
        "get mesh 'b' from Truchas output"
        meshb           = tstorage.mlist[0]
        for meshspace in meshb.mslist:
            tstorage.getValues(meshspace.vlist)    
        ofile = 'tmeshb.ascii'

    "write out Truchas mesh 'b' to an ascii format for component testing"
    if opt.i: writeOutMesh(tstorage,meshb,outfile=ofile)

    "Done loading meshes 'a' and 'b'"

    
    """
    Step 3:
    Change the shape of the XML field variable so it fits on mesh 'a'.
    We do this to give mesh 'a' a field variable to map to mesh 'b'
    """
    for meshspace in mesha.mslist:
        if meshspace.name == thevar.meshspace:
            if thevar.name == 'RHS':
                thevar.shape[0] = 3*meshspace.size
            else:
                thevar.shape[0] = meshspace.size

    """
    Step 4 - THE TEST:
    Perform the mapping from 'a' to 'b'
    """
    try:
    #if(1):
        var = mapVariable(thevar,mesha,meshb,
                          mustmap = 1,
                          maprule = params.mapRule,
                          fp      = fpwatch, 
                          debug   = opt.d)
        var.getData()
        var.str()
    except:
    #else:
        print >> fpwatch, '\nTruchas remap failed.',\
              'Please check your source and target meshes for '
        print >> fpwatch, 'internal consistencies and consistency with each other.\n'
        raise
