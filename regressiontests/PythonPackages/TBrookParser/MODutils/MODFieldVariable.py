

"""

 MODFieldVariable

 -----------------------------------------------------------------------------
  Purpose:

     Instantiates and remaps a variable from one source (mesh) to a
     different mesh.
     The result is a new variable object mapped onto a different mesh.

  Public Interface(s):

    mFV = modFieldVariables(oldvar,oldmesh,newmesh,mapvardir)
    mFV.str()

  Contains:
    class modFieldVariable
        __init__(self,oldvar,oldmesh,newmesh,mapvardir,defval=0.0,debug=0)
        __setMapRules(self,file)
        str(self)

    Unit Test Block 
        uTcleanUp(file)

  Version:
    $ID$

  Author(s): Sharen Cummins (scummins@lanl.gov)
             Erin Iesulauro Barker (eibarker@lanl.gov)
 -----------------------------------------------------------------------------
"""
import os, sys

if __name__=='__main__':
    " Set sys.path for component testing mode "
    thisdir   = os.path.dirname(__file__)
    parsedir  = os.path.abspath(thisdir+'../') 

    sys.path.append(parsedir)

from PYTHONutils import getdir
from MAPutils    import mapVariable
from BASEutils   import baseVariable

class modFieldVariable(baseVariable):

    def __init__(self,oldvar,oldmesh,newmesh,mapvardir,defval=0.0,fp=sys.stdout,debug=0):

        "initialise the variable object"
        baseVariable.__init__(self)

        "Modify the FieldVariable object"

        self.name          = oldvar.name
        self.nickname      = oldvar.nickname
        self.rank          = oldvar.rank
        self.shape         = None
        self.data          = None
        self.type          = oldvar.type
        self.mesh          = oldvar.mesh      # mesh name     - not the object
        self.meshspace     = oldvar.meshspace # meshspace name - not the object

        self.defval        = defval    # target variable value in areas not
                                       #   contained by the source mesh domain

        self.debug         = debug
        self.fp            = fp
        self.oldata        = oldvar.data
        self.oldmesh       = oldmesh
        self.newmesh       = newmesh

        # create and initialize the maprule and mustmap attributes
        defaultsmaps  = getdir(__file__) + '/MAPutils/defaultmaps.txt'
        self.__setMapRules(defaultsmaps)

        try:
            var = mapVariable(oldvar,oldmesh,newmesh,
                              self.mustmap,
                              self.defval,
                              self.maprule,
                              fp=self.fp,
                              debug=self.debug)
            var.getData()
        except:
            print >> self.fp
            print >> self.fp, 'Truchas remap failed for variable %s.'  %(oldvar.name)
            print >> self.fp, 'Please check your source and target meshes for',\
                  'internal consistencies and consistency with each other.'
            print >> self.fp
            sys.exit(1)

        self.data          = var.data
        self.shape         = var.shape
        self.oldintgrl     = var.oldintgrl
        self.newintgrl     = var.newintgrl
        
    def __setMapRules(self,file):

        """
        set variable mapping attributes from file-based rules
        """
        # defaults
        self.maprule       = 'constants_preserving'
        self.mustmap       = 0  # flag to indicate this variable will
                                # be mapped, rather than ignored
        file  = open(file,'r')
        outpt = file.readlines()
        file.close()

        for line in outpt:
            if (   ('PC' in self.name and 'PHASECONC' in line)
                or ('VOF0' in self.name and 'VOF' in line)
                or (self.name in line)):
                self.mustmap    = 1
                if 'conserv' in line:
                    self.maprule = 'conservative'

    """
    def metastr(self):

        print >> self.fp, self.name, self.nickname, self.rank, \
              self.shape, self.type, self.mesh, self.meshspace

    def str(self):

	self.metastr()
        print >> self.fp, self.data
    """

if __name__=='__main__':
    
    def uTcleanUp(file):
        # delete the file
        import os
        cmd = 'rm %s ' %(file)
        os.system(cmd)

    from PYTHONutils import uTestOpts
    dfltdir = '../../../scripts/test_TBrookParse/samples/'
    dflts   = {'f' : dfltdir + 'map_output/map.TBrook.xml',
               'a' : 'Exodus',
               'b' : 'Exodus',
               'm' : 'good-mesh-1.exo'}
    
    opts = uTestOpts('fdabm', dir=dfltdir, defaults=dflts)
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    fpwatch = sys.stdout
    try:
        from XMLgetStorageObject import getStorageObject
        from MODMesh             import modMesh

        # Create the storage object from the provided filename
        tstorage = getStorageObject(opt.f,fp=fpwatch)

        # Obtain a variable on this truchas mesh
        thevars = []
        for var in tstorage.tlist[0].vlist:
            if var.rank > 0:
                thevars.append(var)
    
        tstorage.getValues(thevars)

        wdir            = os.getcwd()
        meshdir         = wdir + '/EXOutils'
        oldmesh_type    = 'HEX'
        oldmesh_msnames = ['CELL', 'VERTEX', 'FACE']
        
        # get mesh 'a' either from Truchas or Exodus
        if opt.a == 'Exodus': 
            # get mesh 'a' from Exodus file
            meshfile = opt.m   
            cmd      = 'cp ' + dfltdir + '%s .' %(meshfile)
            os.system(cmd)

            try:
                mesha  = modMesh(oldmesh_type,oldmesh_msnames,
                                 meshfile,meshdir)

            except:
                print >> fpwatch, '\nProblems occurred in generating the Exodus',\
                      'mesh from %s' %(meshfile)
                print >> fpwatch, 'Exiting this component test\n'

                uTcleanUp(meshfile)
                sys.exit(1)

        else:
            #get mesh 'a' from Truchas output        
            mesha = tstorage.mlist[0]
            for meshspace in mesha.mslist:
                tstorage.getValues(meshspace.vlist)

        mesha.str()            

        # get mesh 'b'..either from an Exodus file or a Truchas file
        if opt.b == 'Exodus':
            # get mesh 'b' from Exodus file
            meshfile = opt.m
            cmd      = 'cp ' + dfltdir + '%s .' %(meshfile)
            os.system(cmd)
            try:
                meshb = modMesh(oldmesh_type,oldmesh_msnames,meshfile,meshdir)
            except: 
                print >> fpwatch
                print >> fpwatch, 'Problems occurred in generating the Exodus mesh'\
                      ' from %s' %(meshfile)
                print >> fpwatch, 'Exiting this component test'
                print >> fpwatch
                uTcleanUp(meshfile)
                sys.exit(1)
        else:
            meshb = tstorage.mlist[0]
            #now obtain the data on the mesh 'b'
            for meshspace in meshb.mslist:
                tstorage.getValues(meshspace.vlist)
    
        try:
           import numpy.oldnumeric as Numeric
        except ImportError:
           import Numeric
        except:
           raise

        """
        Map the variables
        - First change the shape of the truchas XML variables
          so they fit on mesh 'a'
        - Then map them to mesh 'b'
        """
        wdir            = os.getcwd() + '/MAPutils'
        for meshspace in mesha.mslist:
            if meshspace.name == 'CELL':
                mesha.fillCells(meshspace)
            elif meshspace.name == 'VERTEX':
                mesha.fillVertices(meshspace)

        for thevar in thevars:
            if thevar.name == 'RHS':
                thevar.shape[0] = 3*mesha.vertices['N']
            elif thevar.meshspace == 'CELL':
                thevar.shape[0] = mesha.cells['N']
            elif thevar.meshspace == 'VERTEX':
                thevar.shape[0] = mesha.vertices['N']
            
            thevar.data = Numeric.resize(thevar.data,thevar.shape)
            thevar      = modFieldVariable(thevar,mesha,meshb,wdir,fp=fpwatch,debug=0)
            thevar.str()

    except: # uTest overall exception
        print >> fpwatch, "---> Test failed in some aspect <---"
        print >> fpwatch, "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise

