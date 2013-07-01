"""

 MODMeshSpace

 -----------------------------------------------------------------------------
  Purpose:

    Instantiate a modified mesh object and perform the modification.

  Public Interface(s):

   mMS = modMeshSpace(name,NewMesh,scale_factor)
   mMS.str()

  Contains:
    class modMeshSpace
        __init__(self,name,NewMesh,scale_factor)
        str(self)

    Unit Test Block

  Version:
    $ID$

  Author(s): Sharen Cummins (scummins@lanl.gov)
             Erin Iesulauro Barker (eibarker@lanl.gov)
 -----------------------------------------------------------------------------
"""

if __name__=='__main__':
    " Set sys.path for component testing mode "
    import os, sys
    thisdir   = os.path.dirname(__file__)
    # modify dot-dots appropriately
    parsedir  = os.path.abspath(thisdir+'../') 

    sys.path.append(parsedir)

from BASEutils         import baseMeshSpace
from MODMeshVariable   import modMeshVariable

class modMeshSpace(baseMeshSpace):

    def __init__(self,name,NewMesh,scale_factor):

        "initialise the mesh object"
        baseMeshSpace.__init__(self)
	
        "Modify the MeshSpace object"

        self.name           = name
        newvar              = None

        try:
           import numpy.oldnumeric as Numeric
        except ImportError:
           import Numeric
        except:
           raise

        #set variable attributes given the 'NewMesh'

        if name == 'CELL':
            self.size   = NewMesh.ncells
            newvar      = modMeshVariable(name='VERTICES',rank=2,
                                          shape=[NewMesh.ncells,NewMesh.nvc],type='i',
                                          data=NewMesh.connectivity,
                                          mesh = NewMesh.name,meshspace = name)
            self.vlist.append(newvar)
            
            newvar      = modMeshVariable(name='CENTROIDS',rank=2,
                                          shape=[NewMesh.ncells,NewMesh.ndim],type='d',
                                          data=NewMesh.cellcentroids,
                                          mesh = NewMesh.name,meshspace = name)
            newvar.data = scale_factor*newvar.data
            self.vlist.append(newvar)

            newvar      = modMeshVariable(name='BLOCKID',rank=1,
                                          shape=[NewMesh.ncells],type='i',
                                          data=NewMesh.cellblocks,
                                          mesh = NewMesh.name,meshspace = name)
            self.vlist.append(newvar)

            newvar      = modMeshVariable(name='NUMNEIGHBORS',rank=1,
                                          shape=[NewMesh.ncells],type='i',
                                          data=NewMesh.cellnumneighbors,
                                          mesh = NewMesh.name,meshspace = name)
            self.vlist.append(newvar)

            newvar      = modMeshVariable(name='NVC',rank=1,
                                          shape=[1],type='i',
                                          data=NewMesh.nvc,
                                          mesh = NewMesh.name,meshspace = name)
            self.vlist.append(newvar)
            
            self.unpermutation = modMeshVariable(name='UNPERMUTE',rank=1,
                                                 shape=[NewMesh.ncells],type='i',
                                                 data=None,
                                                 mesh = NewMesh.name,meshspace = name)

        if name == 'VERTEX':
            self.size     = NewMesh.nnodes
            newvar        = modMeshVariable(name='COORDS',rank=2,
                                            shape=[NewMesh.nnodes,NewMesh.ndim],
                                            type='d', data=NewMesh.coords,
                                            mesh = NewMesh.name,meshspace = name)
            newvar.data = scale_factor*newvar.data
            self.vlist.append(newvar)
            self.unpermutation = modMeshVariable(name='UNPERMUTE',rank=1,
                                                 shape=[NewMesh.nnodes],type='i',
                                                 data=None,
                                                 mesh = NewMesh.name,meshspace = name)

        if name == 'FACE':
            self.size = 0
            self.unpermutation = modMeshVariable(name='UNPERMUTE',rank=1,
                                                 shape=[NewMesh.nnodes],type='i',
                                                 data=None,
                                                 mesh = NewMesh.name,meshspace = name)
            
        if name != 'FACE' and name != 'NONE':
            #set permute mapping array to identity array 
            self.unpermutation.data = Numeric.array([0],'i')
            self.unpermutation.data = Numeric.resize(self.size,[self.size])
            for i in range(self.size):
                self.unpermutation.data[i] = i+1

    
if __name__=='__main__':


    from PYTHONutils import uTestOpts

    dfltdir = '../../../scripts/test_TBrookParse/samples/'
    """
    Note that the directory is NOT used as part of the default filename.
    This is because the file is copied from the dfltdir prior to use.
    """
    opts = uTestOpts('fsdc', # meshfile, scale_factor, debug
                     defaults = {'f' : 'good-mesh-1.exo',
                                 's' : 1.0,
			         'c' : False},
                     types    = {'s' : 'float'},
		     actions  = {'c' : 'store_true'},
                     dir      = dfltdir)
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    try:

	if opt.c:
	    print "\nCleaning up from a past rest run..."
            import os
            cmd = 'rm %s ' %(opt.f)
            os.system(cmd)
	    print "Cleanup successful!!\n"

        else:
            from EXOutils import getEXOMesh
        
            wdir     = os.getcwd()
            wdir     = wdir + '/EXOutils'
        
            cmd      = 'cp ' + dfltdir + '%s .' %(opt.f)
            os.system(cmd)
        
            #newmesh  = getEXOMesh(opt.f,wdir,clean='clean',debug=1)
            newmesh  = getEXOMesh(opt.f,debug=1)
            newmesh.getData()
        
            space    = modMeshSpace('CELL',newmesh,opt.s)
            space.metastr()
            space    = modMeshSpace('VERTEX',newmesh,opt.s)
            space.metastr()
            space    = modMeshSpace('FACE',newmesh,opt.s)
            space.metastr()
        
	    print '\n--> Note: to clean up, re-execute same command line plus "-c" \n'
        
    except:
        print "---> Test failed in some aspect <---"
        print "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise

