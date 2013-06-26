"""

 MODMesh

 -----------------------------------------------------------------------------
  Purpose:

    Instantiates and provides methods for MODMesh objects.

  Public Interface(s):

    mMsh = MODMesh(oldmesh_type,oldmesh_msnames,meshfile,meshdir)
    mMsh.fillMesh()
    mMsh.fillCells(meshspace)
    mMsh.fillVertices(meshspace)
    mMsh.fillFaces(meshspace)
    mMsh.fillEdges(meshspace)
    mMsh.str()
   
  Contains:
    class modMesh
        __init__(self,oldmesh_type,oldmesh_msnames,meshfile,meshdir,
                 scale_factor=1.0,clean='clean',debug=0)
        metastr(self)

    From baseMesh class, it inherits
        fillMesh(self)
        fillCells(self,meshspace)
        fillVertices(self,meshspace)
        fillFaces(self,meshspace)
        fillEdges(self,meshspace)


    Unit Test Block

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
    # modify dot-dots appropriately
    parsedir  = os.path.abspath(thisdir+'../') 

    sys.path.append(parsedir)

from BASEutils          import baseMesh
from MODMeshSpace       import modMeshSpace
from EXOutils           import getEXOMesh

class modMesh(baseMesh):

    def __init__(self,oldmesh_type,oldmesh_msnames,meshfile,meshdir,
                 scale_factor=1.0,clean='clean',fp=sys.stdout,debug=0):

        "initialise the mesh object"

        baseMesh.__init__(self)

        self.type           = oldmesh_type
	
        self.__clean          = clean
        self.__debug          = debug
        self.__fp             = fp
        self.__scale_factor   = scale_factor

        """
        currently we can only do hex to hex mapping so the
        original mesh must be a hex mesh
        """

        assert self.type == 'HEX' or self.type == 'TET'
        
        """
        get the new mesh by reading in from a meshfile
        - currently we do this from EXOutils
        """
        wdir    = os.getcwd()
        newmesh = getEXOMesh(meshfile,fp=self.__fp,debug=self.__debug)
        newmesh.getData()
        os.chdir(wdir)

        self.name  = newmesh.name  

        "now start defining the attributes for this new mesh object"

        for name in oldmesh_msnames:
            # i.e meshspace names will not change
            thispace  = modMeshSpace(name,newmesh,self.__scale_factor)
            self.mslist.append(thispace)
            if thispace.unpermutation != None:
                self.updata[thispace.name] = thispace.unpermutation.data
            else:
                self.updata[thispace.name] = []

    def metastr(self):

        print >> self.__fp, self.name
        print >> self.__fp, self.type
        print >> self.__fp, self.__scale_factor

        print 'meshspace characteristics'
        for mshspc in self.mslist:
            mshspc.str()

if __name__=='__main__':
	
    def uTcleanUp(file):
        # delete the file
        import os
        cmd = 'rm %s ' %(file)
        os.system(cmd)
 
    # component test for modifying mesh

    from PYTHONutils import uTestOpts
    dfltdir    = '../../../scripts/test_TBrookParse/samples/'

    """
    Note that the directory is NOT used as part of the default filename.
    This is because the file is copied from the dfltdir prior to use.
    """
    opts       = uTestOpts('fdc', 
		           dir      = dfltdir,
                           defaults = {'f' : 'mesh-a_f.exo',
                                       'c' : False},
                           actions  = {'c' : 'store_true'}
			   )
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    fpwatch = open('crap.dat','w') #sys.stdout
    
    try:
	if opt.c:
	    print "Cleaning up from a past test run...",
	    import os
	    files = opt.f +' '+ opt.f+'.gmv' +' crap.dat'
	    cmd = 'rm %s ' % files
	    os.system(cmd)
	    print "Cleanup successful!!"

	else:
            from GMVwriteStorageObject import GMVwriteStorageObject
    
            gmvout          = opt.f + '.gmv'
            file            = dfltdir + opt.f
            cmd             = 'cp %s .' %(file)
            os.system(cmd)
            oldmesh_type    = 'HEX'
            oldmesh_msnames = ['CELL', 'VERTEX', 'FACE']
        
            wdir            = os.getcwd()
            meshdir         = wdir + '/EXOutils'
        
            try:
                modifymesh = modMesh(oldmesh_type, oldmesh_msnames,
                                 opt.f, meshdir, clean='cleanest',
                                 scale_factor=0.05, fp=fpwatch,
                                 debug=1)
            
            except:
                print >> fpwatch, \
    		'\nProblems occurred in generating the mesh from %s' %(opt.f)
                print >> fpwatch, 'Exiting this component test\n'
                uTcleanUp(opt.f)
                raise

            # populate the modified mesh
            modifymesh.fillMesh()
    
            writer = GMVwriteStorageObject(gmvout, 'ascii', modifymesh,
                                       0.0, 0, [], fpwatch = fpwatch, debug=1)
            #modifymesh.str()
    
            
    except:
        print >> fpwatch, "---> Test failed in some aspect <---"
        print >> fpwatch, "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        uTcleanUp(opt.f)
        uTcleanUp(gmvout)
        if opt.d: raise

