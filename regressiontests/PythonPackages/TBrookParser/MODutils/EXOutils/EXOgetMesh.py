import os, sys
thisdir = os.path.abspath(os.path.dirname(__file__))

if __name__=='__main__':
    print "\nFor component test in %s \n" %(__file__)
    parserdir   = thisdir + '/../../../TBrookParser'
    sys.path.append(parserdir)

class getEXOMesh:

    def __init__(self,filename,fp=sys.stdout,debug=0):

        self.filename         = filename
        self.debug            = debug
        self.fp               = fp 

        self.name             = 'ExoMesh'
        self.ncells           = 0
        self.nnodes           = 0
        self.ndim             = 3
        self.nvc              = 8
        self.nfc              = 6
        self.connectivity     = None
        self.coords           = None
        self.cellcentroids    = None
        self.cellnumneighbors = None
        self.cellblocks       = None

        
    def getData(self):

        import os, sys, string
        try:
            import numpy.oldnumeric as Numeric
        except ImportError:
            import Numeric
        except:
            raise

        try:
            import exodus
        except ImportError:
            print >> self.fp, 'error importing exodus.so'
            raise

        self.nnodes, self.ncells = exodus.GetSizes(self.filename)

        if self.debug:
            print >> self.fp
            print >> self.fp, 'stats in EXOgetMesh.py'
	    self.metastr()
            print >> self.fp
        
        self.coords         = Numeric.array([0],'d')
        self.coords         = Numeric.resize(self.coords,[self.nnodes,self.ndim])
        self.connectivity   = Numeric.array([0],'i')
        self.connectivity   = Numeric.resize(self.connectivity,[self.ncells,self.nvc])
        self.cellblocks     = Numeric.array([0],'i')
        self.cellblocks     = Numeric.resize(self.cellblocks,[self.ncells])

        exodus.GetMesh(self.filename,
                       self.connectivity,
                       self.coords,
                       self.cellblocks)

        #catch any exceptions associated with getting the mesh
        AllZero = 'Coords all zero'

        def __nozeroes(coords):
            T = Numeric.sometrue(Numeric.greater(coords,0)
                                 | Numeric.less(coords,0))
            "Accomodate Numpy/Numeric differences"
            try:
                if not T.any(): raise AllZero
            except:
                if not Numeric.greater(T,0): raise AllZero

        try:
            __nozeroes(self.coords)
        except AllZero:
            print >> self.fp 
            print >> self.fp, 'We will not continue extracting exodus mesh from %s' %(self.filename)
            print >> self.fp, 'as the mesh coordinates were not obtained.'
            print >> self.fp 
            sys.exit(1)
        
        #additional geometry fields
        self.cellcentroids    = Numeric.array([0],'d')
        self.cellcentroids    = Numeric.resize(self.cellcentroids,[self.ncells,self.ndim])
        self.cellnumneighbors = Numeric.array([0],'i')
        self.cellnumneighbors = Numeric.resize(self.cellnumneighbors,[self.ncells])

        exodus.GetCellFields(self.ncells,
                             self.nnodes,
                             self.ndim,
                             self.nvc,
                             self.connectivity,
                             self.coords,
                             self.cellcentroids,
                             self.cellnumneighbors)

    def metastr(self):
	print >> self.fp, 'number of cells:    ',self.ncells
	print >> self.fp, 'number of nodes:    ',self.nnodes
	print >> self.fp, 'number of vertices: ',self.nvc
	print >> self.fp, 'number of faces:    ',self.nfc

    def str(self):
	self.metastr()
        print >> self.fp, 'coords'
        print >> self.fp, self.coords
        print >> self.fp, 'connectivity'
        print >> self.fp, self.connectivity
        print >> self.fp, 'cell centroids'
        print >> self.fp, self.cellcentroids
        print >> self.fp, 'cell num neighbors'
        print >> self.fp, self.cellnumneighbors
        print >> self.fp, 'cell block ids'
        print >> self.fp, self.cellblocks


if __name__== '__main__':
    
    from PYTHONutils import uTestOpts

    dfltdir = '../../../../scripts/test_TBrookParse/samples/'
    opts = uTestOpts('fd', dir=dfltdir, defaults={'f': 'mesh-a_f.exo'})
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    fpwatch = sys.stdout
    try:
        wdir     = os.getcwd()
        cmd      = 'cp %s%s .' %(dfltdir,opt.f)
        os.system(cmd)
        try:
            newmesh  = getEXOMesh(opt.f,fp=fpwatch,debug=opt.d)
            newmesh.getData()
            
        except:
            print >> fpwatch 
            print >> fpwatch, 'Problems occurred in generating the Exodus mesh from %s' %(opt.f)
            print >> fpwatch, 'Exiting this component test'
            print >> fpwatch 
            raise

        newmesh.str()

    except:
        print >> fpwatch, "---> Test failed in some aspect <---"
        print >> fpwatch, "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise


