import os, sys
thisdir = os.path.abspath(os.path.dirname(__file__))

if __name__=='__main__':
    print "\nFor component test in %s \n" %(__file__)
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

from PYTHONutils  import unique

class getBlockIndicies:

    def __init__(self,cblocks,id,nnodesmesh,meshconnectivity,fp=sys.stdout,debug=0):

        self.ID               = id
        self.cellblocks       = cblocks
        self.meshconnectivity = meshconnectivity
        self.nnodesmesh       = nnodesmesh
        self.cellids          = None
        self.vtxids           = None
        self.blckconnectivity = None
        self.fp               = fp
        self.debug            = debug
    
    def getData(self):

        try:
           import numpy.oldnumeric as Numeric
        except ImportError:
           import Numeric
        except:
           raise

        try:
            import ensight
        except ImportError:
            print >> self.fp, 'error importing ensight.so'
            raise
        
        nvc = len(self.meshconnectivity[0])

        #now get global cell ids and vertex ids for this block
        self.cellids = []
        self.vtxids  = []
        for j in range(len(self.cellblocks)):
            if self.cellblocks[j] == self.ID:
                self.cellids.append(j+1)
        for j in self.cellids:
            for jj in self.meshconnectivity[j-1]:
                self.vtxids.append(jj)

        self.vtxids = unique(self.vtxids)
        self.vtxids.sort()

        ncells                = len(self.cellids)
        nnodes                = len(self.vtxids)
        ncellsmesh            = len(self.meshconnectivity)
        self.cellids          = Numeric.array(self.cellids,'i')
        self.vtxids           = Numeric.array(self.vtxids,'i')
        ivtx                  = Numeric.array([0],'i')
        ivtx                  = Numeric.resize(ivtx,[self.nnodesmesh])
        self.blckconnectivity = Numeric.array([0],'i')
        
        if ncells > 0:
            self.blckconnectivity = Numeric.resize(self.blckconnectivity,[ncells,nvc])


            ensight.getconnectivity(self.meshconnectivity,
                                    self.cellids,
                                    self.vtxids,
                                    self.blckconnectivity,
                                    ivtx,
                                    ncellsmesh,
                                    ncells,
                                    self.nnodesmesh,
                                    nnodes,
                                    nvc)
        else:
            self.blckconnectivity = Numeric.resize(self.blckconnectivity,[1,nvc])

        if self.debug:
            print >> self.fp, 'stats in getBlockIndicies.py'
            print >> self.fp, self.ID
            print >> self.fp, ncellsmesh
            print >> self.fp, len(self.cellblocks)
            print >> self.fp, len(self.cellids)
            print >> self.fp, len(self.vtxids)
        
        #catch any exceptions associated with getting the cellids and vtxids associated with this block
        AllZero = 'Coords all zero'

        def __nozeroes(cellids):
            T = Numeric.sometrue(Numeric.greater(cellids,0) | Numeric.less(cellids,0))
            "Accomodate Numpy/Numeric differences"
            try:
                if not T.any(): raise AllZero
            except:
                if not Numeric.greater(T,0): raise AllZero
        try:
            if len(self.cellids) > 0 :
                __nozeroes(self.cellids)
        except AllZero:
            print >> self.fp
            print >> self.fp, 'We will not continue extracting block connectivity for block %i ' %(self.ID)
            print >> self.fp, 'as the block cell IDs were not obtained.'
            print >> self.fp
            sys.exit(1)
        
    def str(self):
        print >> self.fp, 'for block:'
        print >> self.fp, self.ID
        print >> self.fp, self.nnodesmesh
        print >> self.fp, self.meshconnectivity
        print >> self.fp, 'global vertex ids'
        print >> self.fp, self.vtxids
        print >> self.fp, 'global cell ids'
        print >> self.fp, self.cellids
        print >> self.fp, 'local block connectivity'
        print >> self.fp, self.blckconnectivity

if __name__== '__main__':

    try:
       import numpy.oldnumeric as Numeric
    except ImportError:
       import Numeric
    except:
       raise
    
    fpwatch  = sys.stdout
    wdir     = os.getcwd()
    filename ='../../../scripts/test_TBrookParse/samples/map_output/map.TBrook.xml'
    
    from XMLgetStorageObject import getStorageObject

    storage   = getStorageObject(filename,fp=fpwatch,debug=1)	
    mesh      = storage.mlist[0] #assume default mesh

    #first get all data values associated with this mesh
    for meshspace in mesh.mslist:
        storage.getValues(meshspace.vlist)

    #having obtained the data values now fill in the following convenient dictionaries
    mesh.fillMesh()

    cblocks        = mesh.cells['BLOCKID']
    meshconn       = mesh.cells['VERTICES']
    blockid        = cblocks[Numeric.argmax(cblocks)]
    nnodesmesh     = mesh.vertices['N']
    try:
        blockinds  = getBlockIndicies(cblocks,blockid,nnodesmesh,meshconn,fp=fpwatch,debug=1)
        blockinds.getData()
    except:
        print >> fpwatch
        print >> fpwatch, 'Problems occurred in creating the indices from meshblock %i in %s' %(blockid,filename)
        print >> fpwatch, 'Exiting this component test'
        print >> fpwatch
        sys.exit(1)

    blockinds.str()


