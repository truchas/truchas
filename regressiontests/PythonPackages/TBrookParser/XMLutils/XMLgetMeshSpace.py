"""

 XMLgetMeshSpace

 -----------------------------------------------------------------------------
  Purpose:
  
     To instantiate and provide methods for a MeshSpace object.
  
  Public Interface(s):
  
     ms = getMeshSpace(srcPointer)
     ms.str()

  Contains:
    class getMeshSpace
        __init__(self,srcPointer,name='CELL',parse=MiniDom)
        str(self)

    Unit Test Block
  
  Author(s): Sharen Cummins (scummins@lanl.gov)
             Erin Iesulauro Barker (eibarker@lanl.gov)
 -----------------------------------------------------------------------------
"""

if __name__=='__main__':
    # Set sys.path for component testing mode
    import os, sys
    thisdir    = os.path.dirname(__file__)
    parserdir  = os.path.abspath(thisdir+'../')
    sys.path.append(parserdir)

from BASEutils      import baseMeshSpace
from XMLutils       import MiniDom, getVariable

class getMeshSpace(baseMeshSpace):

    def __init__(self,srcPointer,name='CELL',parse=MiniDom):

        "Initialize the MeshSpace object"
        baseMeshSpace.__init__(self)

        "now set values"

        self.name  = name
        
        meshnode   = srcPointer['node']
        spacelu    = srcPointer['spacelu']
        elems      = parse().getElements(meshnode,'FILEVAR')
 
        for element in elems:
            varname     = parse().getElemAttribute(element,'Name')
            thispace    = parse().getElemAttribute(element,'Map')
            if thispace == 'NONE':
                thispace = parse().getElemAttribute(element,'INFORMATION')  
            if spacelu.has_key(thispace):
                thispace = spacelu[thispace]
            dformat     = srcPointer['datafmt']
            tfile       = srcPointer['tfile']
            if thispace == name:
                if varname == 'UNPERMUTE':
                    "get actual unpermute mappings"
                    self.unpermutation = getVariable(srcPointer,element)
                    self.unpermutation.getData()
                elif varname == 'PERMUTE':
                    "get actual permute mappings"
                    self.permutation = getVariable(srcPointer,element)
                    self.permutation.getData()
                elif varname == 'N':
                    tmp        = getVariable(srcPointer,element)
                    tmp.getData()
                    self.size = tmp.data[0]
                else :
                    thisvar    = getVariable(srcPointer,element)
                    self.vlist.append(thisvar) #required for fields        

    
if __name__=='__main__':

    'for testing this component'

    from PYTHONutils    import getpath
    from FORTRANutils   import fortransupport
    from XMLutils       import BrookUtilities

    from PYTHONutils import uTestOpts
    opts = uTestOpts('fd')
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    try:
        parse      = MiniDom
        node       = parse().getDoc(opt.f)
        nodes      = parse().getElements(node,'MESH')
        for node in nodes:
            filenode     = parse().getElement(node,'FILE')
            meshfile     = str(parse().getValue(filenode))
            corrmeshfile = getpath(opt.f,meshfile)
        node     = parse().getDoc(corrmeshfile)
        meshnode = parse().getElement(node,'MESH')
        filenode = parse().getElement(meshnode,'FILE')
        file     = str(parse().getValue(filenode))
        corrfile = getpath(opt.f,file)
        file       = BrookUtilities().fixFileName(corrfile)
        tfile      = fortransupport.T_Binary_File(file)
    
        spacelu = {None:'CELL','NONE':'NONE',
                   'CELL':'CELL', 'CELLS':'CELL',
                   'FACES':'FACE', 'FACE':'FACE',
                   'VERTICES':'VERTEX', 'VERTEX':'VERTEX',
                   'NODES':'VERTEX','NODE':'VERTEX',
                   'EDGES':'EDGE','EDGE':'EDGE',
                   'MESH':'MESH'}

        pointer    = {'node':meshnode, 'datafmt':'binary',
                      'tfile':tfile, 'spacelu':spacelu}
        
        for map in ['CELL','VERTEX','FACE']:
            space  = getMeshSpace(pointer,map)
            space.metastr()

    except:
        print "---> Test failed in some aspect <---"
        print "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise
