"""

 XMLgetMesh

 -----------------------------------------------------------------------------
  Purpose:
  
     To instantiate and provide methods for Mesh Objects.
  
  Public Interface(s):

     m = getMesh(filename,spacelu)
     m.fillCells(meshspace)    
     m.fillVertices(meshspace) 
     m.fillFaces(meshspace)    
     m.fillEdges(meshspace)    
     m.str()                   

  Contains:
    class getMesh
        __init__(self,filename,spacelu,parse=MiniDom)
        str(self)                     

    From baseMesh class, it inherits
        fillMesh(self)
        fillCells(self,meshspace)
        fillVertices(self,meshspace)
        fillFaces(self,meshspace)
        fillEdges(self,meshspace)

    Unit Test Block
  
  Author(s): Sharen Cummins (scummins@lanl.gov)
             Erin Iesulauro Barker (eibarker@lanl.gov)
 -----------------------------------------------------------------------------
"""

if __name__=='__main__':
    # Set sys.path for component testing mode
    import os, sys
    thisdir   = os.path.dirname(__file__)
    parsedir  = os.path.abspath(thisdir+'../')
    sys.path.append(parsedir)

from BASEutils          import baseMesh
from PYTHONutils        import getpath,unique,SystemExit
from FORTRANutils       import fortransupport
from XMLutils           import MiniDom,BrookUtilities,getVariable,getMeshSpace

class getMesh(baseMesh):

    def __init__(self,filename,spacelu,parse=MiniDom):

        "initialise the mesh object"
        baseMesh.__init__(self)

        "now start defining the attributes in the mesh object"

        try:
            node            = parse().getDoc(filename)
        except:
            print 'Unable to parse the mesh XML file '
            print filename
            print 
            print 'Please ensure your TBrook.xml and DefaultMesh.xml files'
            print 'live in the same directory.'
            print 
            print 'Also please check the XML file is formed correctly.'
            SystemExit()
            
        meshnode        = parse().getElement(node,'MESH')
        meshfmt         = parse().getElemAttribute(meshnode,'FORMAT')
        filenode        = parse().getElement(meshnode,'FILE')
        datafmt         = parse().getElemAttribute(filenode,'FORMAT')
        file            = str(parse().getValue(filenode))
        
        self.name       = parse().getElemAttribute(meshnode,'LABEL')
        
        assert datafmt == 'binary'

        #first ensure path of binary data file is correct

        #given path of filename, ensure actual path of meshfile is correct
        corrfile        = getpath(filename,file)
        
        file            = BrookUtilities().fixFileName(corrfile)

        try:
            from POSTPROCESSORutils.dataSource import dataSource
            tfile        = dataSource(file,'f95b',debug=0)
        except:
            print
            print 'Cannot read binary mesh look aside file'
            print file
            SystemExit()
            
        elems        = parse().getElements(meshnode,'FILEVAR')
        self.pointer = {'node':meshnode, 'datafmt':datafmt,
                        'tfile':tfile, 'spacelu':spacelu}
        
        tmplist      = []
        for element in elems:
            name        = parse().getElemAttribute(element,'Name')
            if name == 'MESHTYPE':
                tmp       = getVariable(self.pointer,element)
                tmp.getData()
                if (tmp.data[0] == 'T'):
                    self.type = 'TET'
                else:
                    self.type = 'HEX'
            space       = parse().getElemAttribute(element,'Map')
            if spacelu.has_key(space):
                space   = spacelu[space]
            tmplist.append(space.encode())

        tmplist = unique(tmplist)

        for i in range(len(tmplist)):
            thispace  = getMeshSpace(self.pointer,tmplist[i])
            if thispace.name != 'NONE':
                self.mslist.append(thispace)
            if thispace.unpermutation != None:
                self.updata[thispace.name] = thispace.unpermutation.data
            else:
                self.updata[thispace.name] = []
            if thispace.permutation != None:
                self.pdata[thispace.name] = thispace.permutation.data
            else:
                self.pdata[thispace.name] = []

    """
    def str(self):

        print self.name
        print self.type

        print 'meshspace characteristics'
        for mshspc in self.mslist:
            mshspc.str()
    """

if __name__=='__main__':
    
    'for testing this component'

    from PYTHONutils import uTestOpts
    opts = uTestOpts('fd')
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)
    
    try:
        parse        = MiniDom
        node         = parse().getDoc(opt.f)
        node         = parse().getElement(node,'MESH')
        node         = parse().getElement(node,'FILE')
        
        meshfile     = str(parse().getValue(node))
        corrmeshfile = getpath(opt.f,meshfile)
    
        spacelu = {None:'CELL','NONE':'NONE',
                   'CELL':'CELL', 'CELLS':'CELL',
                   'FACES':'FACE', 'FACE':'FACE',
                   'VERTICES':'VERTEX', 'VERTEX':'VERTEX',
                   'NODES':'VERTEX','NODE':'VERTEX',
                   'EDGES':'EDGE','EDGE':'EDGE',
                   'MESH':'MESH'}

        xmlmesh = getMesh(corrmeshfile,spacelu,parse)

        if len(xmlmesh.mslist)>0:
            xmlmesh.metastr()
        else:
            print "No meshspace(s) found."

    except:
        print "---> Test failed in some aspect <---"
        print "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise
        


