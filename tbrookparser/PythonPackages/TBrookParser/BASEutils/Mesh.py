"""

 Mesh

 -----------------------------------------------------------------------------
  Purpose:

     Provides a generic base class for mesh objects.
     This object class is not intended to be directly used.
     It is used as the base class for (XML)getMesh and MODMesh classes.

  Public Interface(s):


  Contents
      class baseMesh
          __init__()
	  fillMesh()
          fillCells(meshspace)    
          fillVertices(meshspace) 
          fillFaces(meshspace)    
          fillEdges(meshspace)    
          str()                   
          debug(verbosity)

  Version:
    $ID$

  Author(s): Larry Cox (ljcox@lanl.gov)
 -----------------------------------------------------------------------------
"""
if __name__=='__main__':
    # Set sys.path for component testing mode
    import os, sys
    thisdir   = os.path.dirname(__file__)
    parsedir  = os.path.abspath(thisdir+'../')
    sys.path.append(parsedir)

class baseMesh:

    def __init__(self):
        "define the mesh object"

        self.name     = None
        self.type     = None # either HEX or TET
        self.mslist   = [] # mesh spaces that live on this mesh
        self.pdata    = {} # permute data associated with each mesh space
        self.updata   = {} # unpermute data associated with each mesh space
        self.cells    = {} # CELLS space variables and their data
        self.edges    = {} # EDGES space variables and their data
        self.faces    = {} # FACES space variables and their data
        self.vertices = {} # VERTICES space variables and their data

    def fillMesh(self):
        """
        fills in data for cells, vertices, faces and edge data, if present,
        for each meshspace on the mesh
        """
        for ms in self.mslist: 
            if ms.name == 'CELL':
                self.fillCells(ms)
            elif ms.name == 'VERTEX':
                self.fillVertices(ms)
            elif ms.name == 'FACE':
                self.fillFaces(ms)
            elif ms.name == 'EDGE':
                self.fillEdges(ms) 

    def fillCells(self,meshspace):
        "fills in data for variables that live on the meshspace mesh.cells" 

        if self.cells == {}:

            for variable in meshspace.vlist:
                self.cells[str(variable.name)] = variable.data
            self.cells['N'] = meshspace.size

    def fillVertices(self,meshspace):
        "fills in data for variables that live on the meshspace mesh.vertices" 

        if self.vertices == {}:

            for variable in meshspace.vlist:
                self.vertices[str(variable.name)] = variable.data
            self.vertices['N'] = meshspace.size

    def fillFaces(self,meshspace):
        "fills in data for variables that live on the meshspace mesh.faces" 

        if self.faces == {}:
            
            for variable in meshspace.vlist:
                self.faces[str(variable.name)] = variable.data
            self.faces['N'] = meshspace.size

    
    def fillEdges(self,meshspace):
        "fills in data for variables that live on the meshspace mesh.edges" 

        if self.edges == {}:

            for variable in meshspace.vlist:
                self.edges[str(variable.name)] = variable.data
            self.edges['N'] = meshspace.size

            
    def metastr(self):
        "prints out meta values of attributes in the mesh object"    

        print self.name
        print self.type
        print '\n meshspace characteristics \n'
        for mshspc in self.mslist:
            mshspc.metastr()

    def str(self):
        "prints meta + data values for all attributes in the mesh object"
        
        self.fillMesh()
        print self.cells
        print self.vertices
        print self.faces
        print self.edges
        self.metastr()

    def debug(self,verbosity=0):

        if verbosity > 0:
            self.str()
        else            :
            self.metastr()

if __name__=='__main__':


    from PYTHONutils import uTestOpts
    opts = uTestOpts('fd')
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    try:
       print "Incomplete testing package"

       "Put testing code here"

    except:
        print "---> Test failed in some aspect <---"
        print "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise
