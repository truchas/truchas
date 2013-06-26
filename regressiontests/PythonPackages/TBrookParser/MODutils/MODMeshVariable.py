"""

 MODMeshVariable

 -----------------------------------------------------------------------------
  Purpose:

    Initialization for a modMeshvariable object.
    Only stores data and provide a simple 'str' method

  Public Interface(s):

    mMV = modMeshVariable(name,rank,shape,type,data,mesh,meshspace)
    mMV.str()
  
  Contains:
    class modMeshVariable
        __init__(self,name,rank,shape,type,data,mesh,meshspace,debug=0)
        str(self)

    ---> NO Unit Test Block

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

from BASEutils import baseVariable

class modMeshVariable(baseVariable):

    def __init__(self,name,rank,shape,type,data,mesh,meshspace,debug=0):

        "initialise the variable object"
        baseVariable.__init__(self)

        "Modify the MeshVariable object"
        
        self.name          = name
        self.nickname      = name
        self.rank          = rank
        self.shape         = shape
        self.data          = data
        self.type          = type
        self.mesh          = mesh      # mesh name      - not the object
        self.meshspace     = meshspace # meshspace name - not the object

        

        
        


        
