"""

 CREATEVariable

 -----------------------------------------------------------------------------
  Purpose:
  
     Creates mapped variable objects.
  
  Public Interface(s):
    v = createVariable(name,nickname,rank,shape,
                       type,data,mesh,meshspace)
    v.str()
    
  Contains:
    class createVariable
        __init__(self,name,nickname,rank,shape,
                 type,data,mesh,meshspace,
                 mustmap=0,maprule='constants_preserving',debug=0)
        str(self)

    ---> NO Unit Test Block
  
  Version:
    $ID$
    
  Author(s): Sharen Cummins (scummins@lanl.gov)
 -----------------------------------------------------------------------------
"""

class createVariable:

    def __init__(self,name,nickname,rank,shape,
                 type,data,mesh,meshspace,
                 mustmap=0,maprule='constants_preserving',debug=0):

        "create a mapped variable object"

        self.name      = name
        self.nickname  = nickname
        self.rank      = rank
        self.shape     = shape
        self.data      = data
        self.type      = type
        self.mesh      = mesh      # mesh name      - not the object itself
        self.meshspace = meshspace # meshspace name - not the object itself
        self.mustmap   = mustmap
        self.maprule   = maprule
        self.oldintgrl = 0.0
        self.newintgl  = 0.0

    def metastr(self):
        print self.name, self.nickname, self.rank, \
              self.shape, self.type, self.mesh, \
              self.meshspace, self.mustmap

    def str(self):

	print self.metastr()
        print self.data
