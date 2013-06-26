"""

 createFileObject

 -----------------------------------------------------------------------------
  Purpose:
  
     Instantiates an empty FileObject
  
  Public Interface(s):
   
     fo = createFileObject

  Contains:
    class createFileObject
        __init__(self)

    ---> NO Unit Test Block
  
  Author(s): 
 -----------------------------------------------------------------------------
"""

class createFileObject:

    def __init__(self):

        "initialise the file object"

        self.file           = None      #name of the file used to create this file object
        self.storage        = None      #storage object created from this file
        self.regions        = []        #list of regions defined to live on the meshes in this file
