"""

 MeshSpace

 -----------------------------------------------------------------------------
  Purpose:

     Provides a generic base class for meshspace object.
     This object class is not intended to be directly used.
     It is used as the base class for (XML)getMeshSpace and MODMeshSpace classes.

  Public Interface(s):

    None

  Contents
      class baseMeshSpace
          __init__()
	  metastr()
	  str()
	  debug(verbosity)

  Version:
    $ID$

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

class baseMeshSpace:

    def __init__(self):
        "define the meshspace object"

        self.name           = None
        self.size           = None
        self.permutation    = None
        self.unpermutation  = None
        self.vlist          = []

    def metastr(self):
        "prints out meta values of attributes in the timesteps object"    

        print '\n meshspace attributes \n'
        print self.name
        print self.size
	print '\n variables in this meshspace :\n'
	if self.permutation   != None: self.permutation.str()
	if self.unpermutation != None: self.unpermutation.str()
        for thisvar in self.vlist:
            thisvar.str()

    def str(self):
        "prints out meta + data values of attributes in the timesteps object"

	self.metastr()
	if self.permutation   != None: self.permutation.str()
	if self.unpermutation != None: self.unpermutation.str()
	for var in self.vlist:
	    var.getData()
	    var.str()
        #fillValues(self.vlist)
        #getValues(self.vlist)
        
    def debug(self,verbosity=0):

        if verbosity > 0:
            self.str()
        else:
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
    
