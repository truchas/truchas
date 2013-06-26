"""

 Variable

 -----------------------------------------------------------------------------
  Purpose:

     Provides a generic base class for variable object.
     This object class is not intended to be directly used.
     It is used as the base class for
     -(XML)getVariable,
     -MODMeshVariable,
     -MODFieldVariable,
     -MAPVariable
     classes.

  Public Interface(s):

    None

  Contents

  Version:
    $ID$

  Author(s): Sharen Cummins (scummins@lanl.gov)
 -----------------------------------------------------------------------------
"""
if __name__=='__main__':
    # Set sys.path for component testing mode
    import os, sys
    thisdir   = os.path.dirname(__file__)
    parsedir  = os.path.abspath(thisdir+'../')
    sys.path.append(parsedir)

class baseVariable:

    def __init__(self):
        "define the variable object"

        self.name          = None
        self.nickname      = None
        self.rank          = None
        self.shape         = None
        self.data          = None
        self.type          = None
        self.mesh          = None    # mesh name      - not the object
        self.meshspace     = None    # meshspace name - not the object


    def getData(self):
        "retrieves data for the variable.data attribute"

    def metastr(self):
        "prints out meta values of attributes in the variable object"    

        print '\n variable attributes \n'
        print self.name, self.nickname, self.rank, self.shape,\
              self.type, self.mesh, self.meshspace

    def str(self):
        "prints out meta + data values of attributes in the timesteps object"

        self.metastr()
	self.getData()
        print self.data

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
