"""

 Probe

 -----------------------------------------------------------------------------
  Purpose:

     Provides a generic base class for probe object.
     This object class is not intended to be directly used.
     It is used as the base class for the (XML)getProbe class.

  Public Interface(s):

    None

  Contents

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

class baseProbe:

    def __init__(self):
        "define the probe object"

        self.name          = None
        self.description   = None
        self.coords        = [] # coords of the probe location 
        self.node          = {} # id, coords of nearest node to this probe
        self.cell          = {} # id, coords of nearest cell to this probe
        self.vlist         = [] # a list of variable objects defined in this probe

    def tossData(self):
        "remove all vlist data from this probe"

    def metastr(self):
        "prints out meta values of attributes in the probe object"    

        print '\n probe attributes \n'
        print self.name
        print self.description
        print self.coords
        print self.node
        print self.cell
        print 'this probe variable attributes'
        for var in self.vlist:
            var.metastr() 
        
    def str(self):
        "prints out meta + data values of attributes in the probe object"
        
        self.metastr()
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
