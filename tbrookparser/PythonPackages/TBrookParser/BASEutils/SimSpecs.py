"""

 SimSpecs

 -----------------------------------------------------------------------------
  Purpose:

     Provides a generic base class for simspecs objects.
     This object class is not intended to be directly used.
     It is used as the base class for (XML)getSimSpecs and MODSimSpecs classes.

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

class baseSimSpecs:

    def __init__(self):
        "define the simspecs object"

        self.fileid      = 'TRF-3' # 3 - no sideset information included
        self.nsspecs     = 0  # number of simulation specifications
        self.sspecs      = [] # list of simulation specifications descriptions
        self.nfeats      = 0  # number of features
        self.feats       = [] # list of feature descriptions
        self.nphases     = 0  # number of phases
        self.ncomps      = 0  # number of components in each phase
        self.csf         = 1  # mesh coordinate scale factor
        self.tbrookvers  = 1  # XML output version number
	self.nspecies    = 0  # number of species concentrations present for diffusion

    def str(self):

        print self.fileid
        print self.nsspecs
        print self.sspecs
        print self.nfeats
        print self.feats
        print self.nphases
        print self.ncomps
        print self.csf, type(self.csf)
	print self.nspecies


    def debug(self):
        """
        for debugging this object
        """

        self.str()

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
