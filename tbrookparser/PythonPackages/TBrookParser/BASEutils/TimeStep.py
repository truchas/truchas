"""

 TimeStep

 -----------------------------------------------------------------------------
  Purpose:

     Provides a generic base class for a timestep object.
     This object class is not intended to be directly used.
     It is used as the base class for (XML)getTimeStep and MODTimeStep classes.

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

class baseTimeStep:

    def __init__(self):
        "define a timestep object"

        self.id         = None
        self.cycle      = None
        self.emseq      = None
        self.emcycle    = None
        self.emfile     = None
        self.time       = None
        self.dt         = None
        self.dtcourant  = None
        self.vlist      = [] 
	self.itercount  = None


    def tossData(self):
        "remove all vlist data from this timestep" 


    def metastr(self):
        "prints out meta values of attributes in the timestep object"    

        print '\n timestep attributes \n'
        print self.id
        print self.cycle
        print self.time
        print self.dt
        print self.dtcourant
        print self.itercount
        print '\n this timestep variables attributes \n'
        for var in self.vlist:
            var.metastr() 

    def str(self):
        "prints out meta + data values of attributes in the timesteps object"

	self.metastr()
        #fillValues(self.vlist)
        #getValues(self.vlist)
	for var in self.vlist: 
	    var.getData()
	    var.str()

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
