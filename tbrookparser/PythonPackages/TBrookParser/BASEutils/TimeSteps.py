"""

 TimeSteps

 -----------------------------------------------------------------------------
  Purpose:

     Provides a generic base class for timesteps object.
     This object class is not intended to be directly used.
     It is used as the base class for (XML)getTimeSteps and MODTimeSteps classes.

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

class baseTimeSteps:

    def __init__(self):
        "define the timesteps object"

        self.idlist     = []
        self.timeslist  = []
        self.tlist      = []
        self.cyclelist  = []
        self.vlist      = {}
        self.clist      = {}


    def metastr(self):
        "prints out meta values of attributes in the timesteps object"    

        print '\n timesteps attributes \n'
        print self.idlist
        print self.timeslist
        print self.cyclelist
        for thistime in self.tlist:
            thistime.metastr()

    def str(self):
        "prints out meta + data values of attributes in the timesteps object"

	self.metastr()
        for thistime in self.tlist:
            for var in thistime.vlist:
		var.getData()
		var.str()
            #getValues(thistime.vlist)
            #fillValues(thistime.vlist)
        

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
