"""

 Aborts

 -----------------------------------------------------------------------------
  Purpose:

     Provides a generic base class for aborts object.
     This object class is not intended to be directly used.
     It is used as the base class for the (XML)getAborts class.

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

class baseAborts:

    def __init__(self):
        "define the aborts object"

        self.lrlist     = [] # linear residual list 
        self.nlrlist    = [] # non-linear residual list 
             

    def metastr(self):
        "prints out meta values of attributes in the aborts object"    

        print '\n aborts attributes \n'

        print 'linear residual aborts:'
        for thislr in self.lrlist:
            thislr.str()
        print 'non-linear residual aborts:'
        for thisnlr in self.nlrlist:
            thisnlr.str()
        
    def str(self):
        "prints out meta + data values of attributes in the aborts object"

        self.metastr()
	for thislr in self.lrlist:
	    thislr.vlist.getData()
	for thisnlr in self.nlrlist:
	    thisnlr.vlist.getData()
	"""
	# -- eib: fillVales isn't defined anywhere!!
        for thislr in self.lrlist:
            fillValues(thislr.vlist)

        for thisnlr in self.nlrlist:
            fillValues(thisnlr.vlist)
	"""


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
