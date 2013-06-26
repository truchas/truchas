"""

 MODSimSpecs

 -----------------------------------------------------------------------------
  Purpose:

     An object class to store static SimSpec information

     Note: This is a candidate for a subClass of a generic SimSpec
           that only HOLDS data.  The function of populating the object, as
           in ../XMLutils/XMLgetSimSpec, is not included.

           The base class would only define empty attributes and
           perhaps the simple str() method as defined in this class.

  Public Interface(s):

    mSS = modSimSpecs(TStorageSpecs)
    mSS.str()

  Contains:
    class modSimSpecs
        __init__(self,TStorageSpecs,debug=0)
        str(self)

    Unit Test Block

  Version:
    $ID$

  Author(s): Sharen Cummins (scummins@lanl.gov)
             Erin Iesulauro Barker (eibarker@lanl.gov)
 -----------------------------------------------------------------------------
"""

if __name__=='__main__':
    " Set sys.path for component testing mode "
    import os, sys
    thisdir   = os.path.dirname(__file__)
    # modify dot-dots appropriately
    parsedir  = os.path.abspath(thisdir+'../') 

    sys.path.append(parsedir)

from BASEutils       import baseSimSpecs

class modSimSpecs(baseSimSpecs):

    def __init__(self,TStorageSpecs,debug=0):

        "initialise the simspecs object"
        baseSimSpecs.__init__(self)

        self.fileid      = TStorageSpecs.fileid
        self.nsspecs     = TStorageSpecs.nsspecs
        self.sspecs      = TStorageSpecs.sspecs
        self.nfeats      = TStorageSpecs.nfeats
        self.feats       = TStorageSpecs.feats
        self.nphases     = TStorageSpecs.nphases
        self.ncomps      = TStorageSpecs.ncomps
        #self.__parse     = TStorageSpecs.parse
        self.__debug     = debug
        
if __name__=='__main__':

    # test modifying simspecs structure
    from PYTHONutils import uTestOpts
    dfltdir = '../../../scripts/test_TBrookParse/samples/'
    opts = uTestOpts('fd',
		     dir=dfltdir,
		     defaults = {'f' : dfltdir + 'map_output/map.TBrook.xml'})
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    try:

        "first get an existing XML storage objec"
        
        from XMLgetStorageObject import getStorageObject
        tstorage     = getStorageObject(opt.f,debug=opt.d)
        
        "then make a modSimSpecs object from its specs"

        modsimspecs  = modSimSpecs(tstorage.specs, debug=opt.d)
        modsimspecs.str()
        
    except:
        print "---> Test failed in some aspect <---"
        print "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise



