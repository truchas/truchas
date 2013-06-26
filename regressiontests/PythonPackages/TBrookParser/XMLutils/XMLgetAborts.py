"""

 XMLgetAborts

 -----------------------------------------------------------------------------
  Purpose:
  
      Provide class and initialization for Aborts object
  
  Public Interface(s):

      a = getAborts(filename,spacelu,thedoc)
      a.str()

  Contains:
    class getAborts
        __init__(self,filename,spacelu,thedoc,parse=MiniDom)
        str(self)

    Unit Test Block
  
  Author(s): Sharen Cummins (scummins@lanl.gov)
             Erin Iesulauro Barker (eibarker@lanl.gov)
 -----------------------------------------------------------------------------
"""
import sys
if __name__=='__main__':
    # Set sys.path for component testing mode
    import os
    thisdir   = os.path.dirname(__file__)
    parsedir  = os.path.abspath(thisdir+'../')
    sys.path.append(parsedir)

from BASEutils      import baseAborts
from XMLutils       import MiniDom
from XMLgetTimeStep import getTimeStep

class getAborts(baseAborts):

    def __init__(self,filename,spacelu,thedoc,fp=sys.stdout,parse=MiniDom):

        "initialise Aborts object"
	baseAborts.__init__(self)

        self.__pointer    = None
        self.__fp         = fp

        lrelems  = parse().getElements(thedoc,'LINEAR_RESIDUAL')

        for elem in lrelems:
            id          = int(parse().getElemAttribute(elem,'SEQ'))
            time        = 0.0

            "create a list of linear residual structures "
            thislr      = getTimeStep(elem,None,filename,spacelu,self.__fp,parse)

            for var in thislr.vlist:
                var.nickname = var.name
                
            self.lrlist.append(thislr)

        nlrelems  = parse().getElements(thedoc,'NONLIN_RESIDUAL')

        for elem in nlrelems:
            id          = int(parse().getElemAttribute(elem,'SEQ'))
            time        = 0.0

            "create a list of nonlinear residual structures "
            thisnlr     = getTimeStep(elem,None,filename,spacelu,self.__fp,parse)

            for var in thisnlr.vlist:
                var.nickname = var.name

            self.nlrlist.append(thisnlr)

    """
    def str(self):

	self.metastr()
	for thislr in self.lrlist:
	    thislr.vlist.getData()
	for thisnlr in self.nlrlist:
	    thisnlr.vlist.getData()
    """

if __name__=='__main__':
    
    'testing this component'

    from PYTHONutils import uTestOpts
    opts = uTestOpts('fd')
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    try:
        parse       = MiniDom
        node        = parse().getDoc(opt.f)
        
        spacelu     = {None:'CELL','NONE':'NONE',
                       'CELL':'CELL', 'CELLS':'CELL',
                       'FACES':'FACE', 'FACE':'FACE',
                       'VERTICES':'VERTEX', 'VERTEX':'VERTEX',
                       'NODES':'VERTEX','NODE':'VERTEX',
                       'EDGES':'EDGE','EDGE':'EDGE',
                       'MESH':'MESH'}

        xmlaborts = getAborts(opt.f,spacelu,node,parse)
        nresid    = len(xmlaborts.lrlist) + len(xmlaborts.nlrlist)
        if nresid:
            xmlaborts.metastr()
        else:
            print "No abort residuals found."
        
    except:
        print "---> Test failed in some aspect <---"
        print "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise

    

    



        

