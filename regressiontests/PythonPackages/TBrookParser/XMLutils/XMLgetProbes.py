"""

 XMLgetProbes

 -----------------------------------------------------------------------------
  Purpose:
  
     To instantiate and provide methods for getProbes objects.
  
  Public Interface(s):
  
     p = getProbes(filename,spacelu,thedoc)
     p.str()

  Contains:
    class getProbes
        __init__(self,filename,spacelu,thedoc,parse=MiniDom,debug=0)
        str(self)

    Unit Test Block
  
  Author(s): Sharen Cummins (scummins@lanl.gov)
 -----------------------------------------------------------------------------
"""
import sys
if __name__=='__main__':
    # Set sys.path for component testing mode
    import os
    thisdir    = os.path.dirname(__file__)
    parserdir  = os.path.abspath(thisdir+'../')
    sys.path.append(parserdir)
    
from XMLutils       import MiniDom
from XMLgetProbe    import getProbe

class getProbes:

    def __init__(self,filename,spacelu,thedoc,fp=sys.stdout,parse=MiniDom,debug=0):

        "initialise Probes object"

        self.nameslist  = []
        self.plist      = []
        self.pointer    = None
        self.fp         = fp

        elems           = parse().getElements(thedoc,'PROBE')

        for element in elems:

            name        = str(parse().getElemAttribute(element,'NAME'))

            "create a list of probes "
            probe       = getProbe(element,filename,spacelu,self.fp,parse)

            self.plist.append(probe)
            self.nameslist.append(name)


    def str(self):

        for thisprobe in self.plist:
            thisprobe.str()
        
if __name__=='__main__':

    'testing this component'

    from PYTHONutils import uTestOpts
    opts = uTestOpts('fd',
		     defaults = {'f':'../../../../regressiontests/sensitivity_stefan/stefan_golden/stefan.TBrook.xml'} )
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)
    
    try:
        parse     = MiniDom
        node      = parse().getDoc(opt.f)
        fnode     = parse().getElement(node,'FILE')
        spacelu   = {None:'CELL','NONE':'NONE',
                     'CELL':'CELL', 'CELLS':'CELL',
                     'FACES':'FACE', 'FACE':'FACE',
                     'VERTICES':'VERTEX', 'VERTEX':'VERTEX',
                     'NODES':'VERTEX','NODE':'VERTEX',
                     'EDGES':'EDGE','EDGE':'EDGE',
                     'MESH':'MESH'}

        xmlprobes = getProbes(opt.f,spacelu,node)
        if xmlprobes.plist:
            print len(xmlprobes.plist), 'probes found in file:',\
                  os.path.basename(opt.f)
            xmlprobes.str()

        else:
            print "No probes found."
        
    except:
        print "---> Test failed in some aspect <---"
        print "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise
