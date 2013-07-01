"""

 XMLgetProbe

 -----------------------------------------------------------------------------
  Purpose:
  
     To instantiate and provide methods for getProbe objects.
  
  Public Interface(s):
  
     p = getProbe(node,TBrookfilename,spacelu)
     p.tossData()
     p.str()

  Contains:
    class getProbe
        __init__(self,node,TBrookfilename,spacelu,parse=MiniDom)
        tossData(self)
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
    thisdir    = os.path.dirname(__file__)
    parserdir  = os.path.abspath(thisdir+'../')
    sys.path.append(parserdir)

from BASEutils      import baseProbe
from XMLutils       import MiniDom, BrookUtilities, getVariable
from PYTHONutils    import getpath
from FORTRANutils   import fortransupport


class getProbe(baseProbe):
        
    def __init__(self,node,TBrookfilename,spacelu,fp=sys.stdout,parse=MiniDom):

        "initialise the probe object"
	baseProbe.__init__(self)

        #self.name        = None
        #self.description = None
        #self.coords      = []
        #self.node        = {} # for id, coords of nearest node to this probe
        #self.cell        = {} # for id, coords of nearest cell to this probe
        #self.vlist       = [] # a list of variable objects defined in this probe
        self.__pointer     = None
        self.__fp          = fp

        try:
            self.name    = str(parse().getElemAttribute(node,'NAME'))
            self.coords.append(float(parse().getElemAttribute(node,'X')))
            self.coords.append(float(parse().getElemAttribute(node,'Y')))
            self.coords.append(float(parse().getElemAttribute(node,'Z')))
        except:
            print >> self.__fp
            print >> self.__fp, "Error this probe is not a standard probe as",
            print >> self.__fp, "a probe name and coordinates are not provided."
            print >> self.__fp, "We will continue to postprocess",
            print >> self.__fp, "but will not include this probe in further diagnostics"
            print >> self.__fp
            sys.exit()

        try:
            descnode         = parse().getElement(node,'DESCRIPTION')
            self.description = str(parse().getValue(descnode))
        except:
            self.description = 'No description for probe %s' %(self.name)
            
        try:
            ndenode              = parse().getElement(node,'NODE')
            self.node['ID']      = int(parse().getElemAttribute(ndenode,'ID'))
            self.node['COORDS']  = []

            a = float(parse().getElemAttribute(ndenode,'X'))
            self.node['COORDS'].append(a)

            a = float(parse().getElemAttribute(ndenode,'Y'))
            self.node['COORDS'].append(a)

            a = float(parse().getElemAttribute(ndenode,'Z'))
            self.node['COORDS'].append(a)
            
            cllnode              = parse().getElement(node,'CELL')
            self.cell['ID']      = int(parse().getElemAttribute(cllnode,'ID'))
            self.cell['COORDS']  = []

            a = float(parse().getElemAttribute(cllnode,'X'))
            self.cell['COORDS'].append(a)

            a = float(parse().getElemAttribute(cllnode,'Y'))
            self.cell['COORDS'].append(a)

            a = float(parse().getElemAttribute(cllnode,'Z'))
            self.cell['COORDS'].append(a)

        except:
            print >> self.__fp
            print >> self.__fp, "No probe cell and node information available"
            print >> self.__fp
            
        try:
            varnodes       = parse().getElements(node,'PROBEVAR')
            for vnode in varnodes:
                varname         = parse().getElemAttribute(vnode,'NAME')

                filenode        = parse().getElement(vnode,'FILE')
                datafmt         = parse().getElemAttribute(filenode,'FORMAT')
                file            = str(parse().getValue(filenode))
                assert datafmt == 'binary'
                """
                given path of filename, ensure actual path of binary
                file is correct
                """
                corrfile        = getpath(TBrookfilename,file) 
                file            = BrookUtilities().fixFileName(corrfile)
           
                Estring = 'import error'
                from POSTPROCESSORutils.dataSource import dataSource
                Estring = 'getting data source'
                tfile   = dataSource(file,'f95b',debug=0)

                if tfile.status > -1 :
                    #look aside file exists..
                    Estring         = 'setting pointer'
                    self.__pointer    = {'node':vnode,'datafmt':datafmt,
                                       'tfile':tfile,'spacelu':spacelu}
                    Estring         = 'getting elements'
                    elems           = parse().getElements(vnode,'FILEVAR')
                    Estring         = 'looping'
                    for element in elems:
                        thisvar     = getVariable(self.__pointer,element,parse)
                    self.vlist.append(thisvar)
                else:
                    print >> self.__fp, 'Caution data file associated with probe',
                    print >> self.__fp, '"%s" and variable "%s"' %(str(self.name),str(varname))
                    print >> self.__fp, 'is not available. Will not include this variable in any',
                    print >> self.__fp, 'postprocessing operations associated with this probe. '
                    print >> self.__fp
        except:
            print >> self.__fp
            print >> self.__fp, 'No data exists in probe %s' %(str(self.name))
            print >> self.__fp, 'Will not include this probe in any',
            print >> self.__fp, 'further postprocessing operations'
            print >> self.__fp


    def tossData(self):
        """
        Release all pointers to data in self.vlist
        """
        for v in self.vlist:
            v.data = None

        if(self.__pointer != None):
            if(self.__pointer.has_key('tfile')):
                if self.__pointer['tfile'] != None:
                    self.__pointer['tfile'].close()
            self.__pointer['tfile'] = None
        return

            
if __name__=='__main__':
    
    'for testing this component'

    from PYTHONutils import uTestOpts
    opts = uTestOpts('fd',
		     defaults = {'f':'../../../../regressiontests/sensitivity_stefan/stefan_golden/stefan.TBrook.xml'} )
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)
    
    try:
        spacelu = {None:'CELL','NONE':'NONE', 'CELL':'CELL',
                   'CELLS':'CELL', 'FACES':'FACE', 'FACE':'FACE',
                   'VERTICES':'VERTEX', 'VERTEX':'VERTEX',
                   'NODES':'VERTEX','NODE':'VERTEX',
                   'EDGES':'EDGE','EDGE':'EDGE', 'MESH':'MESH'}
        parse   = MiniDom
        node    = parse().getDoc(opt.f)
        probes  = parse().getElements(node,'PROBE')
        
        if probes:
            myprobe  = probes[-1] # take last one 
            xmlprobe = getProbe(myprobe,opt.f,spacelu)
            print 'probe:'
	    # eib - swap str for metastr if you want to see the probe data
            #xmlprobe.metastr()
            xmlprobe.metastr()
            
        else:
            print "No probes found."

    except:
        print "---> Test failed in some aspect <---"
        print "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise
