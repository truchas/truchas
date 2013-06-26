"""

 XMLgetTimeStep

 -----------------------------------------------------------------------------
  Purpose:
  
    To instantiate and provide methods for TimeStep objects.
  
  Public Interface(s):

    ts = getTimeStep(node,emnode,TBrookfilename,spacelu)
    ts.tossData()

  Contains:
    class getTimeStep
        __init__(self,node,emnode,TBrookfilename,spacelu,parse=MiniDom)
        tossData(self)

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

from BASEutils      import baseTimeStep
from XMLutils       import MiniDom, BrookUtilities, getVariable
from PYTHONutils    import getpath
from FORTRANutils   import fortransupport 

class getTimeStep(baseTimeStep):
        
    def __init__(self,node,emnode,TBrookfilename,spacelu,fp=sys.stdout,parse=MiniDom):

        "initialise the timestep object"
	baseTimeStep.__init__(self)

        self.__pointer    = None
        self.__fp         = fp

        try:
            self.cycle      = int(parse().getElemAttribute(node,'CYCLE'))
            self.time       = float(parse().getElemAttribute(node,'t'))
            self.dt         = float(parse().getElemAttribute(node,'dt'))
            self.dtcourant  = float(parse().getElemAttribute(node,'dt_courant'))
            self.id         = int(parse().getElemAttribute(node,'ID'))
        except:
            filenode = parse().getElement(node,'FILE')
            tmp      = 'Warning step %s is not a standard timestep as ' %(str(parse().getValue(filenode)))
            print >> self.__fp
            print >> self.__fp, tmp
            print >> self.__fp, "cycle, time, dt, dtcourant, id are not provided."
            print >> self.__fp, "Will continue to postprocess.."
            print >> self.__fp
        try:
            filenode        = parse().getElement(node,'FILE')
            datafmt         = parse().getElemAttribute(filenode,'FORMAT')
            file            = str(parse().getValue(filenode))
            assert datafmt == 'binary'
            #given path of filename, ensure actual path of binary file is correct
            corrfile        = getpath(TBrookfilename,file) 
            file            = BrookUtilities().fixFileName(corrfile)
           
            Estring = 'import error'
            from POSTPROCESSORutils.dataSource import dataSource
            Estring = 'getting data source'
            tfile   = dataSource(file,'f95b',debug=0)
            if tfile.status > -1 :
                #look aside file exists..
                Estring         = 'setting pointer'
                self.__pointer    = {'node':node,'datafmt':datafmt,
                                   'tfile':tfile,'spacelu':spacelu}
                Estring         = 'getting elements'
                elems           = parse().getElements(node,'FILEVAR')
                #check for corruption of the binary lookaside file
                if len(elems) != len(tfile.getInfo(None,None,None,None)):
                    tfile.status = -1
                Estring          = 'looping'
                for element in elems:
                    thisvar      = getVariable(self.__pointer,element,parse)
                    self.vlist.append(thisvar)
                    
            if tfile.status < 0:
                print >> self.__fp, 'Caution data file associated with cycle',
                print >> self.__fp, '%i is corrupted' %(self.cycle)
                print >> self.__fp, 'or not available. Will not include this cycle in any'
                print >> self.__fp, 'postprocessing operations. '
                print >> self.__fp 
                self.cycle     = None
                self.time      = None
                self.dt        = None
                self.dtcourant = None
                self.vlist     = None
                return
        except:
            print >> self.__fp
            print >> self.__fp, 'No data exists in this cycle '
            print >> self.__fp, 'Will not include this cycle in any postprocessing operations'
            print >> self.__fp
            return

        "now place em data in the time step if it exists"
                
        try:
            self.emseq      = int(parse().getElemAttribute(emnode,'SEQ'))
        except:
            "no em data to be placed in this timestep"
        else:
            filenode        = parse().getElement(emnode,'FILE')
            datafmt         = parse().getElemAttribute(filenode,'FORMAT')
            file            = str(parse().getValue(filenode))
            self.emfile     = file
            assert datafmt == 'binary'
            #given path of filename, ensure actual path of binary file is correct
            corrfile        = getpath(TBrookfilename,file)
            file            = BrookUtilities().fixFileName(corrfile)
            tfile           = dataSource(file,'f95b',debug=0)            
            if tfile.status > -1:
                Estring         = 'setting pointer'
                self.__pointer    = {'node':emnode,'datafmt':datafmt,
                                   'tfile':tfile,'spacelu':spacelu}
                Estring         = 'getting elements'
                filelems        = parse().getElements(emnode,'FILEVAR')
                elems           = parse().getElements(emnode,'COIL')
                
                #check for corruption of the binary lookaside file
                if len(filelems) != len(tfile.getInfo(None,None,None,None)):
                    tfile.status = -1

                if tfile.status < 0:
                    print >> self.__fp, 'Caution em data file associated with sequence',
                    print >> self.__fp, '%i is corrupted' %(self.emseq)
                    print >> self.__fp, 'or not available. Will not be able to obtain data',
                    print >> self.__fp, 'for joule heat variables.'
                    print >> self.__fp, 'Will continue to postprocess but restart creation',
                    print >> self.__fp, 'will not be effective'
                    print >> self.__fp
                else:
                    Estring          = 'looping'
                    "consider coil elements"
                    newelems        = []
                    for elem in elems:
                        cname = parse().getElemAttribute(elem,'NAME')
                        vars  = parse().getElements(elem,'FILEVAR')
                        for var in vars:
                            thisvar       = getVariable(self.__pointer,var,parse)
                            thisvar.name  = cname + ' ' + thisvar.name
                            self.vlist.append(thisvar)
                            newelems.append(var)
                    "now consider remaining elements"
                    for filelem in filelems:
                        if not filelem in newelems:
                            thisvar = getVariable(self.__pointer,filelem,parse)
                            self.vlist.append(thisvar)


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
    opts = uTestOpts('fd')
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)
    
    try:
        parse       = MiniDom
        node        = parse().getDoc(opt.f)
        tnode       = parse().getElement(node,'TIMESTEP')
        
        spacelu     = {None:'CELL','NONE':'NONE',
                       'CELL':'CELL', 'CELLS':'CELL',
                       'FACES':'FACE', 'FACE':'FACE',
                       'VERTICES':'VERTEX', 'VERTEX':'VERTEX',
                       'NODES':'VERTEX','NODE':'VERTEX',
                       'EDGES':'EDGE','EDGE':'EDGE',
                       'MESH':'MESH'}

        #emnode      = parse().getElement(node,'JOULEHEAT')
        emnode      = None
        xmltimestep = getTimeStep(tnode,emnode,opt.f,spacelu)

        print 'timestep:'
        xmltimestep.metastr()

    except:
        print "---> Test failed in some aspect <---"
        print "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise
       

