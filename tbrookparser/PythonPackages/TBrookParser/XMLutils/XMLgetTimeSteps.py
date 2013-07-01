"""

 XMLgetTimeSteps

 -----------------------------------------------------------------------------
  Purpose:
  
     <purpose>
  
  Public Interface(s):
  
    ts = getTimeSteps(filename,spacelu,thedoc)
    ts.__getPreviousEMElement(element,cycle,emelems,parse)
    ts.__appendLastTimeStep(thedoc,filename,spacelu,parse)
    ts.__removeLastIncompleteTimeStep()
    ts.str()

  Contains:
    class getTimeSteps
        __init__(self,filename,spacelu,thedoc,parse=MiniDom)
        __getPreviousEMElement(self,element,cycle,emelems,parse)
        __appendLastTimeStep(self,thedoc,filename,spacelu,parse)
        __removeLastIncompleteTimeStep(self)
        str(self)

    Unit Test Block
  
  Author(s): Sharen Cummins (scummins@lanl.gov)
             Erin Iesulauro Barker (eibarker@lanl.gov)
 -----------------------------------------------------------------------------
"""

if __name__=='__main__':
    # Set sys.path for component testing mode
    import os, sys
    thisdir    = os.path.dirname(__file__)
    parserdir  = os.path.abspath(thisdir+'../')
    sys.path.append(parserdir)

import sys
from BASEutils      import baseTimeSteps
from XMLutils       import MiniDom, getTimeStep
from PYTHONutils    import unique, getpath

class getTimeSteps(baseTimeSteps):

    def __init__(self,filename,spacelu,thedoc,fp=sys.stdout,parse=MiniDom):

        "initialise TimeSteps object"
	baseTimeSteps.__init__(self)

        self.__pointer    = None
        self.__fp         = fp

        self.__elems      = parse().getElements(thedoc,'TIMESTEP')
        self.__emelems    = parse().getElements(thedoc,'JOULEHEAT')

        "create a list of timesteps "

        for element in self.__elems:
            id          = int(parse().getElemAttribute(element,'ID'))
            time        = float(parse().getElemAttribute(element,'t'))
            cycle       = int(parse().getElemAttribute(element,'CYCLE'))

            emelement   = self.__getPreviousEMElement(element,cycle,
                                                    self.__emelems,parse)
            
            timestep    = getTimeStep(element,emelement,filename,spacelu,self.__fp,parse)

            if len(self.__emelems) > 0 and timestep.cycle != None :
                print >> self.__fp
                print >> self.__fp, '***************************************************'
                print >> self.__fp, 'For cycle %i using em sequence %i' %(cycle,timestep.emseq)
                print >> self.__fp, 'and em file %s' %(timestep.emfile)
                print >> self.__fp, '***************************************************'
                print >> self.__fp

            if timestep.cycle != None:

                "cycle data exists so incorporate it into the storage object.."
                
                self.tlist.append(timestep)
                self.idlist.append(id)
                self.timeslist.append(time)
                self.cyclelist.append(cycle)

                "create a list of variables in each cycle"

                self.vlist[id] = []
                velems         = timestep.vlist

                for velem in velems:
                    self.vlist[id].append(velem)
                    "initialise the clist dictionary"
                    self.clist[velem.name] = []


        self.__appendLastTimeStep(thedoc,filename,spacelu,parse)

        "create a list of cycles for each variable -  i.e the clist dictionary"

        cnt = 0
        for id in self.idlist:
            thiscycle = self.cyclelist[cnt]
            for velem in self.vlist[id]:
                self.clist[velem.name].append(thiscycle)
            cnt += 1

        """
        Check if all variables are written out in this timestep. If not,
        remove this last timestep completely from the storage structure.
        """

        self.__removeLastIncompleteTimeStep()

    def __getPreviousEMElement(self,element,cycle,emelems,parse):

        "find relevant em data associated with this cycle"
        emelement   = None
            
        """
        Get latest preceding emelement to this timestep element.
        Do this by obtaining the parent CYCLE node to JOULE HEAT.
        Compare this cycle to the truchas timestep cycle.
        """
            
        count   = 0
        emcycle = 0
        while emcycle <= cycle :
            if count > 0:
                emparent = parse().getElemParent(emelems[count])
                emcycle  = int(parse().getElemAttribute(emparent,"ID"))
                if emcycle > cycle : break
            count   += 1
            if count > len(emelems) - 1 : break
        if len(emelems) == 1:
            emelement = emelems[0]
        if len(emelems) > 1:
            assert count > 0
            emelement = emelems[count-1] 

        return emelement


    def __appendLastTimeStep(self,thedoc,filename,spacelu,parse):

        """
        check if last successful timestep has been outputted -
        if so LASTSTEP flag exists append this to the timesteps object
        """
        
        elems               = parse().getElements(thedoc,'LASTSTEP')
        if len(elems) > 0:

            errstring  = '\n \n Multiple LASTSTEP flags exist in %s \n' %(filename)
            errstring += '      There should only be one! Aborting!! \n'

            assert len(elems) == 1, errstring
            
            filenode   = parse().getElement(elems[0],'FILE')
            filefmt    = parse().getElemAttribute(filenode,'FORMAT')

            errstring  = '\n \n LASTSTEP lookaside file should be xml format. \n' 
            errstring += '      Instead it is %s ! Aborting!! \n' %(filefmt)

            assert filefmt == 'xml', errstring

            laststpfile     = str(parse().getValue(filenode))
            corrlaststpfile = getpath(filename,laststpfile)
            node            = parse().getDoc(corrlaststpfile)
            tnode           = parse().getElement(node,'LASTSTEP')
            emnode = self.__getPreviousEMElement(self.__elems[-1],self.cyclelist[-1],self.__emelems,parse)
            "LJC Note: the syntax s[len(s)-1] is the same as s[-1] and is the last element"
            
            "now obtain the timestep object"
            timestep            = getTimeStep(tnode,emnode,
                                              corrlaststpfile,spacelu,self.__fp,parse)
            id                  = max(self.idlist) + 1
            timestep.id         = id
            time                = float(parse().getElemAttribute(tnode,'t'))
            cycle               = int(parse().getElemAttribute(tnode,'CYCLE'))

            """
            now append to our lists -
            if cycle is unique from standard timestep cycles
            """
            if cycle not in self.cyclelist:
                self.tlist.append(timestep)
                self.idlist.append(id)
                self.timeslist.append(time)
                self.cyclelist.append(cycle)

            self.vlist[id] = []
            velems         = timestep.vlist

            for velem in velems:
                self.vlist[id].append(velem)


    def __removeLastIncompleteTimeStep(self):

        if len(self.idlist) > 0:
            endindex = len(self.idlist)-1
            endid    = self.idlist[endindex]    

            """
            LJC Note: the syntax x[0:-1] or x[:-1]
            gives all but the last item/char in a list/string
            and is equivalent to x[0:len(x)-1] (and much shorter)
            """
            if len(self.tlist[endindex].vlist) < len(self.tlist[0].vlist):
                self.tlist     = self.tlist[0:-1]
                self.idlist    = self.idlist[0:-1]
                self.timeslist = self.timeslist[0:-1]
                self.cyclelist = self.cyclelist[0:-1]
                del self.vlist[endid]
                for velemname in self.clist:
                    self.clist[velemname] = self.clist[velemname][0:-1]

    """
    def str(self):
	
        for thistime in self.tlist:
            thistime.str()
    """
        
if __name__=='__main__':
    
    'testing this component'

    from PYTHONutils import uTestOpts
    opts = uTestOpts('fd')
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    fp   = sys.stdout
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

        xmltimesteps = getTimeSteps(opt.f,spacelu,node,fp,parse)
        print 'timesteps'
        xmltimesteps.metastr()

    except:
        print "---> Test failed in some aspect <---"
        print "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise


    



        

