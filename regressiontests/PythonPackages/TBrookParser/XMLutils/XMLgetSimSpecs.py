

"""

 XMLgetSimSpecs

 -----------------------------------------------------------------------------
  Purpose:
  
     To instantiate and provide methods for SimSpecs objects.  A
     SimSpecs object contains information about the simulation
     including such things as:
       'CODE', 'LIBRARIES', 'BUILDARCHITECTURE', 'BUILDDATETIME',
       'BUILDFLAGS', 'BUILDHOST', 'RUNARCHITECTURE', 'RUNHOST',
       'RUNDATE'
  
  Public Interface(s):
  
     s = getSimSpecs(filename)
     s.str()

  Contains:
     class getSimSpecs
         __init__(self,filename,parse=MiniDom,debug=0)
         __getThisSpec(self,theDoc,specName)
         __getThePhysics(self,theDoc,physicsName)
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

from PYTHONutils     import SystemExit
from XMLutils        import MiniDom
from BASEutils       import baseSimSpecs

class getSimSpecs(baseSimSpecs):

    def __init__(self,filename,parse=MiniDom,fp=sys.stdout,debug=0):

        "initialize the simspecs object"
        baseSimSpecs.__init__(self)

        self.__parse = parse
        self.__fp    = fp     
        
        try:
            theDoc         = parse().getDoc(filename,self.__fp)
            self.thedoc    = theDoc
        except:
            print >> self.__fp
            print >> self.__fp, 'An error has occurred in the parsing the TBrook.xml file. '
            print >> self.__fp, 'Please check the XML file is formed correctly.'
            SystemExit(self.__fp)
        
        variable           = parse().getElement(theDoc,'TruchasData')
        tmp                = parse().getElement(variable,'Variable')
        tmp                = parse().getElement(variable,'DataValues')
        self.tbrookvers    = float(parse().getValue(tmp))
        
        try:
            variable  = parse().getElement(theDoc,'PROGRAMSPECIFICATIONS')
            
            specslist = ['CODE','LIBRARIES','BUILDARCHITECTURE',
                         'BUILDDATETIME','BUILDFLAGS','BUILDHOST',
                         'RUNARCHITECTURE','RUNHOST','RUNDATE']

            for i in specslist:
                tmp          = self.__getThisSpec(variable,i)
                self.sspecs.append(tmp)
                self.nsspecs = self.nsspecs + 1

            """
            special treatment for the 'RUNPROCESSORS' simulation specification
            - convert to a string of length 256
            """
            tmp              = parse().getElement(variable,'RUNPROCESSORS')
            tmp              = parse().getValue(tmp)
            newtmp           = {}
            newtmp['value']  = str(tmp)
            newtmp['length'] = 256
            self.sspecs.append(newtmp)
            self.nsspecs     = self.nsspecs + 1

            variable         = parse().getElement(theDoc,'SIMULATIONINFORMATION')
            tmp              = self.__getThisSpec(variable,'TITLE')
            self.sspecs.append(tmp)
            self.nsspecs     = self.nsspecs + 1


            """
            Given the PHYSICS incorporated in the original
            simulation construct the #features to be recorded in the
            restart file.
            """
            physics        = {}
            physics        = self.__getThePhysics(variable,'PHYSICS')
            if physics.has_key('heat_conduction'):
                if physics['heat_conduction'] == 'T':
                    self.feats.append('temperature')
                    self.nfeats             = self.nfeats + 1
            if physics.has_key('fluid_flow'):
                if physics['fluid_flow'] == 'T':
                    self.feats.append('fluid_flow')
                    self.nfeats             = self.nfeats + 1
            if physics.has_key('electromagnetics'):
                if physics['electromagnetics'] == 'T':
                    self.feats.append('joule_heat')
                    self.nfeats             = self.nfeats + 1
            if physics.has_key('solid_mechanics'):
                if physics['solid_mechanics'] == 'T':
                    self.feats.append('solid_mechanics')
                    self.nfeats             = self.nfeats + 1

            alloyinfo      = parse().getElement(variable,'ALLOYINFORMATION')
            self.nphases   = int(parse().getElemAttribute(alloyinfo,'NPHASES'))
            self.ncomps    = int(parse().getElemAttribute(alloyinfo,'NCOMPONENTS'))

            if self.nphases > -1 or self.ncomps > -1:
                self.feats.append('alloy')
                self.nfeats += 1

	    try:
                speciesinfo    = parse().getElement(variable,'SPECIESINFO')
                self.nspecies  = int(parse().getElemAttribute(speciesinfo,'NSPECIES'))
	        if self.nspecies > 0:
		    self.feats.append('species')
		    self.nfeats += 1
            except:
	        "no species information"

            try:
                csf         = parse().getElement(variable,'MESHSCALEFACTOR')
                self.csf    = float(parse().getValue(csf).encode())
            except:
                self.csf    = 1.0

            try:
                sensitivity = parse().getElement(variable,'SENSITIVITY')
                self.feats.append('sensitivity')
		#EIB: Sensitivity feature wasn't be written out in restart header
		#     If including these variables causes trouble commented out
		#     the following line and sensitivity variables won't be read be
		#     restart_driver.F90
		self.nfeats = self.nfeats + 1
            except:
                "no sensitivity parameters"
                
        except:
            print >> self.__fp
            print >> self.__fp, 'Warning: No program specifications supplied in the XML file'
            print >> self.__fp, 'The parser will only be able to generate the mesh and',\
                  'abort variables.'
            print >> self.__fp

    def __getThisSpec(self,theDoc,specName):

        #returns a formatted string associated with a particular specification

        parse          = self.__parse
        node           = parse().getElement(theDoc,specName)
        strl           = int(parse().getElemAttribute(node,'TSTRLEN'))
        valu           = parse().getValue(node).encode()
        spec           = {}
        spec['value']  = valu
        spec['length'] = strl
        
        return spec

    def __getThePhysics(self,theDoc,physicsName):
        
        #Returns a string containing information about what PHYSICS 
        #have been incorporated in the original simulation. 
        
        parse     = self.__parse
        variables = parse().getElements(theDoc,physicsName)
        physics = {}
        for variable in variables:
            name  = parse().getElemAttribute(variable,'NAME').lower()
            value = parse().getElemAttribute(variable,'VALUE')
            physics[str(name)] = str(value)

        return physics
        
if __name__=='__main__':

    'test creating XML simspecs structure'

    from PYTHONutils import uTestOpts
    opts = uTestOpts('fd')
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    try:    
        xmlsimspecs = getSimSpecs(opt.f)
        xmlsimspecs.debug()
    except:
        print "---> Component test failed in some aspect <---"
        print "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise

