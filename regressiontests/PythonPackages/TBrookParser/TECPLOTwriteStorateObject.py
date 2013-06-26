

"""

 TECPLOTwriteStorateObject

 -----------------------------------------------------------------------------
  Purpose:

     Defines a subclass of writeStorageObject for the writing of TECPLOT files from
     Truchas simulation files.

     Handles writing of ascii TECPLOT files from a given
     storage object and other information (mesh, timestep...).

     Writing of binary files is not supported.

  Public Interface(s):

    TECwSO = TECPLOTwriteStorageObject(fileName, outputfmt, mesh, t, cycle, bVariables)
    TECwSO.asciiWrite()
    TECwSO.writeHeader()
    TECwSO.writeVariables()
    TECwSO.finalizeAndClose()

  Contains:
    class TECPLOTwriteStorageObject(writeStorageObject)
        __init__(self, fileName, outputfmt, mesh, t, cycle, bVariables,
                 seq_no   = 0,
                 bFlags   = None,
                 bRegions = None,
                 dformat  = '%20.5e',
                 fpwatch  = sys.stdout,
                 debug    = 0)
        asciiWrite(self)
        writeHeader(self)
        writeVariables(self)
        finalizeAndClose(self)

    Unit Test Block

  Version:
    $ID$

  Author(s): Sharen Cummins (scummins@lanl.gov)
 -----------------------------------------------------------------------------
"""

#
# TECPLOT writer
# Written by Sriram Swaminarayan August 30, 2005

from writeStorageObject import writeStorageObject
import os, sys

class TECPLOTwriteStorageObject(writeStorageObject):

    def __init__(self, fileName, outputfmt, mesh, t, cycle, bVariables,
                 seq_no   = 0,
                 bFlags   = None,
                 bRegions = None,
                 dformat  = '%20.5e',
                 fpwatch  = sys.stdout,
                 debug    = 0):

        writeStorageObject.__init__(self,fpwatch,debug)
	self.fileName = fileName
        self.vars     = bVariables
        self.t        = t
        self.cycle    = cycle
        self.mesh     = mesh
	self.bFlags   = bFlags
	self.bRegions = bRegions
        self.debug    = debug
        self.seq_no   = seq_no
        self.mapLookup  = {None:0, 'CELL':0, 'NODE':1, 'VERTEX':1, 'FACE':2}
        self.typeLookup = {'i':'d', 'd':'f', 'l':'d', 'c':'s'}
        self.typeLookup2= {'i':'LONGINT', 'f':'SINGLE', 'd':'DOUBLE', 'l':'BYTE'}

        # special variables only for us
        self.cellVars = []
        self.nodeVars = []

        self.fpwatch  = fpwatch
        self.dformat  = dformat
	self.iformat  = '%7i'

        if (outputfmt == 'ascii'):
            self.asciiWrite()
        else:
            print "\nSorry, only tecplot ascii write is implemented.\n"

	return


    def asciiWrite(self):
        # Overwrite the ascii write function
        sys.stdout.flush()

        self.module_name = 'asciiwriter'
        if(self.seq_no==0):
            self.FP = open(self.fileName, "w")
            print >> self.fpwatch
            print >> self.fpwatch,'Writing ascii '+ self.fileName,
            print >> self.fpwatch
        else:
            print >> self.fpwatch
            print >> self.fpwatch,'  Appending sequence %d to ascii '%(self.seq_no)+ self.fileName,
            print >> self.fpwatch
            self.FP = open(self.fileName, "a")

        # Write the header
        self.writeHeader()

	# Write a case file
	self.writeCase()

        # Write the variables
        self.writeVariables()

        # Write the flags
        self.writeFlags(self.bFlags)
        self.writeGroups(self.bFlags)

        # Close the file
        self.finalizeAndClose()

        # return success
        return



    #define all the (b)inary routines specific to TECPLOT
    def writeHeader(self):
        try:
           import numpy.oldnumeric as Numeric
        except ImportError:
           import Numeric
        except:
           raise
	#for the tecplot ascii writer just write a simple header


        if (self.debug): print >> self.fpwatch, 'Header',
        sys.stdout.flush()

        fp = self.FP;
        if(self.seq_no==0):
            tmp = self.dformat%(self.t)
            fp.writelines([
                '# '+self.fileName+\
                ' written by Truchas OutputParser, time = ' + tmp + '\n\n', #%13.6f\n\n'%self.t,
                'Title="Truchas :'+self.fileName[self.fileName.rfind('/')+1:]+'"\n\n',
                ])

            fp.write('DATAPACKING=BLOCK\n')


            dtVars = ''
        else:
            dtVars = 'VARSHARELIST=([1,2,3]=1)\nCONNECTIVITYSHAREZONE=1\n'
            dtVars = 'VARSHARELIST=([1,2,3]=1)\n'

        cellVars = []
        nodeVars = []
        dtVars  += 'DT=(DOUBLE,DOUBLE,DOUBLE'


        for x in self.vars:
            # only rank 1 and above will be written
            writevar = self.writeThisVar(x) and x.rank > 0 and x.rank < 3
             
            """ 
            if(x.rank < 1 ): continue
            elif(x.name.startswith('M_')             or
                 x.name.startswith('Z_VF')           or
                 x.name.startswith('RHS')            or
                 x.name.startswith('TOTAL_STRAIN')   or
                 x.name.startswith('PLASTIC_STRAIN') or
                 x.name.startswith('ELASTIC_STRAIN') or
                 x.name.startswith('ELASTIC_STRESS') or
                 x.name.startswith('PLASTIC_STRESS') or
                 x.name.startswith('TOTAL_STRESS')
                 ): continue
            elif(x.name.startswith('Z_VF')): continue
            elif(x.rank>2):
                print >> self.fpwatch, 'Ignoring '+ x.name + ': Rank 3 and above unimplemented in tecplot'
                continue
            """

            if writevar:
                if self.mapLookup.has_key(x.meshspace) and self.typeLookup2.has_key(x.type):
                    myType = self.mapLookup[x.meshspace]
                    if myType == 0:
                        cellVars.append(x)
                    elif myType == 1:
                        nodeVars.append(x)
                    else:
                        print >> self.fpwatch, 'Variable: %s ignored because it has unknown map: %s'\
                              %(x.nickname,x.meshspace)
                        continue
                    dtVars += ', %s'%self.typeLookup2[x.type]

                    #Append rank 2  solid mechanics quantities
                    if (x.rank > 1):
                        for n in range(x.shape[1]-1):
                            dtVars += ', %s'%self.typeLookup2[x.type]

            else:
                """
                print >> self.fpwatch, 'Variable: %s ignored because it has unknown type: %s'\
                      %(x.nickname,x.meshspace)
                """
                continue

        dtVars += ')\n'

        nnVars = 3  # x, y,and z always present

        varlist = 'VARIABLES = "x" "y" "z"'

        for x in nodeVars:
            if x.rank == 1:
                varlist += ' "%s"'%x.nickname
                nnVars  += 1
            else:
                for n in range(x.shape[1]):
                    varlist += ' "%s_%d"'%(x.nickname,n+1)
                    nnVars  += 1


        ncVars = nnVars
        for x in cellVars:
            if x.rank == 1:
                varlist += ' "%s"'%x.nickname
                ncVars  += 1
            else:
                for n in range(x.shape[1]):
                    varlist += ' "%s_%d"'%(x.nickname,n+1)
                    ncVars  += 1


        varlist=varlist+"\n"
        if(self.seq_no == 0): fp.write(varlist)

        self.cellVars = cellVars
        self.nodeVars = nodeVars

        cellVars=None
        nodeVars=None


        v = self.mesh.vertices
        c = self.mesh.cells


        hd = 'FEBRICK'
        if (self.mesh.type == 'HEX'):
            hd = 'FEBRICK'
        elif (self.mesh.type == 'TET'):
            hd = 'FETETRAHEDRON'


        fp.write("""
        ZONE
        T="%s, time = %13.6f"
        ZONETYPE=%s
        N=%d
        E=%d
        DATAPACKING=BLOCK
        """%(self.mesh.name,self.t,hd,v['N'],c['N']))

        tmp = self.dformat%(self.t)
        fp.write('AUXDATA TIME=' + '"' + tmp + '"' + '\n') #"%-13.6f"\n'%(self.t))
        fp.write('AUXDATA CYCLE="%d"\n'%(self.cycle))
        fp.write('\n')

        varlocs = 'VARLOCATION=([1-%d]=NODAL'%nnVars
        if ( ncVars != nnVars):
            if(ncVars-nnVars>1):
                varlocs += ',[%d-%d]=CELLCENTERED'%(nnVars+1,ncVars)
            else:
                varlocs += ',%d=CELLCENTERED'%ncVars
        varlocs += ')\n'

        fp.write(varlocs)
        fp.write(dtVars)

        # write the x, y, and z arrays
        if(self.seq_no==0):
            vt     = Numeric.transpose(v['COORDS'])
            self.writeArray(2,vt,'',self.dformat,fp)

        c  = None
        v  = None
        vt = None

	return

    def writeVariables(self):

        try:
           import numpy.oldnumeric as Numeric
        except ImportError:
           import Numeric
        except:
           raise

        # Write out the cells vertices.  Need total number of cells, and also
        # total number of vertices for all cells.
        fp     = self.FP

        for x in self.nodeVars:
            tlu = self.typeLookup2[x.type]
            if(tlu == 'LONGINT' or tlu == 'BYTE' or tlu == 'SHORTINT'):
                fmt = self.iformat
            else:
                fmt = self.dformat
            if(x.rank == 1):
                vt = x.data
            else:
                vt = Numeric.transpose(x.data)
            self.writeArray(x.rank,vt,'',fmt,fp,nperline=10)

        for x in self.cellVars:
            tlu = self.typeLookup2[x.type]
            if(tlu == 'LONGINT' or tlu == 'BYTE' or tlu == 'SHORTINT'):
                fmt = self.iformat
            else:
                fmt = self.dformat
            if(x.rank == 1):
                vt = x.data
                self.writeArray(x.rank,vt,'',fmt,fp,nperline=10)
            else:
                vt    = Numeric.transpose(x.data)
                comps = []
                for i in vt:
                    comps.append('')
                self.writeCompsArray(vt,comps,fmt,fp,nperline=10)
                
        return

    def finalizeAndClose(self):
        try:
           import numpy.oldnumeric as Numeric
        except ImportError:
           import Numeric
        except:
           raise
        # Note that connectivity goes here
        fp = self.FP
        c = self.mesh.cells
        n = 8
        if(self.mesh.type == 'TET'):
            n = 4
        if(self.debug): print 'nnodes =', c['N']

        self.writeArray(2,c['VERTICES'],'\n', self.iformat,fp,nperline=8)
        textposn = 36 - self.seq_no*4
        textposn = 4
        fp.write('TEXT X=1, Y=%i, T="t=%-13.6f cycle=%d", F=TIMES, H=3, ZN=%d\n'\
                 %(textposn,self.t,self.cycle,self.seq_no+1))
        fp.write('\n\n')

        self.FP.close()
        self.FP = None

if __name__ =='__main__':

    'for component testing TECPLOT write component'
    dfltdir  = '../../scripts/test_TBrookParse/samples/'

    from PYTHONutils import uTestOpts

    # file, debug, output prefix, clean
    opts = uTestOpts('fdoc',
                     defaults = {'c' : False,
                                 'o' : 'TECPLOTwriteUnitTest'},
                     actions  = {'c' : 'store_true'},
                     dir      = dfltdir)
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    fpwatch = sys.stdout
    
    try:

        if opt.c:
            print >> fpwatch, "Cleaning up from a past test run...",
            cmd = 'rm -f %s' %(dfltdir + opt.o + '.ascii')
            os.system(cmd)
            print >> fpwatch, "Cleanup successful!!"

        else:
            from XMLgetStorageObject import getStorageObject

            "create a storage object and obtain its mesh"
            storage   = getStorageObject(opt.f)	
            mesh      = storage.mlist[0] #assume default mesh

            "get the data values for the mesh fields"
            for meshspace in mesh.mslist:
                meshspace.vlist = storage.getValues(meshspace.vlist)

            "fill in actual data values"
            mesh.fillMesh()

            print >> fpwatch, "All %d timesteps are written to a single file."%\
			      (len(storage.tlist))

	    # ASCII test
            itimeout  = 0
            for x in storage.tlist:
                t         = storage.tlist[itimeout].time
                vars      = storage.tlist[itimeout].vlist
                storage.getValues(vars)

                fileout  = dfltdir + opt.o + '.ascii'
                writer   = TECPLOTwriteStorageObject(fileout, 'ascii', mesh, t,
                                                     itimeout, vars,
                                                     seq_no=itimeout, 
						     fpwatch=fpwatch, 
						     debug=opt.d)
                itimeout = itimeout + 1

	    # Binary test not avaiable
	    #    If a test is added, add binary to clean up
	    print '\nWriting Binary ...  NOT TESTING BINARY, NOT AVAILABLE\n'

            print >> fpwatch, '*'*60
            print >> fpwatch, \
		'--> Note: to clean up, re-execute same command line plus "-c"'

    except:
        print >> fpwatch, "---> Test failed in some aspect <---"
        print >> fpwatch, "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise
