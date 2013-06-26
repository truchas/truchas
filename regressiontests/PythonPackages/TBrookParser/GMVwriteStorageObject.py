

"""

 GMVwriteStorageObject

 -----------------------------------------------------------------------------
  Purpose:

     Defines a subclass of writeStorageObject for the writing of GMV
     graphics files from Truchas simulation output.

  Public Interface(s):

    GwSO = GMVwriteStorageObject(fileName, outputfmt, mesh, t, cycle, bVariables)
    GwSO.bImport()
    GwSO.writebHeader()
    GwSO.writebMesh()
    GwSO.writebVariables()
    GwSO.writebFlags(flags)
    GwSO.writebGroups(flags)
    GwSO.bfinalizeAndClose()
    GwSO.writeHeader()
    GwSO.writeMesh()
    GwSO.writeFlags(flags)
    GwSO.writeGroups(flags)
    GwSO.writeCase()
    GwSO.writeVariables()
    GwSO.finalizeAndClose()

  Contains:
    class GMVwriteStorageObject(writeStorageObject)
        __init__(self, fileName, outputfmt, mesh, t, cycle, bVariables,
                 seq_no   = 0,
                 bFlags   = None,
                 bRegions = None,
                 dformat  = '%20.5e',
                 fpwatch  = sys.stdout,
                 debug    = 0)
        bImport(self)
        writebHeader(self)
        writebMesh(self)
        writebVariables(self)
        writebFlags(self, flags)
        writebGroups(self, flags)
        bfinalizeAndClose(self)
        writeHeader(self)
        writeMesh(self)
        writeFlags(self, flags)
        writeGroups(self, flags)
        writeCase(self)
        writeVariables(self)
        finalizeAndClose(self)

    Unit Test Block

  Version:
    $ID$

  Author(s): Sharen Cummins (scummins@lanl.gov)
             Erin Iesulauro Barker (eibarker@lanl.gov)
 ----------------------------------------------------------------------------
"""

from writeStorageObject import writeStorageObject
import os, sys, re
from struct import pack
try:
    import numpy.oldnumeric as Numeric
except ImportError:
    import Numeric
except:
    raise

class GMVwriteStorageObject(writeStorageObject):

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
        self.mesh     = mesh
        self.bFlags   = bFlags
        self.bRegions = bRegions
        self.debug    = debug
        self.seq_no   = 0
        self.mapLookup  = {None:0, 'CELL':0, 'NODE':1, 'VERTEX':1, 'FACE':2}
        self.typeLookup = {'i':'d', 'd':'f', 'l':'d', 'c':'s'}
        self.fpwatch    = fpwatch
        self.dformat    = dformat

        if (outputfmt == 'ascii'):
            self.FP     = open(fileName, 'w')
            self.asciiWrite()
        if (outputfmt == 'binary'):
            self.FP     = open(fileName, 'wb')
            self.binaryWrite()

        return

    #------------------------------------------------------------------------
    #define all the (b)inary routines specific to GMV
    #------------------------------------------------------------------------
    def bImport(self):
        # no-op
        pass

    #------------------------------------------------------------------------
    def writebHeader(self):

        self.FP.write("gmvinput")
        self.FP.write("ieeei4r8")

    #------------------------------------------------------------------------
    def writebMesh(self):

        mesh = self.mesh

        # Write out cell co-ords
        sys.stdout.flush()
        if (self.debug): print >> self.fpwatch, '  Coords'
        v      = mesh.vertices
        nnodes = v['N']
        vt     = Numeric.transpose(v['COORDS'])
        xc     = Numeric.array(vt[0])
        yc     = Numeric.array(vt[1])
        zc     = Numeric.array(vt[2])
        self.FP.write('nodes   ')
        self.FP.write(pack('i',nnodes))
        if self.debug:
            xcstr = xc.tostring()
            print >> self.fpwatch, "  x nodes string length for %d nodes: %d (should be %d)" \
                  %(nnodes,len(xcstr),8*nnodes)
        self.FP.write(xc.tostring())
        self.FP.write(yc.tostring())
        self.FP.write(zc.tostring())

        # Write out the cells vertices.
        sys.stdout.flush()
        if (self.debug): print >> self.fpwatch, '  Cell Vertices'

        c      = mesh.cells
        ncells = c['N']
        cVerts = c['VERTICES']
        self.FP.write('cells   ')
        self.FP.write(pack('i',ncells))

        if (mesh.type == 'HEX'):
            s  = 'hex\0\0\0\0\0'
            nv = 8
        elif (mesh.type == 'TET'):
            s  = 'tet\0\0\0\0\0'
            nv = 4

        for element in cVerts:
            cv = Numeric.array(element)
            self.FP.write(s)
            self.FP.write(pack('i',nv))
            self.FP.write(cv.tostring())

        #Write out material ids using mesh block information
        if mesh.cells.has_key('BLOCKID'):

            cblocks  = mesh.cells['BLOCKID']
            ncblocks = cblocks[Numeric.argmax(cblocks)]

            self.FP.write("material")
            self.FP.write(pack('2i',ncblocks,0))

            for i in range(ncblocks):
                mat_name = 'mat_%i' %(i+1)
                mat_name = '%-8s' %(mat_name[:min(8,len(mat_name))])
                self.FP.write(re.sub(' ','\0',mat_name))

            self.FP.write(cblocks.tostring())

        # Write variable header
        self.FP.write('variable')

        #write out centroid data - appropriate for solid mechanics visualisation
        if (self.debug): print >> self.fpwatch, '  Centroids'

        x  = mesh.cells['CENTROIDS']
        vt = Numeric.transpose(x)
        cx = Numeric.array(vt[0])
        self.FP.write('Centx\0\0\0')
        self.FP.write(pack('i',0))
        self.FP.write(cx.tostring())
        cy = Numeric.array(vt[1])
        self.FP.write('Centy\0\0\0')
        self.FP.write(pack('i',0))
        self.FP.write(cy.tostring())
        cz = Numeric.array(vt[2])
        self.FP.write('Centz\0\0\0')
        self.FP.write(pack('i',0))
        self.FP.write(cz.tostring())

        #write out number of neighbours data for meshspace = CELL
        if (self.debug): print >> self.fpwatch, '  Neighbors'

        x = mesh.cells['NUMNEIGHBORS']
        x = Numeric.array(x,'d')
        self.FP.write('Ngbrs\0\0\0')
        self.FP.write(pack('i',0))
        self.FP.write(x.tostring())

    #------------------------------------------------------------------------
    def writebVariables(self):

        # Write variable data
	print >> self.fpwatch, '  Variables:'
        for x in self.vars:
            if (self.debug):
                print >> self.fpwatch, '    ',x.name,x.rank,x.shape,x.meshspace,x.nickname,self.writeThisVar(x)
            if (x.rank == 1):
                #print out variables that live on a mesh only
                if (len(x.mesh)):
                    #print out variables belong to a predefined meshspace 
                    if (self.mapLookup.has_key(x.meshspace)) :
                        writevar = self.writeThisVar(x)
                        if writevar:
                            vname = x.nickname
                            dtype = self.mapLookup[x.meshspace]
                            vname = '%-8s' %(vname[:min(8,len(vname))])
                            self.FP.write(re.sub(' ','\0',vname))
                            self.FP.write(pack('i',dtype))
                            self.FP.write(Numeric.array(x.data).tostring())

            # write out rank=2 data
            elif (x.rank == 2):
                writevar = self.writeThisVar(x)
                if writevar:
                    #provide a generic name for each component except for velocity
                    comps        = []
                    if (x.nickname != 'Velocity'):
                        for i in range(x.shape[1]):
                            dtype    = self.mapLookup[x.meshspace]
                            thiscomp = '%s_%s %i' %(x.nickname,i+1,dtype)
                            comps.append(thiscomp)
                    if (x.name == 'sigma'):
                        comps     = ['sigxx 0 \n', 'sigyy 0 \n', 'sigzz 0 \n',
                                      'sigxy 0 \n', 'sigxz 0 \n', 'sigyz 0 \n']
                    if (x.name == 'epsilon'):
                        comps     = ['epsxx 0 \n', 'epsyy 0 \n', 'epszz 0 \n',
                                      'epsxy 0 \n', 'epsxz 0 \n', 'epsyz 0 \n']
                    if (x.name == 'e_plastic'):
                        comps     = ['eplxx 0 \n', 'eplyy 0 \n', 'eplzz 0 \n',
                                      'eplxy 0 \n', 'eplxz 0 \n', 'eplyz 0 \n']
                    if (x.name == 'epstherm'):
                        comps     = ['epsthxx 0 \n', 'epsthyy 0 \n', 'epsthzz 0 \n',
                                     'epsthxy 0 \n', 'epsthxz 0 \n', 'epsthyz 0 \n']
                    if (x.name == 'epspc'):
                        comps     = ['epspcxx 0 \n', 'epspcyy 0 \n', 'epspczz 0 \n',
                                     'epspcxy 0 \n', 'epspcxz 0 \n', 'epspcyz 0 \n']
                    if (x.name == 'Displacement'):
                        comps     = ['Dx 1 \n', 'Dy 1 \n', 'Dz 1 \n']
                    vt        = Numeric.transpose(x.data)
                    dtype     = self.mapLookup[x.meshspace]
                    for comp in comps:
                        ind   = comps.index(comp)
                        c     = Numeric.array(vt[ind])
                        scomp = "%-8s" %(comp[:min(8,len(comp))])
                        scomp = re.sub(' ','\0',scomp[:8])
                        self.FP.write(scomp)
                        self.FP.write(pack('i',dtype))
                        self.FP.write(c.tostring())

        #write variable tail
        self.FP.write('endvars ')

        # Look for, and write velocities
        for x in self.vars:
            if (x.rank == 2 and x.name == 'Z_VC' and x.nickname == 'Velocity'):
                dtype = self.mapLookup[x.meshspace]
                vt    = Numeric.transpose(x.data)
                xv    = Numeric.array(vt[0])
                yv    = Numeric.array(vt[1])
                zv    = Numeric.array(vt[2])
                self.FP.write('velocity')
                self.FP.write(pack('i',dtype))
                self.FP.write(xv.tostring())
                self.FP.write(yv.tostring())
                self.FP.write(zv.tostring())

    #------------------------------------------------------------------------
    def writebFlags(self, flags):

        self.FP.write('flags\0\0\0')

        if (self.mesh.cells.has_key('PARTITIONS')):

            cpart = self.mesh.cells['PARTITIONS']
            npe   = cpart[Numeric.argmax(cpart)]

            self.FP.write("cellpart")
            self.FP.write(pack('2i',npe,0))
            
            for i in range(npe):
                flag_name = 'P%i' %(i+1)
                flag_name = '%-8s' %(flag_name[:min(8,len(flag_name))])
                self.FP.write(re.sub(' ','\0',flag_name))

            self.FP.write(cpart.tostring())

        if (self.mesh.vertices.has_key('PARTITIONS')):

            vpart = self.mesh.vertices['PARTITIONS']
            npe   = vpart[Numeric.argmax(vpart)]

            self.FP.write("vertpart")
            self.FP.write(pack('2i',npe,1))
            
            for i in range(npe):
                flag_name = 'P%i' %(i+1)
                flag_name = '%-8s' %(flag_name[:min(8,len(flag_name))])
                self.FP.write(re.sub(' ','\0',flag_name))

            self.FP.write(vpart.tostring())            

        try:
            n = self.mesh.cells['N']
            for f in flags.keys():
                values  = Numeric.array([1],'i')
                values  = Numeric.resize(values,[n])
                if (flags[f] != None):
                    if (len(flags[f])):
                        for x in flags[f]:
                            values[x-1] = 2
                nf    = 2
                dtype = 0
                f     = '%-8s' %(f[:min(8,len(f))])
                self.FP.write(re.sub(' ','\0',f))
                self.FP.write(pack('2i', 2, dtype))
                self.FP.write('off\0\0\0\0\0')
                self.FP.write('on\0\0\0\0\0\0')
                self.FP.write(values.tostring())
        except:
            "continue..."

        self.FP.write('endflag\0')

    #------------------------------------------------------------------------
    def writebGroups(self, flags):

        if (flags == None):
            return
        elif (len(flags) == 0):
            return
        self.FP.write('groups\0\0')
        for f in flags.keys():
            if (flags[f] != None):
                if (len(flags[f])):
                    values = Numeric.array(flags[f],'i')
                    ne     = len(flags[f])
                    dtype  = 0
                    f = '%-8s' %(f[:min(8,len(f))])
                    self.FP.write(re.sub(' ','\0',f))
                    self.FP.write(pack('2i',dtype,ne))
                    self.FP.write(values.tostring())

        self.FP.write('endgrp\0\0')

    #------------------------------------------------------------------------
    def bfinalizeAndClose(self):

        if (self.debug): print >> self.fpwatch, '  Closing'
        self.FP.write('probtime')
        self.FP.write(pack('d',self.t))
        self.FP.write('endgmv  ')
        self.FP.close()

    #------------------------------------------------------------------------
    # define the ASCII routines
    #------------------------------------------------------------------------
    def writeHeader(self):
        #for the gmv ascii writer just write a simple header

        print >> self.FP, 'gmvinput ascii'

    #------------------------------------------------------------------------
    def writeMesh(self):

        #write the mesh structure

        fp = self.FP

        if (self.debug): print >> self.fpwatch, '  Vertices',
        sys.stdout.flush()

        v = self.mesh.vertices
        print >> fp, 'nodes', v['N']

        if (self.debug): print >> self.fpwatch, '  Coords',
        sys.stdout.flush()

        self.writeArray(1,v['COORDS'][:,0],'',self.dformat,fp)
        print >> fp
        self.writeArray(1,v['COORDS'][:,1],'',self.dformat,fp)
        print >> fp
        self.writeArray(1,v['COORDS'][:,2],'',self.dformat,fp)
        print >> fp

        # Write out the cells vertices.  Need total number of cells, and also
        # total number of vertices for all cells.

        if (self.debug): print >> self.fpwatch, '  Cell Vertices',
        sys.stdout.flush()

        c      = self.mesh.cells
        cVerts = c['VERTICES']
        totalLength = 0
        print >> fp, 'cells', c['N']

        hd = 'hex'
        nv = str(8)
        if (self.mesh.type == 'HEX'):
            hd = 'hex'
            nv = 8
        elif (self.mesh.type == 'TET'):
            hd = 'tet'
            nv = 4

        hdr = hd + ' ' + str(nv)

        self.writeArray(2,cVerts, hdr, ' %7i', fp, nperline=8)

        #Write out material ids using mesh block information
        if self.mesh.cells.has_key('BLOCKID'):

            cblocks  = self.mesh.cells['BLOCKID']
            ncblocks = cblocks[Numeric.argmax(cblocks)]

            print >> fp,  'material %i %i' %(ncblocks,0)

            for i in range(ncblocks):
                mat_name = 'mat_%i' %(i+1)
                print >> fp, mat_name

            self.writeArray(1,cblocks,'',self.iformat,fp,nperline=8)

        print >> fp, 'variables'

        #write out centroid data - appropriate for solid mechanics visualisation

        if (self.debug): print '  Centroids'
        x          = self.mesh.cells['CENTROIDS']
        vt         = Numeric.transpose(x)
        comps      = ['Centx 0 \n', 'Centy 0 \n', 'Centz 0 \n']
        self.writeCompsArray(vt,comps,self.dformat,fp)

        #write out number of neighbours data for meshspace = CELL

        if (self.debug): print >> self.fpwatch, '  Neighbors'
        x = self.mesh.cells['NUMNEIGHBORS']
        self.writeArray(1,x, 'Ngbrs 0 \n', self.iformat, fp, nperline=10)
        print >> fp

        return

    #------------------------------------------------------------------------
    def writeFlags(self, flags):

        fp = self.FP
        print >> fp, 'flags'

        if (self.mesh.cells.has_key('PARTITIONS')):
            cpart = self.mesh.cells['PARTITIONS']
            npe   = cpart[Numeric.argmax(cpart)]
            f     = 'cellpartn %d %d' %(npe,0)
            print >> fp, f
            for i in range(npe):
                flag_name = 'P%i' %(i+1)
                f         = '%-8s' %(flag_name[:min(8,len(flag_name))])
                print >> fp, f
            s = ''
            for x in cpart:
                s += '%1d '%x
                if(len(s)>80):
                    print >> fp, s
                    s = ''
            print >> fp, s

        if (self.mesh.vertices.has_key('PARTITIONS')):
            vpart = self.mesh.vertices['PARTITIONS']
            npe   = vpart[Numeric.argmax(vpart)]
            f     = 'vellpartn %d %d' %(npe,1)
            print >> fp, f
            for i in range(npe):
                flag_name = 'P%i' %(i+1)
                f         = '%-8s' %(flag_name[:min(8,len(flag_name))])
                print >> fp, f
            s = ''
            for x in vpart:
                s += '%1d '%x
                if(len(s)>80):
                    print >> fp, s
                    s = ''
            print >> fp, s

        try:
            n = self.mesh.cells['N']
            for f in flags.keys():
                values  = Numeric.array([1],'i')
                values  = Numeric.resize(values,[n])
                if(flags[f]!=None):
                    if(len(flags[f])):
                        for x in flags[f]:
                            values[x-1]=2

                print >>fp, f, 2, 0
                print >>fp, 'off on'
                s = ''
                for x in values:
                    s += '%1d '%x
                    if(len(s)>80):
                        print >> fp, s
                        s = ''
                print >> fp, s
                #print >> fp
        except:
            "continue...."
            
        print >> fp, 'endflag'

    #------------------------------------------------------------------------
    def writeGroups(self, flags):

        if(flags == None) :
            return
        if(len(flags) == 0) :
           return

        fp = self.FP

        print >> fp, 'groups'
        n = self.mesh.cells['N']

        for f in flags.keys():
            values = []
            if(flags[f]!=None):
                if(len(flags[f])):
                    values = Numeric.array(flags[f],'i')

            print >>fp, f, 0, len(values)
            s = ''
            for x in values:
                s = s +  '%d '%(x)
                if(len(s)>80):
                    print >> fp, s
                    s = ''
            print >>fp, s
            #print >> fp
        print >>fp, 'endgrp'

    #------------------------------------------------------------------------
    def writeCase(self):

        pass

    #------------------------------------------------------------------------
    def writeVariables(self):

        fp     = self.FP

	print >> self.fpwatch, '  Variables:'
        for x in self.vars:
            if (self.debug):
                print >> self.fpwatch, '    ',x.name,x.rank,x.shape,x.meshspace,x.nickname,self.writeThisVar(x)
		#print >> self.fpwatch, x.data
            if(x.rank == 1):
                #print out variables that live on a mesh only
                if (len(x.mesh)):
                    #print out variables belong to a predefined meshspace
                    if (self.mapLookup.has_key(x.meshspace)) :
                        writevar = self.writeThisVar(x)
                        if writevar:
                            hdr = "%s %d \n" %(x.nickname, self.mapLookup[x.meshspace])
                            self.writeArray(x.rank,x.data,hdr,self.dformat,fp)
                            print >> fp

            if (x.rank == 2):
                #print solid mechanics quantities
                writevar = self.writeThisVar(x)
                if writevar:
                    comps         = []
                    #provide a generic name for each component except for velocity
                    if (x.nickname != 'Velocity'):
                        for i in range(x.shape[1]):
                            dtype     = self.mapLookup[x.meshspace]
                            thiscomp  = '%s_comp%s %i \n' %(x.nickname,i+1,dtype)
                            comps.append(thiscomp)
                    if (x.name == 'sigma'):
                        comps     = ['sigxx 0 \n', 'sigyy 0 \n', 'sigzz 0 \n',
                                      'sigxy 0 \n', 'sigxz 0 \n', 'sigyz 0 \n']
                    if (x.name == 'epsilon'):
                        comps     = ['epsxx 0 \n', 'epsyy 0 \n', 'epszz 0 \n',
                                      'epsxy 0 \n', 'epsxz 0 \n', 'epsyz 0 \n']
                    if (x.name == 'e_plastic'):
                        comps     = ['eplxx 0 \n', 'eplyy 0 \n', 'eplzz 0 \n',
                                      'eplxy 0 \n', 'eplxz 0 \n', 'eplyz 0 \n']
                    if (x.name == 'epstherm'):
                        comps     = ['epsthxx 0 \n', 'epsthyy 0 \n', 'epsthzz 0 \n',
                                     'epsthxy 0 \n', 'epsthxz 0 \n', 'epsthyz 0 \n']
                    if (x.name == 'epspc'):
                        comps     = ['epspcxx 0 \n', 'epspcyy 0 \n', 'epspczz 0 \n',
                                     'epspcxy 0 \n', 'epspcxz 0 \n', 'epspcyz 0 \n']
                    if (x.name == 'Displacement'):
                        comps     = ['Dx 1 \n', 'Dy 1 \n', 'Dz 1 \n']
                    vt         = Numeric.transpose(x.data)
                    self.writeCompsArray(vt,comps,self.dformat,fp)
        print >> fp, 'endvars'

        for x in self.vars:
            if(x.rank == 2 and x.name == 'Z_VC' and x.nickname == 'Velocity'):
                print >> fp, 'velocity', 0
                self.writeArray(1,x.data[:,0],'',self.dformat,fp)
                print >> fp
                self.writeArray(1,x.data[:,1],'',self.dformat,fp)
                print >> fp
                self.writeArray(1,x.data[:,2],'',self.dformat,fp)
                print >> fp
        return

    #------------------------------------------------------------------------
    def finalizeAndClose(self):

        s = self.dformat%self.t
        print >> self.FP, 'probtime ' + s
        print >> self.FP, 'endgmv'

# end of classGMVwriteStorageObject
#------------------------------------------------------------------------

if __name__ =='__main__':

    'component testing GMV write component'

    dfltdir = '../../scripts/test_TBrookParse/samples/'

    from PYTHONutils import uTestOpts

    # file, debug, output prefix, binary, ascii, clean
    opts = uTestOpts('fdoc',
                     defaults = {'c' : False,
                                 'o' : 'GMVwriteUnitTest'},
                     actions  = {'c' : 'store_true'},
                     dir      = dfltdir)
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    fpwatch    = sys.stdout
    try:
        if opt.c:
            print >> fpwatch, "Cleaning up from a past test run...",
            cmd = 'rm -f %s' %(dfltdir + opt.o + '.ascii')
            os.system(cmd)
            cmd = 'rm -f %s' %(dfltdir + opt.o + '.binary')
            os.system(cmd)
            print >> fpwatch, "Cleanup successful!!"
        else:

            from XMLgetStorageObject import getStorageObject

            storage   = getStorageObject(opt.f)
            mesh      = storage.mlist[0] #assume default mesh

            " first get all data values associated with this mesh"
            for meshspace in mesh.mslist:
                storage.getValues(meshspace.vlist)


            "fill the mesh with is actual coorindate values"
            mesh.fillMesh()


            itimeout  = 0
            t         = storage.tlist[itimeout].time
            vars      = storage.tlist[itimeout].vlist
            storage.getValues(vars)

            # ASCII test write
            print >> fpwatch, '*'*60
            print >> fpwatch, 'Testing ascii write...'
            fileout  = dfltdir + opt.o + '.ascii'
            writer   = GMVwriteStorageObject(fileout, 'ascii', mesh, t,
                                             itimeout, vars, fpwatch = fpwatch,
					     debug = opt.d)
            print >> fpwatch, 'GMV file %s\n  written to %s' % \
			      (opt.o+'.ascii',dfltdir)

            # Binary test write
            print >> fpwatch,'*'*60
            print >> fpwatch,'Testing binary write...'
            fileout  = dfltdir+ opt.o + '.binary'
            writer   = GMVwriteStorageObject(fileout, 'binary', mesh, t,
                                             itimeout, vars, fpwatch=fpwatch,
					     debug=opt.d)
            print >> fpwatch,'GMV file %s\n  written to %s' % \
			     (opt.o+'.binary',dfltdir)

            print >> fpwatch,'*'*60
            print >> fpwatch, \
	        '--> Note: to clean up, re-execute same command line plus "-c"'

    except:
        print >> fpwatch,"---> Test failed in some aspect <---"
        print >> fpwatch,"\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise

