

"""

 RESTARTwriteStorageObject

 -----------------------------------------------------------------------------
  Purpose:

     Defines a subclass of writeStorageObject for the writing of Truchas
     restart files.
     
     Handles writing of binary and ascii restart files from a given
     storage object and other information (mesh, timestep...).

  Public Interface(s):

    RwSO = RESTARTwriteStorageObject(fileName, outputfmt, storage, mesh,
                                     itimeStep, t, bVariables)
    RwSO.writebHeader()
    RwSO.writebMesh()
    RwSO.writebVariables()
    RwSO.writeHeader()
    RwSO.writeMesh()
    RwSO.writeVariables()

  Contains:
    class RESTARTwriteStorageObject(writeStorageObject)
        __init__(self, fileName, outputfmt, storage, mesh,
                 itimeStep, t, bVariables,
                 bFlags   = None,
                 bRegions = None,
                 dformat  = '%20.5e',
                 debug    = 0)
        writebHeader(self)
        writebMesh(self)
        writebVariables(self)
        writeHeader(self)
        writeMesh(self)
        writeVariables(self)

    Unit Test Block

  Version:
    $ID$

  Author(s): Sharen Cummins (scummins@lanl.gov)
 -----------------------------------------------------------------------------
"""

from writeStorageObject import writeStorageObject
from FORTRANutils       import fortransupport	
import os, sys

class RESTARTwriteStorageObject(writeStorageObject):

    def __init__(self, fileName, outputfmt, storage, mesh,
                 itimeStep, t, bVariables,
                 bFlags   = None,
                 bRegions = None,
                 dformat  = '%20.5e',
                 fpwatch  = sys.stdout,
                 debug    = 0):

        writeStorageObject.__init__(self,fpwatch,debug)
        self.storage   = storage
        self.vars      = bVariables
        self.t         = t
        self.itimeStep = itimeStep
        self.mesh      = mesh
        self.header    = None
        self.fileName  = fileName
        self.seq_no    = 0
        self.debug     = debug
        self.bFlags    = bFlags
        sys.stdout.flush()	
        byteSwap       = 0
        self.fpwatch   = fpwatch
        self.dformat   = dformat

        if (outputfmt == 'binaryswap'):
            byteSwap = 1
            
        if (outputfmt == 'ascii'):
            self.FP = open(fileName,'w')
            self.asciiWrite()
        else:
            self.FP        = fortransupport.F95Binary()
            self.FP.open(fileName,'w',byteSwap=byteSwap)
            self.binaryWrite()

    #define all the (b)inary routines specific to RESTART

    def writebHeader(self):

        #Write restart header information to file

        strsize8    = 8
        strsize32   = 32
        strsize128  = 128
        strsize256  = 256
        strsize1024 = 1024

        specs       = self.storage.specs

        S = self.blankfill(specs.fileid,strsize8)
        self.setNandWrite(S,1)

        self.setNandWrite(specs.nfeats,1)

        for i in range(specs.nfeats):
            S = self.blankfill(specs.feats[i],strsize32)
            self.setNandWrite(S,1)

        self.setNandWrite(specs.nsspecs,1)

        for i in range(specs.nsspecs):
            S = self.blankfill(str(specs.sspecs[i]['value']),strsize1024)
            self.setNandWrite(S,1)

        thistimestep = self.storage.tlist[self.itimeStep]

        # Write out time, time step and cycle.
        L = thistimestep.time
        self.setNandWrite(L,1)
        L = thistimestep.dt
        self.setNandWrite(L,1)
        L = thistimestep.cycle
        self.setNandWrite(L,1)

        # Write out mesh size parameters

        L = self.mesh.cells['N']
        self.setNandWrite(L,1)
        L = self.mesh.vertices['N']
        self.setNandWrite(L,1)

        return

    def writebMesh(self):

        mesh = self.mesh
        try:
           import numpy.oldnumeric as Numeric
        except ImportError:
           import Numeric
        except:
           raise

        # Write out mesh connectivity, vertex info and mesh geometry info
        fp  = self.FP

        tmp = Numeric.transpose(mesh.cells['VERTICES'])
        for imv in range(len(tmp)):
            fp.write(data=tmp[imv],endRecord=1,debug=self.debug)

        ncblock = 0
        if mesh.cells.has_key('BLOCKID'):
            ncblock = 1
        self.setNandWrite(ncblock,1)
        if ncblock > 0:
            fp.write(data=mesh.cells['BLOCKID'],endRecord=1)

        tmp        = Numeric.transpose(mesh.vertices['COORDS'])
        for imv in range(len(tmp)):
            fp.write(data=tmp[imv],endRecord=1,debug=self.debug)

        return


    def writebVariables(self):

        try:
           import numpy.oldnumeric as Numeric
        except ImportError:
           import Numeric
        except:
           raise
        import string
        fp    = self.FP
        debug = self.debug

        #print out ZONE variables for standard truchas restart file

        # The Zone Writer
	if debug: print '  writing Zone variables '
        self.findAndWrite('Z_RHO',debug=debug)
        self.findAndWrite('Z_TEMP',debug=debug)
        self.findAndWrite('Z_ENTHALPY',debug=debug)
        self.findAndWrite('Z_P',debug=debug)
        self.findAndWrite('Z_VC',transpose=1,byComponent=1,debug=debug)
        self.findAndWrite('Face_Vel',transpose=1,byComponent=1,debug=debug)

        #special treatment for Z_VF

        #for x in self.vars:
        #    if 'Z_VF' in x.name:
        #        tdata = Numeric.transpose(x.data)

        #for j in range(len(tdata[0])):
        #    for i in range(len(tdata)):
        #        self.FP.write(tdata[i][j], 1)

        # the Matl Writer

	if debug: print >> self.fpwatch, '  Matl Writer'
        nmat     = 0
        novof    = 1
        for i in self.vars:
            if 'VOF' in i.name and 'M_VOF' not in i.name:
                nmat    = nmat + 1
                novof   = novof*0

        if novof: #if only 1 material set vof array = 1
            nmat = 1
            vf   = Numeric.array([1],'d')
            vf   = Numeric.resize(vf,[self.mesh.cells['N']])

        self.setNandWrite(nmat,1)

        if novof:
            self.FP.write(vf,1)
        else:
            for i in self.vars:
                if 'VOF' in i.name and 'M_VOF' not in i.name:
		    if debug: print >> self.fpwatch, '  writing ',i.name
                    self.FP.write(i.data, 1)

	if debug: print >> self.fpwatch, '  writing out features: ',self.storage.specs.feats
        # Sensitivity Writer..

        if ('sensitivity' in self.storage.specs.feats):
	    if debug: print >> self.fpwatch, '  Sensitivity Writer'
            for i in self.vars:
                if 'nsens_design_variables' in i.name:
                    ndesign    = i.data[0]
                    self.setNandWrite(ndesign,1)
                if 'nsens_functions' in i.name:
                    nfunctions = i.data[0]

            if ndesign > 0:
                temp_sens       = []
                enthalpy_sens   = []
                vol_frac_sens   = []
                for x in self.vars:
                    if 'temp_sens_' in x.name:
                        tdata         = Numeric.transpose(x.data)
                        temp_sens     = temp_sens + [tdata]
                    if 'enthalpy_sens_' in x.name:
                        tdata         = Numeric.transpose(x.data)
                        enthalpy_sens = enthalpy_sens + [tdata]
                    if 'volume_fraction_sens_' in x.name:
                        tdata         = Numeric.transpose(x.data)
                        vol_frac_sens = vol_frac_sens + [tdata]

                #assert len of temp_sens, enthalpy_sens = ndesign
                assert ndesign == len(temp_sens)
                assert ndesign == len(enthalpy_sens)

                #assert len of volume_frac_sens arrays = ndesign*nmat
                assert len(vol_frac_sens) == ndesign*nmat

                self.setNandWrite(int(len(vol_frac_sens)),1)

                for i in range(len(temp_sens)):
                    self.writeByComponent(temp_sens[i],1)
                for i in range(len(enthalpy_sens)):
                    self.writeByComponent(enthalpy_sens[i],1)
                for i in range(len(vol_frac_sens)):
                    self.writeByComponent(vol_frac_sens[i],1)

            self.setNandWrite(nfunctions,1)

            if nfunctions > 0:
                self.findAndWrite('sens_function_values',debug=debug)
                sens_grad = []
                for x in self.vars:
                    if 'sens_function_gradient_' in x.name:
                        sens_grad     = sens_grad + [x.data]

                #assert nfunctions = len of sens_grad
                assert nfunctions == len(sens_grad)
                for i in range(len(sens_grad)):
                    self.FP.write(sens_grad[i],endRecord=1)


        # Write out microstructure parameter derived types.
        # first check if they exist
        if ('alloy' in self.storage.specs.feats):

	    if debug: print >> self.fpwatch, '  alloy Writer'
            self.findAndWrite('LIQUIDUS_TEMP',debug=debug)
            self.findAndWrite('EUTECTIC_FRACTION',debug=debug)
            self.findAndWrite('SEC_ARM_SPACING',debug=debug)
            self.setNandWrite(self.storage.specs.nphases,1)
            self.setNandWrite(self.storage.specs.ncomps,1)

            phase_conc    = []

            for x in self.vars:
                if ('PC(' in x.name):
                    phase_conc = phase_conc + [x.data]

            for isph in range(len(phase_conc)):
                self.FP.write(phase_conc[isph], 1)

        # Write out solid mechanics derived types, check if they exist.
        if ('solid_mechanics' in self.storage.specs.feats):

	    if debug: print >> self.fpwatch, '  solid_mechanics Writer'
            # Note the calisthenics to get SMech_IP and SMech as it is in the code
            total_strain        = []
            elastic_stress      = []
            plastic_strain      = []
            plastic_strain_rate = []
            sigma               = []
            epsilon             = []
            eplastic            = []
            for x in self.vars:
                tdata = Numeric.transpose(x.data)
                if  ('TOTAL_STRAIN_' in x.name):
                    total_strain        = total_strain + [tdata]
                elif ('ELASTIC_STRESS_' in x.name):
                    elastic_stress      = elastic_stress + [ tdata]
                elif ('PLASTIC_STRAIN_' in x.name and not 'RATE' in x.name):
                    plastic_strain      = plastic_strain + [tdata]
                elif ('PLASTIC_STRAIN_RATE_' in x.name):
                    plastic_strain_rate = plastic_strain_rate + [tdata]
                elif (x.name == 'sigma'):
                    sigma               = tdata
                    if debug: print >> self.fpwatch, '    ',x.name,x.shape
                elif (x.name == 'epsilon'):
                    epsilon             = tdata
                    if debug: print >> self.fpwatch, '    ',x.name,x.shape
                elif (x.name == 'e_plastic'):
                    eplastic            = tdata
                    if debug: print >> self.fpwatch, '    ',x.name,x.shape
                elif (x.name == 'epsdot'):
                    epsdot              = tdata
                    if debug: print >> self.fpwatch, '    ',x.name,x.shape
                elif (x.name == 'Displacement'):
                    ndim                = x.shape[1]
                    if debug: print >> self.fpwatch, '    ',x.name,x.shape
                elif (x.name == 'RHS'):
                    rhs                 = tdata
                    if debug: print >> self.fpwatch, '    ',x.name,x.shape

            # Must do an error check to ensure that
            # all the above are same dimension

            assert len(total_strain)        == len(elastic_stress)
            if (len(total_strain)) > 0:
                assert len(total_strain[0]) == len(elastic_stress[0])
            assert len(total_strain)        == len(plastic_strain)
            if (len(total_strain)) > 0:
                assert len(total_strain[0]) == len(plastic_strain[0])
            assert len(plastic_strain_rate) == len(elastic_stress)

            for ismech in range(len(total_strain)):
                self.writeByComponent(total_strain[ismech], 2)
                self.writeByComponent(elastic_stress[ismech], 2)
                self.writeByComponent(plastic_strain[ismech], 2)
                self.writeByComponent(plastic_strain_rate[ismech], 1)
            self.writeByComponent(epsilon, 2)
            self.writeByComponent(sigma, 2)
            self.writeByComponent(eplastic, 2)
            self.writeByComponent(epsdot, 1)

            #special treatment RHS printed out for each dimension
            for i in range(ndim):
                tmp  = rhs[i::ndim]
                self.writeByComponent(tmp, 1)

            self.findAndWrite('epstherm',transpose=1,byComponent=1,debug=debug)
            self.findAndWrite('epspc',transpose=1,byComponent=1,debug=debug)
            self.findAndWrite('Displacement',transpose=1,byComponent=1,debug=debug)

	# Write out species data:
        #if ('species' in self.storage.specs.feats):
	if ('species' in self.storage.specs.feats):
	    if debug: print >> self.fpwatch, '  species Writer'
	    self.setNandWrite(self.storage.specs.nspecies,1)
	    for x in self.vars:
		if ('phi' in x.name):
                    self.findAndWrite(x.name,debug=debug)


        # Write out electromagnetics variables, check if they exist.
        if ('joule_heat' in self.storage.specs.feats):

	    if debug: print >> self.fpwatch, '  joule_heat Writer'
            self.findAndWrite('FREQ',debug=debug)
            self.findAndWrite('UHFS',debug=debug)

            #find ncoils,nmu,nsigma,njoule

            import string
            ncoils = 0
            nmu    = 0
            nsigma = 0
            njoule = 0
            for x in self.vars:
                if string.find(str(x.name),'COIL') >= 0:
                    if string.find(str(x.name),'CURRENT') >= 0:
                        ncoils = ncoils + 1
                if string.find(str(x.name),'MU') >= 0:
                    nmu    = x.shape[0]
                if string.find(str(x.name),'SIGMA') >= 0:
                    nsigma = x.shape[0]
                if string.find(str(x.name),'JOULE') >= 0:
                    njoule = x.shape[0]

            #write out ncoils
            if debug: print >> self.fpwatch,  '  ncoils'
            self.setNandWrite(ncoils,1)

            #write out characteristics of each coil

            for ncoil in range(ncoils):
                thiscoil = ncoil + 1
                coilchar = 'COIL ' + str(thiscoil) + ' CURRENT'
                self.findAndWrite(coilchar,debug=debug)
                coilchar = 'COIL ' + str(thiscoil) + ' CENTER'
                self.findAndWrite(coilchar,debug=debug)
                coilchar = 'COIL ' + str(thiscoil) + ' LENGTH'
                self.findAndWrite(coilchar,debug=debug)
                coilchar = 'COIL ' + str(thiscoil) + ' RADIUS'
                self.findAndWrite(coilchar,debug=debug)
                coilchar = 'COIL ' + str(thiscoil) + ' NTURNS'
                self.findAndWrite(coilchar,debug=debug)

            #write out nmu
            if debug: print >> self.fpwatch, '  writing nmu'
            self.setNandWrite(nmu,1)

            #write out mu
            self.findAndWrite('MU',debug=debug)

            #write out nsigma
            if debug: print >> self.fpwatch, '  writing nsigma'
            self.setNandWrite(nsigma,1)

            #write out sigma
            self.findAndWrite('SIGMA',debug=debug)

            #write out njoule
            if debug: print >> self.fpwatch, '  writing njoule'
            self.setNandWrite(njoule,1)

            #write out joule
            self.findAndWrite('JOULE',debug=debug)

    #now define ascii routines particular to RESTART

    def writeHeader(self):

        """Write restart header information to file."""

        sfmt  = '%%-%ss' % 256
        ifmt  = '%%%s' % '7i'
        dfmt  = self.dformat

        specs = self.storage.specs

        S     = sfmt % specs.fileid + '\n'
        S    += ifmt % specs.nfeats + '\n'

        for i in range(specs.nfeats):
            S += specs.feats[i] + '\n'

        S     += ifmt % specs.nsspecs + '\n'

        for i in range(specs.nsspecs):
            S += specs.sspecs[i]['value'] + '\n'

        "write out time step information "

        S     += dfmt % self.storage.tlist[self.itimeStep].time
        S     += dfmt % self.storage.tlist[self.itimeStep].dt
        S     += ifmt % self.storage.tlist[self.itimeStep].cycle + '\n'

        "write out mesh size parameters "

        S     += ifmt % self.mesh.cells['N']
        S     += ifmt % self.mesh.vertices['N'] + '\n'

        print >> self.FP, S

    def writeMesh(self):

        fp      = self.FP
        try:
           import numpy.oldnumeric as Numeric
        except ImportError:
           import Numeric
        except:
           raise

        dformat = self.dformat
        iformat = self.iformat
        debug   = self.debug
        mesh    = self.mesh
	
        # Write out the cells vertices.
        if debug: print >> self.fpwatch, '  Cell Vertices',
        sys.stdout.flush()
        c      = mesh.cells
        cVerts = Numeric.transpose(c['VERTICES'])
        for i in range(len(cVerts)):
            print >> fp, 'cell-vertices %i' %(i)
            self.writeArray(1,cVerts[i],'',iformat, fp, nperline=10)

        #write out block_ids for BRANCH_EM
        if debug: print >> self.fpwatch, '  Block IDs',
        sys.stdout.flush()
        ncblock     = 0
        if c.has_key('BLOCKID'):
            ncblock = 1
        print >> fp, 'ncblock'
        print >> fp, ncblock
        if ncblock > 0:
            cblocks = c['BLOCKID']
            print >> fp, 'block-data'
            self.writeArray(1,cblocks,'',iformat,fp,nperline=10)

        sys.stdout.flush()
        v           = mesh.vertices
        if debug: print >> self.fpwatch, '  Coords',
        sys.stdout.flush()
        vt          = Numeric.transpose(v['COORDS'])
        for i in range(len(vt)):
            print >> fp, 'Coords %i' %(i)
            self.writeArray(1,vt[i],'',dformat,fp)


    def writeVariables(self):

        s       = ''
        fp      = self.FP
        dformat = self.dformat
        iformat = self.iformat
        debug   = self.debug

        try:
           import numpy.oldnumeric as Numeric
        except ImportError:
           import Numeric
        except:
           raise
        import string

	if debug: print '  writing Zone variables '
        for x in self.vars:
            if (x.name.startswith('Z_') and not('OLD' in x.name)):
                if debug: print >> self.fpwatch, '    ',x.name
                print >> fp, x.name
                if (x.rank > 1):
                    x.data = Numeric.transpose(x.data)
                if x.type == 'i': format = iformat
                if x.type == 'd': format = dformat
                #special treatment for Z_VF
                if 'Z_VF' in x.name:
                    for j in range(len(x.data[0])):
                        for i in range(len(x.data)):
                            s    = str(j) + str(i)
                            if debug: print >> self.fpwatch, '    ',s
                            print >> fp, s
                            self.writeArray(1,x.data[i][j],'',format,fp,10)
                else:
                    self.writeArray(x.rank,x.data,'',format,fp,10)

        #face velocity now...

        for x in self.vars:
            if (x.name == 'Face_Vel'):
                if debug: print >> self.fpwatch, '    ',x.name
                print >> fp, x.name
                x.data = Numeric.transpose(x.data)
                self.writeArray(x.rank,x.data,'',dformat,fp,10)

        #matl quantities now

	if debug: print >> self.fpwatch, '  Matl Writer'
        nmat     = 0
        novof    = 1
        for i in self.vars:
            if 'VOF' in i.name and 'M_VOF' not in i.name:
                nmat   += 1
                novof   = novof*0

        if novof: #if only 1 material set vof array = 1
            nmat = 1
            vf   = Numeric.array([1],'d')
            vf   = Numeric.resize(vf,[self.mesh.cells['N']])

        print >> fp, nmat

        if novof and len(self.vars):
	    if debug: print >> self.fpwatch, '    E>>a ',i.name
            print >> fp, i.name	
            self.writeArray(1,vf,'',dformat,fp,10)
        else:
            for i in self.vars:
                if 'VOF' in i.name and 'M_VOF' not in i.name:
		    if debug: print >> self.fpwatch, '    E>>b ',i.name
		    print >> fp, i.name		
                    self.writeArray(i.rank,i.data,'',dformat,fp,10)

	if debug: print >> self.fpwatch, '  writing out features: ',self.storage.specs.feats
        # Sensitivity Writer..

        if ('sensitivity' in self.storage.specs.feats):
	    if debug: print >> self.fpwatch, '  Sensitivity Writer'
            for i in self.vars:
                if 'nsens_design_variables' in i.name:
                    if debug: print >> self.fpwatch, '    n_sens_design_variables %i' %(i.data[0])
                    print >> fp, 'n_sens_design_variables %i' %(i.data[0])
                    ndesign    = i.data[0]
                if 'nsens_functions' in i.name:
                    nfunctions = i.data[0]

            if ndesign > 0:
                temp_sens       = []
                enthalpy_sens   = []
                vol_frac_sens   = []
                for x in self.vars:
                    if 'temp_sens_' in x.name:
                        tdata         = Numeric.transpose(x.data)
                        temp_sens     = temp_sens + [tdata]
                    if 'enthalpy_sens_' in x.name:
                        tdata         = Numeric.transpose(x.data)
                        enthalpy_sens = enthalpy_sens + [tdata]
                    if 'volume_fraction_sens_' in x.name:
                        tdata         = Numeric.transpose(x.data)
                        vol_frac_sens = vol_frac_sens + [tdata]

                #assert len of temp_sens, enthalpy_sens = ndesign
                assert ndesign == len(temp_sens)
                assert ndesign == len(enthalpy_sens)

                #assert len of volume_frac_sens arrays = ndesign*nmat
                assert len(vol_frac_sens) == ndesign*nmat

                for i in range(len(temp_sens)):
                    if debug: print >> self.fpwatch, '    temp_sens %i' %(i+1)
                    print >> fp, 'temp_sens %i' %(i+1)
                    self.writeArray(1,temp_sens[i],'',format,fp,10)
                for i in range(len(enthalpy_sens)):
                    if debug: print >> self.fpwatch, '    enthalpy_sens %i' %(i+1)
                    print >> fp, 'enthalpy_sens %i' %(i+1)
                    self.writeArray(1,enthalpy_sens[i],'',format,fp,10)
                for i in range(len(vol_frac_sens)):
                    if debug: print >> self.fpwatch, '    vol_frac_sens %i' %(i+1)
                    print >> fp, 'vol_frac_sens %i' %(i+1)
                    self.writeArray(1,vol_frac_sens[i],'',format,fp,10)

            if debug: print >> self.fpwatch, '    nsens_functions %i' %(nfunctions)
            print >> fp, 'nsens_functions %i' %(nfunctions)

            if nfunctions > 0:
                sens_grad = []
                for x in self.vars:
                    if 'sens_function_values' in x.name:
                        if debug: print >> self.fpwatch, '   ',x.name
                        print >> fp, x.name, x.rank
                        self.writeArray(x.rank,x.data,'',format,fp,10)

                    if 'sens_function_gradient_' in x.name:
                        sens_grad     = sens_grad + [x.data]

                #assert nfunctions = len of sens_grad
                assert nfunctions == len(sens_grad)

                for i in range(len(sens_grad)):
                    s        = 'sens_grad %i' %(i+1)
                    if debug: print >> self.fpwatch, '   ',s
                    print >> fp, s
                    self.writeArray(1,sens_grad[i],'',format,fp,10)

        #for the remaining quantities assume float format

        format = dformat

        # Write out microstructure parameter derived types.
        # first check if they exist
        if ('alloy' in self.storage.specs.feats):

	    if debug: print >> self.fpwatch, '  alloy Writer'
            for x in self.vars:
                if (x.name == 'LIQUIDUS_TEMP'):
                    if debug: print >> self.fpwatch, '    ',x.name
                    print >> fp, x.name
                    self.writeArray(x.rank,x.data,'',format,fp,10)
            for x in self.vars:
                if (x.name == 'EUTECTIC_FRACTION'):
                    if debug: print >> self.fpwatch, '    ',x.name
                    print >> fp, x.name
                    self.writeArray(x.rank,x.data,'',format,fp,10)
            for x in self.vars:
                if (x.name == 'SEC_ARM_SPACING'):
                    if debug: print >> self.fpwatch, '    ',x.name
                    print >> fp, x.name
                    self.writeArray(x.rank,x.data,'',format,fp,10)

            print >> fp, self.storage.specs.nphases
            print >> fp, self.storage.specs.ncomps

            for x in self.vars:
                if ('PC(' in x.name):
                    if debug: print >> self.fpwatch, '    ',x.name
                    print >> fp, x.name, x.rank
                    self.writeArray(x.rank,x.data,'',format,fp,10)

        # Write out solid mechanics data :
        if ('solid_mechanics' in self.storage.specs.feats):

	    if debug: print >> self.fpwatch, '  solid_mechanics Writer'
            total_strain        = []
            elastic_stress      = []
            plastic_strain      = []
            plastic_strain_rate = []
            sigma               = []
            epsilon             = []
            eplastic            = []

            for x in self.vars:
                tdata = Numeric.transpose(x.data)
                if  ('TOTAL_STRAIN_' in x.name):
                    total_strain        = total_strain + [tdata]
                elif ('ELASTIC_STRESS_' in x.name):
                    elastic_stress      = elastic_stress + [ tdata]
                elif ('PLASTIC_STRAIN_' in x.name and not 'RATE' in x.name):
                    plastic_strain      = plastic_strain + [tdata]
                elif ('PLASTIC_STRAIN_RATE_' in x.name):
                    plastic_strain_rate = plastic_strain_rate + [tdata]
                elif (x.name == 'sigma'):
                    sigma               = tdata
                    if debug: print >> self.fpwatch, '    ',x.name,x.shape
                elif (x.name == 'epsilon'):
                    epsilon             = tdata
                    if debug: print >> self.fpwatch, '    ',x.name,x.shape
                elif (x.name == 'e_plastic'):
                    eplastic            = tdata
                    if debug: print >> self.fpwatch, '    ',x.name,x.shape
                elif (x.name == 'epsdot'):
                    epsdot              = tdata
                    if debug: print >> self.fpwatch, '    ',x.name,x.shape
                elif (x.name == 'epstherm'):
                    epstherm            = tdata
                    if debug: print >> self.fpwatch, '    ',x.name,x.shape
                elif (x.name == 'epspc'):
                    epspc               = tdata
                    if debug: print >> self.fpwatch, '    ',x.name,x.shape
                elif (x.name == 'Displacement'):
                    ndim                = x.shape[1]
                    disp                = tdata
                elif (x.name == 'RHS'):
                    rhs                 = tdata

            for ismech in range(len(total_strain)):
                s        = 'total_strain %i' %(ismech)
                if debug: print >> self.fpwatch, '    ',s
                print >> fp, s
                self.writeArray(2,total_strain[ismech],'',format,fp,10)
                s        = 'elastic_stress %i' %(ismech)
                if debug: print >> self.fpwatch, '    ',s
                print >> fp, s
                self.writeArray(2,elastic_stress[ismech],'',format,fp,10)
                s        = 'plastic_strain %i' %(ismech)
                if debug: print >> self.fpwatch, '    ',s
                print >> fp, s
                self.writeArray(2,plastic_strain[ismech],'',format,fp,10)
                s        = 'plastic_strain_rate %i' %(ismech)
                if debug: print >> self.fpwatch, '    ',s
                print >> fp, s
                self.writeArray(1,plastic_strain_rate[ismech],'',format,fp,10)

            s        = 'epsilon '
            if debug: print >> self.fpwatch, '    ',s
            print >> fp, s
            self.writeArray(2, epsilon, '',format,fp,10)
            s        = 'sigma '
            if debug: print >> self.fpwatch, '    ',s
            print >> fp, s
            self.writeArray(2, sigma, '',format,fp,10)
            s        = 'eplastic '
            if debug: print >> self.fpwatch, '    ',s
            print >> fp, s
            self.writeArray(2, eplastic, '',format,fp,10)
            s        = 'epsdot '
            if debug: print >> self.fpwatch, '    ',s
            print >> fp, s
            self.writeArray(1, epsdot, '',format,fp,10)

            #special treatment RHS printed out for each dimension

            for i in range(ndim):
                s        = 'rhs %i ' %(i)
                if debug: print >> self.fpwatch, '    ',s
                print >> fp, s
                tmp  = rhs[i::ndim]
                self.writeArray(1, tmp, '',format,fp,10)

            s        = 'epstherm '
            if debug: print >> self.fpwatch, '    ',s
            print >> fp, s
            self.writeArray(2, epstherm, '',format,fp,10)

            s        = 'epspc '
            if debug: print >> self.fpwatch, '    ',s
            print >> fp, s
            self.writeArray(2, epspc, '',format,fp,10)

            s        = 'Displacement '
            if debug: print >> self.fpwatch, '    ',s
            print >> fp, s
            self.writeArray(2, disp, '',format,fp,10)

	# Write out species data:
        #if ('species' in self.storage.specs.feats):
	if ('species' in self.storage.specs.feats):
	    if debug: print >> self.fpwatch, '  species Writer'
            print >> fp, 'nspecies %i' %(self.storage.specs.nspecies)
 	    for x in self.vars:
	        if ('phi' in x.name):
		    if debug: print >> self.fpwatch, '    ',x.name
		    print >> fp, x.name
                    if (x.rank > 1):
                        x.data = Numeric.transpose(x.data)
                    if x.type == 'i': format = iformat
                    if x.type == 'd': format = dformat
                    self.writeArray(x.rank,x.data,'',format,fp,10)

        # Write out electromagnetics variables, check if they exist.

        if ('joule_heat' in self.storage.specs.feats):

	    if debug: print >> self.fpwatch, '  joule_heat Writer'
            for x in self.vars:
                if x.name == 'FREQ':
                    if debug:
                        print >> self.fpwatch, '    ',x.name
                        #print >> self.fpwatch, x.data
                    print >> fp, x.name
                    print >> fp, x.data[0]
                if x.name == 'UHFS':
                    if debug:
                        print >> self.fpwatch, '    ',x.name
                        #print >> self.fpwatch, x.data
                    print >> fp, x.name
                    print >> fp, x.data[0]

            #find ncoils,nmu,nsigma,njoule

            import string
            ncoils = 0
            for x in self.vars:
                if string.find(str(x.name),'COIL') >= 0:
                    if string.find(str(x.name),'CURRENT') >= 0:
                        ncoils = ncoils + 1
                if string.find(str(x.name),'MU') >= 0:
                    nmu = x.shape[0]
                if string.find(str(x.name),'SIGMA') >= 0:
                    nsigma = x.shape[0]
                if string.find(str(x.name),'JOULE') >= 0:
                    njoule = x.shape[0]

            #write out ncoils

            if debug:
                print >> self.fpwatch, '  ncoils'
                #print >> self.fpwatch, ncoils
            print >> fp, 'ncoils'
            print >> fp, ncoils

            #write out characteristics of each coil

            for ncoil in range(ncoils):
                coilchar = {}
                thiscoil = ncoil + 1
                coilchar['CURRENT'] = 'COIL ' + str(thiscoil) + ' CURRENT'
                coilchar['CENTER']  = 'COIL ' + str(thiscoil) + ' CENTER'
                coilchar['LENGTH']  = 'COIL ' + str(thiscoil) + ' LENGTH'
                coilchar['RADIUS']  = 'COIL ' + str(thiscoil) + ' RADIUS'
                coilchar['NTURNS']  = 'COIL ' + str(thiscoil) + ' NTURNS'

            for x in self.vars:
                if x.name == coilchar['CURRENT']:
                    if debug:
                        print >> self.fpwatch, '    ',x.name
                        #print >> self.fpwatch, x.data[0]
                    print >> fp, x.name
                    print >> fp, x.data[0]
                if x.name == coilchar['CENTER']:
                    if debug:
                        print >> self.fpwatch, '    ',x.name
                        #print >> self.fpwatch, x.data[0]
                    print >> fp, x.name
                    self.writeArray(x.rank,x.data,'',format,fp)
                if x.name == coilchar['LENGTH']:
                    if debug:
                        print >> self.fpwatch, '    ',x.name
                        #print >> self.fpwatch, x.data[0]
                    print >> fp, x.name
                    print >> fp, x.data[0]
                if x.name == coilchar['RADIUS']:
                    if debug:
                        print >> self.fpwatch, '    ',x.name
                        #print >> self.fpwatch, x.data[0]
                    print >> fp, x.name
                    print >> fp, x.data[0]
                if x.name == coilchar['NTURNS']:
                    if debug:
                        print >> self.fpwatch, '    ',x.name
                        #print >> self.fpwatch, x.data[0]
                    print >> fp, x.name
                    print >> fp, x.data[0]

            #write out nmu
            if debug: print >> self.fpwatch, '  nmu'
            print >> fp, 'nmu'
            print >> fp, nmu

            #write out mu

            for x in self.vars:
                if x.name == 'MU':
                    if debug: print >> self.fpwatch, x.name
                    print >> fp, x.name
                self.writeArray(x.rank,x.data,'',format,fp)

            #write out nsigma

            if debug: print >> self.fpwatch, '  nsigma'
            print >> fp, 'nsigma'
            print >> fp, nsigma

            #write out sigma

            for x in self.vars:
                if x.name == 'SIGMA':
                    if debug: print >> self.fpwatch, x.name
                    print >> fp, x.name
                self.writeArray(x.rank,x.data,'',format,fp)

            #write out njoule

            if debug: print >> self.fpwatch, '  njoule'
            print >> fp, 'njoule'
            print >> fp, njoule

            #write out joule

            for x in self.vars:
                if x.name == 'JOULE':
                    if debug: print >> self.fpwatch, '    ',x.name
                    print >> fp, x.name
                self.writeArray(x.rank,x.data,'',format,fp)


if __name__ =='__main__':

    'for component testing RESTART write component'

    dfltdir = '../../scripts/test_TBrookParse/samples/'
    
    from PYTHONutils import uTestOpts

    # file, debug, output prefix, binary, ascii, clean
    opts = uTestOpts('fdoc', 
                     defaults = {'c' : False,
                                 'o' : 'restartUnitTest'},
                     actions  = {'c' : 'store_true'},
                     dir      = dfltdir)
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    fpwatch    = sys.stdout
    
    try:

        if opt.c:
            print >> fpwatch, "Cleaning up from a past test run...",
            cmd = 'rm -f %s' %(dfltdir+opt.o+'.binary')
            os.system(cmd)
            cmd = 'rm -f %s' %(dfltdir+opt.o+'.ascii')
            os.system(cmd)
            print >> fpwatch, "Cleanup successful!!"

        else:
            from XMLgetStorageObject import getStorageObject
            storage   = getStorageObject(opt.f)	
            mesh      = storage.mlist[0] #assume default mesh
            
            #now get the data values for the mesh fields
            for meshspace in mesh.mslist:
                meshspace.vlist = storage.getValues(meshspace.vlist)

            # fill in actual data values
            mesh.fillMesh()
            
            itimeout  = 1
            t         = storage.tlist[itimeout].time
            vars      = storage.tlist[itimeout].vlist
            vars      = storage.getValues(vars)

	    # ASCII test
            print >> fpwatch, '*'*60
            print >> fpwatch, 'Testing ascii write...'
            fileout = dfltdir + opt.o + '.ascii'
            writer  = RESTARTwriteStorageObject(fileout, 'ascii',
                                                storage, mesh, itimeout, t,
                                                vars, debug=opt.d)
            print >> fpwatch, 'Restart file %s\n  written to %s' % \
			      (opt.o+'.ascii',dfltdir)

	    # Binary test
            print >> fpwatch, '*'*60
            print >> fpwatch, 'Testing binary write...'
            fileout = dfltdir + opt.o + '.binary'
            writer  = RESTARTwriteStorageObject(fileout, 'binary',
                                                storage, mesh, itimeout, t,
                                                vars, debug=opt.d)
            print >> fpwatch, 'Restart file %s\n  written to %s' % \
			      (opt.o+'.binary',dfltdir)
                
            print >> fpwatch, '*'*60
            print >> fpwatch, \
		'--> Note: to clean up, re-execute same command line plus "-c"'

    except:
        print >> fpwatch, "---> Test failed in some aspect <---"
        print >> fpwatch, "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise



