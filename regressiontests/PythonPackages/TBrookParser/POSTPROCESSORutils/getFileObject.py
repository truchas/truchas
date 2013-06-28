

"""

 getFileObject

 -----------------------------------------------------------------------------
  Purpose:
  
    Instantiate a FileObject.
     
    Provides methods for choosing an available mesh and
    for printing FileObject information
  
  Public Interface(s):
  
    fo = getFileObject(userinput,default_file)
    m  = fo.chooseMesh(storage)
    fo.str()

  Contains:
    class getFileObject
        __init__(self,userinput,default_file,mapfromfiles=[],
                 files_storages={},stepids=[0],mapcount=0,debug=0)
        chooseMesh(self,storage)
        str(self)

    Unit Test Block
  
  Author(s): Sharen Cummins (scummins@lanl.gov)
 -----------------------------------------------------------------------------
"""
import os, sys, string
if __name__=='__main__':
    " Set sys.path for component testing mode "
    thisdir   = os.path.dirname(__file__)
    # modify dot-dots appropriately
    parsedir  = os.path.abspath(thisdir+'../') 

    sys.path.append(parsedir)

class getFileObject:

    def __init__(self,userinput,default_file,mapfromfiles=[],
                 files_storages={},stepids=[0],mapcount=0,
		 fpwatch=sys.stdout,debug=0,options=None):

        "initialise the file object"
	if debug: print 'Debug - welcome to getFileObject'

        self.file     = None # name of the file used to create this file object
        self.storage  = None # storage object created from this file
        self.regions  = []   # list of regions defined on the meshes in this file
        self.input    = userinput # user input mechanism
        self.debug    = debug
        self.fpwatch  = fpwatch
        
        "now start defining the attributes in the file object"

        self.file = default_file
        f         = self.file

        try:
            if string.find(f,"xml") >= 0:
                # its an xml file
		if self.debug: print '  Debug - in getFileObject: found xml file ',f
                from XMLgetStorageObject import getStorageObject
                self.storage   = getStorageObject(f,fp=self.fpwatch,
				                  debug=self.debug)
		if self.debug: print '    Debug - in getFileObject: got SO'
                themesh        = self.chooseMesh(self.storage)
		if self.debug: print '    Debug - in getFileObject: chose Mesh'
            elif string.find(f,"exo") >= 0:
                # its an exodus file
		if self.debug: print '  Debug - in getFileObject: found exo file ',f
                from MODStorageObject import modStorageObject
		if options==None:
                    f     = self.input.usetc(
			     'Name of file to be read:', default_file)
                L         = string.split(f,'/')
                exofile   = L[len(L)-1:]
                M2        = L[0:len(L)-1]
                M2        = string.join(M2,'/')
                self.file = 'mapping' + str(mapcount)
                mapcount += 1

                for f2 in mapfromfiles:
                    if string.find(f2,"xml") >= 0:
                        #create new filename 
                        L             = string.split(f2,'/')
                        xmlfile       = L[len(L)-1:]
                        M             = L[0:len(L)-1]
                        M             = string.join(M,'/')
                if len(mapfromfiles) > 1:
                    #for multiple mesh mappings the newly created mapped object
                    #is placed in Exodus file directory
                    if len(M2):
                        self.file  = M2 + '/' + self.file  
                else:
                    #for single mesh mappings the newly created mapped object
                    #is placed in Truchas output directory
                    if len(M):
                        self.file  = M + '/' + self.file  
                    
                #check if a storage object has already 
		#been created from this file..
                #if so then append the storage object with 
		#data from the newly chosen stepid

                if files_storages.has_key(self.file):
                    self.storage = files_storages[self.file][0]

                tstorages = []
                tmeshes   = []
                for f2 in mapfromfiles:
                    from XMLgetStorageObject import getStorageObject
                    tstorage  = getStorageObject(f2,fp=self.fpwatch,
				                 debug=self.debug)
                    tmesh     = self.chooseMesh(tstorage)
                    tstorages.append(tstorage)
                    tmeshes.append(tmesh)
                f3            = os.path.abspath(f)
                print >> self.fpwatch

		if options == None:
                    prompt = 'Specify a coordinate scale factor to scale ' + \
			     'the new mesh geometry:'
                    csf    = self.input.usetf(prompt,tstorages[0].specs.csf)
	        else:
		    if options['scale'] == None:
			options['scale'] = tstorages[0].specs.csf
		    csf    = options['scale']

                self.storage  = modStorageObject(tstorages,tmeshes,self.storage,
				                 f3,csf,stepids,self.fpwatch,
						 self.debug)
                themesh       = self.storage.mlist[0]
                
                #now obtain the following mesh storage structures for 'themesh'
                #themesh.cells
                #themesh.vertices
                #themesh.faces
                #themesh.edges

                for meshspace in themesh.mslist:
                    self.storage.getValues(meshspace.vlist)
                        
                themesh.fillMesh()

            else:
                print >> self.fpwatch
                print >> self.fpwatch,'Currently we can only map from XML files'
                print >> self.fpwatch,'Restart file will not be created'
                print >> self.fpwatch
                return

            #now obtain default regions having loaded in this file
            from getRegionObject import getRegionObject
            region            = getRegionObject(self.input,
			                        file=self.file,
			                        storage=self.storage,
                                                mesh=themesh,
						selector=['Default'],
						debug=self.debug)
            self.regions.append(region)
        except:
            print >> self.fpwatch
            print >> self.fpwatch, 'Unable to open file '
            print >> self.fpwatch, 'Filename:', self.file
            print >> self.fpwatch, 'Either the file does not exist, ',
	    print >> self.fpwatch, 'or it is not an XML or Exodus II file.'
            print >> self.fpwatch
            sys.exit(1)

    def chooseMesh(self,storage):

        import sys
        #for XML files with multiple meshes the user 
	#must choose which mesh to load

	if self.debug: print '    Debug - getFileObject.chooseMesh'
        count = 0
        if len(storage.mlist) > 1:
            meshes     = ''
            for m in storage.mlist:
                meshes = meshes + str(m.name)
            prompt     = 'Available meshes: '+ meshes + \
			 ' Which mesh do you want to load? '
            meshchoice = self.input.usetc(prompt,str(storage.mlist[0].name))
            try:
                cnt    = 0
                for m in storage.mlist:
                    if meshchoice == m.name:
                        count     = cnt
                    cnt    = cnt + 1
            except:
                print >> self.fpwatch
                print >> self.fpwatch, 'Invalid entry.  Must pick one of: ',
		print >> self.fpwatch,  meshes
                print >> self.fpwatch, 'Aborting this region definer.'
                print >> self.fpwatch
                sys.exit(1)

        mesh = storage.mlist[count]

	if self.debug: print '    Debug - getFileObject.chooseMesh: storage.getValues'
        for meshspace in mesh.mslist:
            storage.getValues(meshspace.vlist)
                        
	if self.debug: print '    Debug - getFileObject.chooseMesh: mesh.fillMesh'
        mesh.fillMesh()

        storage.mlist[count] = mesh
        
	if self.debug: print '    Debug - getFileObject.chooseMesh: done'
        return storage.mlist[count]
        
    def str(self):

        print >> self.fpwatch, 'loaded file'
        print >> self.fpwatch, self.file
        print >> self.fpwatch, 'storage'
        print >> self.fpwatch, self.storage.str()
        print >> self.fpwatch, 'regions'
        for reg in self.regions:
            reg.str()
        print >> self.fpwatch, 'finished'

if __name__=='__main__':
    
    'for testing this component'

    from PYTHONutils import uTestOpts
    dir  = '../../../../tools/scripts/test_TBrookParse/samples/'
    file = 'map_output/map.TBrook.xml'
    opts = uTestOpts('fd', defaults={'f': dir+file})
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    fpwatch = sys.stdout

    try:
        import usubs
        inBuffer   = ''
        input      = usubs.input(inBuffer)
        
        fileobject = getFileObject(input, opt.f, fpwatch=fpwatch, debug=opt.d)
        fileobject.str()

    except:
        print >> fpwatch, "---> Test failed in some aspect <---"
        print >> fpwatch, "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise
    





