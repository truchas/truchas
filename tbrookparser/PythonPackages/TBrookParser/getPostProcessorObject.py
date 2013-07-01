"""

 getPostProcessorObject

 -----------------------------------------------------------------------------
  Purpose:

    Instantiates a postprocessor object.

  Public Interface(s):

    Note: some of these may be intended for internal use only.

    PPO = getPostProcessorObject(def_file,def_exofile,userinput)
    PPO.defaults()
    PPO.load(file)        
    PPO.listfiles()        
    PPO.listregions()      
    PPO.listvariables()    
    PPO.probe(X)          
    PPO.query(X)          
    PPO.stat(X)           
    PPO.timeseries(X)     
    PPO.write(X,choice)   
    PPO.restart(Xs)       
    PPO.createregion()     
    PPO.deleteregions()    
    PPO.help(commands)    
    PPO.str()

  Contains:
    class getPostProcessorObject
        __init__(self,def_file,def_exofile,userinput,debug=0)
        defaults(self)
        load(self,file)
        listfiles(self)
        listregions(self)
        listvariables(self)
        probe(self,X)
        query(self,X)
        stat(self,X)
        timeseries(self,X)
        write(self,X,choice)
        restart(self,Xs)
        createregion(self)
        deleteregions(self)
        help(self,commands)
        str(self)
        __chooseSteps(X)
        __chooseRestartSteps(X,type)

    Unit Test Block

  Version:
    $ID$

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

from POSTPROCESSORutils        import getFileObject,getRegionObject,usubs
from writeStorageObject        import writeStorageObject
from ENSIGHTwriteStorageObject import ENSIGHTwriteStorageObject 
from GMVwriteStorageObject     import GMVwriteStorageObject
from RESTARTwriteStorageObject import RESTARTwriteStorageObject 
from TECPLOTwriteStorateObject import TECPLOTwriteStorageObject 
from VTKwriteStorageObject     import VTKwriteStorageObject 
from PYTHONutils               import getdir

class getPostProcessorObject:

    def __init__(self,def_file,def_exofile,userinput,fp=sys.stdout,debug=0):

        #initialise the postprocessor object

        self.files          = []     #list of files loaded 
        self.files_storages = {}     #a list of storage objects for each file
        self.files_regions  = {}     #a list of region objects for each file
        self.input          = userinput   #user input mechanism
        self.def_file       = def_file    #default file if not specified
        self.def_exofile    = def_exofile #default Exodus II file 
	                                  #    if not specified
        self.mapcount       = 0           #number of times mapping process is 
	                                  #    performed
        self.fp             = fp
        self.debug          = debug

        #now start defining the attributes and functions in the 
	#postprocessor object

        self.varformat = 3
        self.precision = '5E'
        self.coordsfmt = '%12.'+self.precision
        self.fieldfmt  = '%20.'+self.precision

	#flag added to skip IO opertations
	self.noIO = 0
	self.options = {}

    def defaults(self):
        
        """
        developers can define defaults...
        """

        oupt = self.input.usetc('Choose what default you would like to set ','precision')

        if oupt == 'precision':
            self.precision = self.input.usetc('Precision format for all diagnostic output ','5e')
            self.coordsfmt = '%12.' + self.precision
            self.fieldfmt  = '%20.' + self.precision
        else:
            print >> self.fp
            print >> self.fp, 'Currently only precision can be set as a default' 
            print >> self.fp
            return

    def load(self,file):

	if self.debug: 
	    print 'Debug > in getPPObject.load'
	    print '    Debug > in getPPObject.load: calling getFileObject'
        X = getFileObject(self.input,file,fpwatch=self.fp,debug=self.debug,options=self.options)
        if not (X.file in self.files):
            self.files.append(X.file)
        
        if self.files_storages.has_key(X.file):
            self.files_storages[X.file].append(X.storage)
        else:
            self.files_storages[X.file] = [X.storage]

        self.files_regions[X.file]      = X.regions
	if self.debug: print '    Debug > in getPPObject.load: succesfully loaded fileObject'

        return X

    def listfiles(self):

        print >> self.fp
        print >> self.fp, 'Loaded files are:'
        print >> self.fp
        for file in self.files:
            print >> self.fp, file
        print >> self.fp 

        return


    def listregions(self):

        print >> self.fp
        print >> self.fp, 'Regions are:'
        print >> self.fp
        for thisfile in self.files:
            for region in self.files_regions[thisfile]:
                print >> self.fp, region.describe
        print >> self.fp 

        return


    def listvariables(self):


        for thisfile in self.files:

            print >> self.fp, '*'*60,'\nFor file:', thisfile, '*'*60
            print >> self.fp 

            names_for_cycle = {}
            times_for_cycle = {}
            cycles_for_name = {}
            cyclenumbers    = []
            
            for storage in self.files_storages[thisfile]:
                for step in storage.tlist:
                    cyc                              = step.cycle
                    t                                = step.time
                    cyclenumbers                     = cyclenumbers + [cyc]
                    names_for_cycle[cyc]             = []
                    times_for_cycle[cyc]             = step.time
                    
                    print >> self.fp, '---------------------------------------------'
                    print >> self.fp, 'Cycle:',cyc,',   Time',t
                    print >> self.fp, '---------------------------------------------'
                    
                    for x in step.vlist:
                        names_for_cycle[cyc]        = names_for_cycle[cyc] + \
                                                      [x.nickname]
                        if (not cycles_for_name.has_key(x.name)):
                            cycles_for_name[x.name] = []
                        cycles_for_name[x.name]     = cycles_for_name[x.name] +\
					              [cyc]
                    names_for_cycle[cyc].sort()

                    l     = 0
                    for x in names_for_cycle[cyc]:
                        l = max(l,len(x))
                    s = '%%-%ds  '%l
    
                    t = ''
                    for x in names_for_cycle[cyc]:
                        xn = s%x
                        if(len(t) + len(xn) > 80):
                            print >> self.fp, t
                            t = ''
                        t = t + xn
                    print >> self.fp, t + '\n' # takes care of last t

                print >> self.fp
                print >> self.fp, 'Variable information'
                print >> self.fp
                    
                keys   = cycles_for_name.keys()
                keys.sort()
                l      = 0
                for x in keys:
                    l  = max(l,len(x))
                    xt = [cycles_for_name[x][0]]
                    for i in range(len(cycles_for_name[x])-1):
                        if(cycles_for_name[x][i] != cycles_for_name[x][i+1]):
                            xt                   += [cycles_for_name[x][i+1]]
                    cycles_for_name[x] = xt
                s      = '%%-%ds exists in cycles: '%l
                for x in keys:
                    t = s%x
                    print >> self.fp, t, cycles_for_name[x]
            
        return

    def probe(self,X):

        import sys
        try:
           import numpy.oldnumeric as Numeric
        except ImportError:
           import Numeric
        except:
           raise

        iformat         = '%15i'

        this_storage    = X.storage
        count           = 0
        prbelist        = ''
        for prbe in this_storage.plist:
            s           = '%+25s' %(str(prbe.name))
            s           = s + ','
            prbelist    = prbelist + s
            count       = count + 1
            if (count%self.varformat == 0):
                prbelist = prbelist + '\n'
        prbelist = prbelist[0:len(prbelist)-1]
        prompt   = 'Choose which probe? Please enclose the choice' + ' in double quotes. \n [%s]' %(prbelist)
        prbename = self.input.usetc(prompt,str(this_storage.plist[0].name))

        noprbe   = 1
        for k in this_storage.plist:
            if (k.name == prbename):
                theprobe = k
                noprbe   = noprbe*0

        if noprbe:
            print >> self.fp
            print >> self.fp, 'Your chosen probe is not in this file'
            print >> self.fp, 'Aborting this probe.'
            print >> self.fp
            return

        prompt          = '\n Output to "screen" or "file"? '
        oupt            = self.input.usetc(prompt, 'screen')
        if oupt == 'file':
            prompt      = '\n Filename to output results to?'
            fname       = self.input.usetc(prompt, 'probe.txt')
            fp          = open(fname,'w')
        else:
            fp          = self.fp
            
        count           = 0
        varlist         = ''
        for var in theprobe.vlist:
            s           = '%+25s' %(str(var.nickname))
            s           = s + ','
            varlist     = varlist + s
            count       = count + 1
            if (count%self.varformat == 0):
                varlist = varlist + '\n'
        varlist     = varlist[0:len(varlist)-1]
        prompt      = '\n Probe which variable(s)? \n [%s]' %(varlist)
        thevars     = '%s,%s' % \
		      (theprobe.vlist[0].nickname,theprobe.vlist[1].nickname)
        varnames    = self.input.usetc(prompt,thevars)
        varnames    = varnames.split(',')

        if (len(varnames)<1):
            print >> self.fp
            print >> self.fp, 'No variable specified in this probe diagnostic'
            print >> self.fp
            return

        thevars = []
        
        for k in theprobe.vlist:

            for j in varnames:

                if (k.nickname == j):

                    thevars.append(k)

        thevars   = this_storage.getValues(thevars)

        tolerance = 1.0e-14
        mintime   = 1.0e+14
        maxtime   = -1.0e+14
        times     = []
        
        for var in thevars:
            #add a small tolerance to the time column of each variable's 
	    #data array
            tmp     = var.data[:,1]
            tmp2    = Numeric.array([tolerance],var.type)
            tmp2    = Numeric.resize(tmp2,Numeric.shape(tmp))
            tmp3    = tmp + tmp2
            times.append(tmp3)
        
            mintime = min(tmp[Numeric.argmin(tmp)],mintime)
            maxtime = max(tmp[Numeric.argmax(tmp)],maxtime)
 
        strmintime = self.coordsfmt%mintime
        strmintime = strmintime[0:len(strmintime)]
        strmaxtime = self.coordsfmt%maxtime
        strmaxtime = strmaxtime[0:len(strmaxtime)]
        
        prompt  = '\n Time range to analyze the probe? [%s->%s]' %(strmintime,strmaxtime)
        trange  = '%s->%s' %(strmintime,strmaxtime)
        trange  = self.input.usetc(prompt,trange)        

        if '[' == trange[0]:
            trange  = trange[1:]
        if ']' == trange[-1]:
            trange  = trange[:-1]

        trange  = trange.split('->')

        tolerance2 = 0.000001*(float(trange[1])-float(trange[0]))
        mintime    = float(trange[0])-tolerance2
        maxtime    = float(trange[1])+tolerance2

        thedatas = []
        for var in thevars:
            index   = thevars.index(var)
            x       = Numeric.choose(Numeric.less_equal(times[index],mintime),
                                     (times[index],0))
            y       = Numeric.choose(Numeric.greater_equal(x,maxtime),(x,0))
            nz      = Numeric.nonzero(y)
            if index == 0:
                #store all columns except for the cycle column
                thedatas.append(var.data[nz[0]:nz[-1]+1,1:var.shape[1]])
            else:
                #store all columns except for the cycle and time columns
                thedatas.append(var.data[nz[0]:nz[-1]+1,2:var.shape[1]])

        alldata = Numeric.concatenate(thedatas,1)

        st      = '\n'
        st     += '#Probe Name               : %s \n' %(theprobe.name)
        st     += '#Probe Description        : %s \n' %(theprobe.description)
        coords  = '['
        for i in theprobe.coords:
            coords += self.coordsfmt%i
        st     += '#Probe Coordinates        : %s] \n' %(coords)

        coords  = '['
        for i in theprobe.cell['COORDS']:
            coords += self.coordsfmt%i
        st     += '#Closest Cell Index       : %i \n' %(theprobe.cell['ID'])
        st     += '#Closest Cell Coordinates : %s] \n' %(coords)

        coords  = '['
        for i in theprobe.node['COORDS']:
            coords += self.coordsfmt%i
        st     += '#Closest Node Index       : %i \n' %(theprobe.node['ID'])
        st     += '#Closest Node Coordinates : %s] \n ' %(coords)
        
        #format title appropriately for scalar, vector, tensor variables
        theshape = Numeric.shape(alldata)[1]+1

        names    = []
        for k in thevars:
            if k.shape[1]-2 > 1:
                for j in range(k.shape[1]-2):
                    names.append(str(k.nickname) + str(j+1))
            else:
                names.append(str(k.nickname))
                
        tformat  = '%+15s'
        title    = '#Time'
        st      += tformat%(title)
        tformat  = ' %+19s'
        for k in names:
            st  += tformat%(k)

        tformat  = ' %+21s \n'
        title    = 'Cycle'
        st      += tformat%(title)

        myShape   = Numeric.array(Numeric.shape(alldata),'i')
        myCycles  = thevars[0].data[nz[0]:nz[-1]+1,0:1].astype('i')
        lencycles = int(len(myCycles))

        "Accomodate Numpy/Numeric differences"
        try:    # Numpy
           myData = alldata.flat
           myData = Numeric.reshape(myData,myShape)
        except: # Numeric
           myData = Numeric.array(alldata)

        try:
            Estring     = 'asciiwriter import/writeProbe error'
            import asciiwriter
            asciiwriter.writeProbe(len(myShape), myShape, myData, myCycles, lencycles, st, fp, iformat, self.fieldfmt, 0)
        except:
            print >> self.fp
            print >> self.fp, 'Ascii C writer failed....'
            print >> self.fp, 'Writing out probe arrays in Python.  ',  
	    print >> self.fp, 'This will be slow!'
            print >> self.fp

            fp.write(st)
            format = self.fieldfmt
            s      = ''
            count  = 0
            for v in alldata:
                for j in v:
                    s    += format%j
                thiscycle = int(myCycles[count][0])
                s     += iformat%(thiscycle)
                s     += '\n'
                count += 1
            fp.write(s)
 
            
    def query(self,X):

        import sys
        try:
           import numpy.oldnumeric as Numeric
        except ImportError:
           import Numeric
        except:
           raise
        
        from writeStorageObject import writeStorageObject
        
        this_storage    = X.storage
        thisfile        = X.file

        prompt          = 'Output to "screen" or "file"? '
        oupt            = self.input.usetc(prompt, 'screen')
        if oupt == 'file':
            prompt      = 'Filename to output results to?'
            fname       = self.input.usetc(prompt, 'query.txt')
            fp          = open(fname,'w')
        else:
            fp          = self.fp
            
        count           = 0
        varlist         = '\n'
        for var in this_storage.tlist[0].vlist:
            s           = '%+25s' %(str(var.nickname))
            s          += ','
            varlist    += s
            count      += 1
            if (count%self.varformat == 0):
                varlist = varlist + '\n'
        varlist     = varlist[0:len(varlist)-1]
        prompt      = 'Query which variable? (%s)' %(varlist)
        varname     = self.input.usetc(prompt,
                                 str(this_storage.tlist[0].vlist[1].nickname))

        if (len(varname)<1):
            return # No data provided by user 

        novar               = 1
        face                = 0
        for k in this_storage.tlist[0].vlist:
            if(k.name == varname or k.nickname == varname):
                thevar       = k
                novar        = novar*0
                if k.rank == 3:
                    print >> self.fp, 'This variable is face-centred vector with %i components' % \
				      (k.shape[-1])
                    facelist = range(k.shape[-2]+1)
                    facelist = str(facelist[1:])
                    prompt = 'Specify which face you would like to query out of ' + facelist
                    face     = self.input.usetc(prompt, '1')
                    face     = int(face)
        if novar:
            print >> self.fp
            print >> self.fp, 'Your chosen variable is not in this file'
            print >> self.fp, 'Aborting this query.'
            print >> self.fp
            return
        
        regionlist     = ''
        rlist          = []
        for region in self.files_regions[X.file]:
            #check region meshspaces to ensure they are
            #consistent with the variable meshspace
            invalidregion = 1
            for mspace in region.meshspace:
                if type(mspace) == list:
                    for i in mspace:
                        if i == thevar.meshspace:
                            invalidregion = invalidregion*0
                else:
                    if mspace == thevar.meshspace:
                        invalidregion     = invalidregion*0
            if not invalidregion:
                regionlist = regionlist + ' ' + region.name + ' '
                rlist      = rlist + [str(region.name)]
        
        prompt     = 'For which region? Available regions:  [ '+regionlist+']'
        rname      = self.input.usetc(prompt, rlist[-1])

        if(len(rname)<1):
            return # No data provided by user
        else:
            for region in self.files_regions[X.file]:
                if region.name == rname:
                    thisregion = region
        try:
            tmp = thisregion
        except:
            print >> self.fp
            print >> self.fp, 'Region %s not in the list' %(rname)
            print >> self.fp, 'Will abort querying this region'
            print >> self.fp
            return

        lowcycle   = this_storage.tlist[0].cycle
        highcycle  = this_storage.tlist[-1].cycle
        
        lower      = self.input.useti('Low Cycle number?',lowcycle)
        upper      = self.input.useti('High Cycle number?',highcycle)

        idx2       = []
        cnt        = 0
        for i in this_storage.timesteps.clist[thevar.name]:
            if (i <= upper and i >= lower):
                idx2 += [cnt]
            cnt += 1

        count       = 0
        thisdformat = self.fieldfmt[0]+str(max(int(self.fieldfmt[1:3])-10,14)) \
		      + self.fieldfmt[3:]
        timeformat  = self.fieldfmt[0] + str(int(self.fieldfmt[1:3])+5) \
		      + self.fieldfmt[3:]
        
        for i in idx2:

            thistime = timeformat%(this_storage.tlist[i].time)
            s = '#' + '-' * 115 +'\n' 
            s += '#' + '%+15s %+25s %+25s' %('VAR', 'TIME', 'CYCLE')
            s += '\n'
            s += '#' + '%+15s %s %23d '\
                 %(varname, thistime, this_storage.tlist[i].cycle)
            s += '\n' + '#' + '-' * 115 + '\n'

            t  = s

            if count == 0:
                #find co-ordinates for the variables
                if thevar.meshspace == 'CELL':
                    #if the variable lives on the cell space,
                    #use CELL CENTROIDS for positions
                    posns   = thisregion.mesh.cells['CENTROIDS']
                if thevar.meshspace == 'VERTEX':
                    #if the variable lives on the vertex space,
                    #use VERTEX COORDS for positions
                    posns   = thisregion.mesh.vertices['COORDS']
                indices = thisregion.indices[thevar.meshspace][0]

            thisid = this_storage.tlist[i].id

            # now update region indices if the region
            # has a 'VAR' selector associated with it
            if 'VAR' in thisregion.selector and thisid > 0:
                thisregion.getIndices('VAR',cycleid=thisid)
                indices = thisregion.indices[thevar.meshspace][thisid]
            
            count += 1

            for k in this_storage.tlist[i].vlist:

                if (k.name == varname or k.nickname == varname):

                    this_storage.getValues([k])
                    if face == 0:
                        st = '%s %+10s %+25s %+30s\n'\
                             %('#','Index', 'Position',  'Value(s)')
                    else:
                        st = '%s %+10s %+25s %+50s %d\n'\
                             %('#','Index','Position','Value(s) at Face: ',face)

                    t     +=  st
                    myindices  = Numeric.array(indices,'i')
                    writeStorageObject(fpwatch = self.fp,
                                       debug   = self.debug
				       ).writeFieldAndRegion(face, k, myindices,
						             posns, t, fp,
                                                             coordsfmt = self.coordsfmt,
                                                             fieldfmt  = thisdformat) 

        if oupt == 'file':
            fp.close()

    def stat(self,X):

        try:
           import numpy.oldnumeric as Numeric
        except ImportError:
           import Numeric
        except:
           raise

        this_storage    = X.storage
        thisfile        = X.file

        prompt          = 'Output to "screen" or "file"? '
        oupt            = self.input.usetc(prompt, 'screen')
        if oupt == 'file':
            prompt      = 'Filename to output results to?'
            fname       = self.input.usetc(prompt, 'stat.txt')
            fp          = open(fname,'w')

        count           = 0
        varlist         = '\n'
        for var in this_storage.tlist[0].vlist:
            s           = '%+25s' %(str(var.nickname))
            s           = s + ','
            varlist     = varlist + s
            count       = count + 1
            if (count%self.varformat == 0):
                varlist = varlist + '\n'
                
        varlist     = varlist[0:len(varlist)-1]
        prompt      = 'Query statistics for which variable? (%s)' %(varlist)
        varname     = self.input.usetc(prompt,
                                 str(this_storage.tlist[0].vlist[1].nickname))
        
        if (len(varname)<1):
            return # No data provided by user 

        novar                = 1
        for k in this_storage.tlist[0].vlist:
            if(k.name == varname or k.nickname == varname):
                thevar       = k
                novar        = novar*0
                if k.rank > 1:
                    dims     = k.shape[-1]
                else:
                    dims     = 1
                dimlist      = range(dims)
                if k.rank == 3:
                    print >> self.fp, 'This variable is face-centred vector with %i components' % \
				      (k.shape[-1])
                    facelist = range(k.shape[-2]+1)
                    facelist = str(facelist[1:])
                    prompt = 'Specify which face you would like to query statistics out of ' + facelist
                    face     = self.input.usetc(prompt, '1')
                    face     = int(face)
        if novar:
            print >> self.fp
            print >> self.fp, 'Your chosen variable is not in this file'
            print >> self.fp, 'Aborting this query.'
            print >> self.fp
            return

        regionlist        = ''
        rlist             = []
        for region in self.files_regions[X.file]:
            #check region meshspaces to ensure
            #they are consistent with the variable meshspace
            invalidregion = 1
            for mspace in region.meshspace:
                if type(mspace) == list:
                    for i in mspace:
                        if i == thevar.meshspace:
                            invalidregion = invalidregion*0
                else:
                    if mspace == thevar.meshspace:
                        invalidregion     = invalidregion*0
            if not invalidregion:
                regionlist = regionlist + ' ' + region.name + ' '
                rlist      = rlist + [str(region.name)]
        
        prompt     = 'For which region? Available regions:  [ '+regionlist+']'
        rname      = self.input.usetc(prompt, rlist[-1])

        if(len(rname)<1):
            return # No data provided by user
        else:
            for region in self.files_regions[X.file]:
                if region.name == rname:
                    thisregion = region
        try:
            tmp = thisregion
        except:
            print >> self.fp
            print >> self.fp, 'Region %s not in the list' %(rname)
            print >> self.fp, 'Will abort finding statistics on this region'
            print >> self.fp
            return

        lowcycle   = this_storage.tlist[0].cycle
        highcycle  = this_storage.tlist[-1].cycle

        lower      = self.input.useti('Low Cycle number?',lowcycle)
        upper      = self.input.useti('High Cycle number?',highcycle)

        idx2       = []
        cnt        = 0
        for i in this_storage.timesteps.clist[thevar.name]:
            if (i <= upper and i >= lower):
                idx2 += [cnt]
            cnt += 1

        s     = '\n'
        s2    = '\n'
        count = 0

        timeformat  = self.fieldfmt[0] + str(int(self.fieldfmt[1:3])+5) \
		      + self.fieldfmt[3:]

        for i in idx2:

            s += '-' * 115 +'\n'
            s += '%+10s %+20s %+25s' %('VAR', 'TIME', 'CYCLE')
            s += '\n'
            thistime = timeformat%this_storage.tlist[i].time
            s += '%+10s %s %20d '\
                 %(varname, thistime, this_storage.tlist[i].cycle)
            s += '\n' + '-' * 115 + '\n'

            if count == 0:
                #find co-ordinates for the variables only do this once
                if thevar.meshspace == 'CELL':
                    #if the variable lives on the cell space,
                    #use CELL CENTROIDS for positions
                     posns   = thisregion.mesh.cells['CENTROIDS']
                if thevar.meshspace == 'VERTEX':
                    #if the variable lives on the vertex space,
                    #use VERTEX COORDS for positions
                    posns   = thisregion.mesh.vertices['COORDS']
                indices     = thisregion.indices[thevar.meshspace][0]

            thisid = this_storage.tlist[i].id

            #now update region indices if the region
            #has a 'VAR' selector associated with it
            if 'VAR' in thisregion.selector and thisid > 0:
                thisregion.getIndices('VAR',cycleid=thisid)
                indices = thisregion.indices[thevar.meshspace][thisid]

            count      += 1
            thiscformat = '    [%s,%s,%s   ]'%(self.coordsfmt,self.coordsfmt,
			                       self.coordsfmt) 
            thisdformat = self.fieldfmt[0] \
			  + str(max(int(self.fieldfmt[1:3])-10,10)) \
			  + self.fieldfmt[3:]

            for k in this_storage.tlist[i].vlist:
                
                if (k.name == varname or k.nickname == varname):
                    this_storage.getValues([k])
                    subset      = Numeric.array([0],k.type)
                    subset      = Numeric.resize(subset,[len(indices)])
                    map         = Numeric.array([0],'i')
                    map         = Numeric.resize(map,[len(indices)])

                    if k.rank == 1:
                        dim   = 1
                        cnt   = 0
                        for j in indices:
                            subset[cnt] = k.data[j-1]
                            map[cnt]    = j
                            cnt        += 1
                        
                        mini  = Numeric.argmin(subset)
                        maxi  = Numeric.argmax(subset)
                        vmin  = subset[mini]
                        vmax  = subset[maxi]
                        size  = len(subset)
                        ident = Numeric.array([1],'d')
                        ident = Numeric.resize(ident,[size])
                        ave   = Numeric.sum(subset)/size
                        stdev = (subset-ave*ident)**2
                        stdev = (Numeric.sum(stdev)/size)**(0.5)

                        xpos   = posns[mini-1][0]
                        ypos   = posns[mini-1][1]
                        zpos   = posns[mini-1][2]
                        pos    = thiscformat %(xpos,ypos,zpos)
                        miniat = '%10i %+15s ' %(map[mini], pos)
                        
                        xpos   = posns[maxi-1][0]
                        ypos   = posns[maxi-1][1]
                        zpos   = posns[maxi-1][2]
                        pos    = thiscformat %(xpos,ypos,zpos)
                        maxiat = '%10i %+15s ' %(map[maxi], pos)

                        s2  = 'COMP                  : %s \n'%('%10i'%(dim))
                        s2 += 'AVE                   : %s \n'%(thisdformat%(ave))
                        s2 += 'MIN                   : %s \n'%(thisdformat%(vmin))
                        s2 += 'MIN at (Index,[Posn]) : %s \n'%(miniat)
                        s2 += 'MAX                   : %s \n'%(thisdformat%(vmax))
                        s2 += 'MAX at (Index,[Posn]) : %s \n'%(maxiat)
                        s2 += 'STDEV                 : %s \n'%(thisdformat%(stdev))

                        s  += s2
                        
                    else: #Dealing with vectors, require further dimensions
                    
                        s2 = '\n'
                        for dim2 in dimlist:
                            dim  = dim2 + 1
                        
                            if k.rank == 2:

                                cnt   = 0
                                for j in indices:
                                    subset[cnt] = k.data[j-1,dim-1]
                                    map[cnt]    = j
                                    cnt        += 1

                                mini  = Numeric.argmin(subset)
                                maxi  = Numeric.argmax(subset)
                                vmin  = subset[mini]
                                vmax  = subset[maxi]
                                size  = len(indices)
                                ident = Numeric.array([1],'d')
                                ident = Numeric.resize(ident,[size])
                                ave   = Numeric.sum(subset)/size
                                stdev = (subset-ave*ident)**2
                                stdev = (Numeric.sum(stdev)/size)**(0.5)
                                
                            if k.rank == 3:

                                cnt   = 0
                                for j in indices:
                                    subset[cnt] = k.data[j-1,face-1,dim-1]
                                    map[cnt]    = j
                                    cnt        += 1

                                mini  = Numeric.argmin(subset)
                                maxi  = Numeric.argmax(subset)
                                vmin  = subset[mini]
                                vmax  = subset[maxi]
                                size  = len(subset)
                                ident = Numeric.array([1],'d')
                                ident = Numeric.resize(ident,[size])
                                ave   = Numeric.sum(subset)/size
                                stdev = (subset-ave*ident)**2
                                stdev = (Numeric.sum(stdev)/size)**(0.5)

                            xpos   = posns[mini-1][0]
                            ypos   = posns[mini-1][1]
                            zpos   = posns[mini-1][2]
                            pos    = thiscformat %(xpos,ypos,zpos)
                            miniat = '%10i %+15s ' %(map[mini], pos)
                        
                            xpos   = posns[maxi-1][0]
                            ypos   = posns[maxi-1][1]
                            zpos   = posns[maxi-1][2]
                            pos    = thiscformat %(xpos,ypos,zpos)
                            maxiat = '%10i %+15s ' %(map[maxi], pos)

                            s2 += 'COMP                  : %s \n'%('%10i'%(dim))
                            s2 += 'AVE                   : %s \n'%(thisdformat%(ave))
                            s2 += 'MIN                   : %s \n'%(thisdformat%(vmin))
                            s2 += 'MIN at (Index,[Posn]) : %s \n'%(miniat)
                            s2 += 'MAX                   : %s \n'%(thisdformat%(vmax))
                            s2 += 'MAX at (Index,[Posn]) : %s \n'%(maxiat)
                            s2 += 'STDEV                 : %s \n'%(thisdformat%(stdev))
                            s2 += '\n'

                        s += s2

                    s += '\n'
                    
                
        if oupt == 'file':
            print >> fp, s
        else:
            print >> self.fp, s
                                                                                        
        if oupt == 'file':fp.close()

    def timeseries(self,X):

        varformat    = 3
        this_storage = X.storage

        count           = 0
        varlist         = '\n'
        vlist           = []
        for var in this_storage.tlist[0].vlist:
            if var.rank > 0:
                s           = '%+30s' %(str(var.nickname))
                s          += ','
                varlist     = varlist + s
                vlist.append(var.nickname)
                count      += 1
                if (count%varformat == 0):
                    varlist = varlist + '\n'
        varlist     = varlist[0:len(varlist)-1]
        prompt      = 'Time series for which variable? (%s)' %(varlist)
        varname     = self.input.usetc(prompt,
                                 str(this_storage.tlist[0].vlist[1].nickname))

        if (len(varname)<1):
            print >> self.fp
            print >> self.fp, 'No variable specified..try again'
            print >> self.fp
            return

        if not (varname in vlist):
            print >> self.fp
            print >> self.fp, 'This variable is not available in the XML output... try again'
            print >> self.fp
            return

        novar               = 1
        face                = 0
        for k in this_storage.tlist[0].vlist:
            if(k.name == varname or k.nickname == varname):
                thevar       = k
                novar        = novar*0
                if k.rank == 3:
                    print >> self.fp, 'This variable is face-centred vector with %i components' %(k.shape[-1])
                    facelist = range(k.shape[-2]+1)
                    facelist = str(facelist[1:])
                    prompt   = 'Specify which face you would like to analyze out of ' + facelist
                    face     = self.input.usetc(prompt, '1')
                    face     = int(face)

        #Check the variable lives on a meshspace.
        # - if so inquire as to which index to analyze
        # - otherwise just print out the variable

        nomeshspace = 0
        if thevar.meshspace == 'CELL':
            for meshspace in this_storage.mlist[0].mslist:
                if meshspace.name == 'CELL':
                    if meshspace.size != thevar.shape[0]:
                        nomeshspace = 1

        if not nomeshspace:
            for meshspace in this_storage.mlist[0].mslist:
                if thevar.meshspace == meshspace.name:
                    ninds     = int(meshspace.size)
                    nindsh    = int(0.5*ninds)
                    prompt    = 'Index to be analyzed in the time series [1 -> %s]' %(str(ninds))
                    index     = self.input.useti(prompt,str(nindsh))
                    #check index is within the correct bounds
                    if index > ninds or index < 1:
                        print >> self.fp
                        print >> self.fp, 'Index not within the correct bounds'
                        print >> self.fp, 'Aborting this timeseries choice..'
                        return
                        
        else:
            index = 'Not relevant'

        #Now inquire about which cycles to print out...

        lowcycle   = this_storage.tlist[0].cycle
        highcycle  = this_storage.tlist[-1].cycle
        
        lower      = self.input.useti('Low Cycle number?',lowcycle)
        upper      = self.input.useti('High Cycle number?',highcycle)

        idx2       = []
        cnt        = 0
        for i in this_storage.timesteps.clist[thevar.name]:
            if (i <= upper and i >= lower):
                idx2 = idx2 + [cnt]
            cnt += 1

        #Header..

        s = '#' + '-' * 115 +'\n'
        s +='#' + '%+10s %+20s %+30s (Index = %s)' \
             %('CYCLE', 'TIME', varname, str(index))
        s += '\n' + '#' + '-' * 115 + '\n'

        for i in idx2:
        
            for var in this_storage.tlist[i].vlist:
                if (var.name == varname or var.nickname == varname):
                    thevar = var

            this_storage.getValues([thevar])

            tmp = ''
            if thevar.rank == 1:
                if nomeshspace :
                    #e.g sensitivity variables or em variables
                    # - i.e variables that don't live on a meshspace
                    for j in thevar.data:
                        tmp += self.fieldfmt %(j)
                else:
                    tmp    = self.fieldfmt %(thevar.data[index-1])
            if thevar.rank == 2:
                for j in thevar.data[index-1]:
                    tmp += self.fieldfmt %(j)
            if thevar.rank == 3:
                cnt = 0
                for j in thevar.data[index-1]:
                    cnt += 1
                    if cnt == face:
                        for f in j:
                            tmp += self.fieldfmt %(f)

            thistime = self.fieldfmt %(this_storage.tlist[i].time)
            s += '%10d %s %+30s '\
                 %(this_storage.tlist[i].cycle, thistime, tmp)
            s += '\n'

        #Now inquire about the name of the txt output file
        prompt          = 'Output to "screen" or "file"? '
        oupt            = self.input.usetc(prompt, 'screen')
        if oupt == 'file':        
            prompt      = 'Filename to output results to?'
            if nomeshspace:
                defname = '%s_%s_to_%s.txt' \
                          %(thevar.nickname,str(lower),str(upper))
            else:
                defname = '%s_%s_%s_to_%s.txt' \
                          %(thevar.nickname,str(index),str(lower),str(upper))
            fname       = self.input.usetc(prompt, defname)
            fp          = open(fname,'w')
            print >> fp, s
            fp.close()
        else:
            print >> self.fp, s


    def write(self,X,choice):

        import os,string

        try:
           import numpy.oldnumeric as Numeric
        except ImportError:
           import Numeric
        except:
           raise        

        this_storage = X.storage
        thisfile     = X.file

        #now check which meshes to write too

        m      = []   
        meshes = X.storage.mlist
        if(len(meshes)):
            for mesh in meshes:
                m = m + [str(mesh.name)]

        if(len(m)): #we have meshes in the file
            meshlist = ''
            for mesh in meshes:
                meshlist = meshlist + str(mesh.name) + ' '
            prompt = 'Available meshes: '+ meshlist + ' Which mesh do you want to output? '
            c      = self.input.usetc(prompt, m[-1])

            try:
                for mesh in meshes:
                    if c == mesh.name: thismesh = mesh
            except:
                print >> self.fp
                print >> self.fp, 'Sorry, invalid entry.  Must pick one of: ',meshlist
                print >> self.fp, 'Aborting graphics writer.'
                print >> self.fp
                return
        else:
            print >> self.fp, 'No meshes in file.'
            print >> self.fp, 'A graphics file will not be produced. \n'
            return

        #now get which cycles you want to output

        idx2 = self.__chooseSteps(X)

        #get choice of ascii or binary

        out     = ['ascii', 'binary']
        if choice == 'tecplot':
            out = ['ascii']

        prompt = 'Available output formats : ' + str(out) + ' Which format do you want to output in? '
        co     = self.input.usetc(prompt, out[-1])
        try:
            for output in out:
                if co == output: outputfmt = output
        except:
            print >> self.fp
            print >> self.fp, 'Sorry, invalid entry.  Must pick one of: ',out
            print >> self.fp, 'Aborting graphics writer.'
            print >> self.fp
            return
        iend         = string.find(X.file,".TBrook.xml")
        thefile      = X.file[0:iend]

        if choice == 'ensight':
            iend         = string.find(X.file,".TBrook.xml")
            ibeg         = string.rfind(X.file,"/") + 1
            fileout      = X.file[0:iend]
            fileroot     = fileout
            if co == 'binary':
                print >> self.fp, '*'*60
                print >> self.fp, 'EnSight binary dumps written in EnSight gold format.'
                print >> self.fp, '*'*60
        else:
            # generate output file name
            prompt   = 'Output graphics file name:'
            if len(idx2) == 0:
                fileout = '%s.%smesh'%(thefile,choice)
            elif len(idx2) == 1 :
                if choice == 'vtk':
                    fileout = '%s.%d.%s'%(thefile,idx2[0],choice)
                else:
                    fileout = '%s.%s.%d'%(thefile,choice,idx2[0])
            else:
                fileout  = '%s'%(thefile)
                prompt   = 'Root of graphics output file names:'

            fileout  = self.input.usetc(prompt, fileout)
            fileroot = fileout

        # Get the partition information
        flags = {}
        for r in X.regions:
            if r.indices.has_key('CELL'):
                flags[r.name]     = r.indices['CELL'][0]
                if 'VAR' in r.selector and thisid > 0:
                    r.getIndices('VAR',cycleid=thisid)
                    flags[r.name] = r.indices['CELL'][thisid]
            else:
                flags[r.name] = []

        if len(idx2) == 0:
            #write only the mesh to the graphics file
            if choice == 'ensight':
                ENSIGHTwriteStorageObject(fileout, outputfmt,
                                          thismesh, 0.0, 0,
                                          dformat = self.fieldfmt,
                                          bFlags  = flags,
                                          fpwatch = self.fp,
                                          debug   = self.debug)
            elif choice == 'gmv':
                GMVwriteStorageObject(fileout, outputfmt,
                                      thismesh, 0.0, 0,
                                      dformat = self.fieldfmt,
                                      bFlags  = flags,
                                      fpwatch = self.fp,
                                      debug   = self.debug) 
            elif choice == 'tecplot':
                TECPLOTwriteStorageObject(fileout, outputfmt,
                                          thismesh, 0.0, 0,
                                          dformat = self.fieldfmt,
                                          bFlags  = flags,
                                          fpwatch = self.fp,
                                          debug   = self.debug)
            elif choice == 'vtk':
                VTKwriteStorageObject(fileout, outputfmt,
                                      thismesh, 0.0, 0,
                                      dformat = self.fieldfmt,
                                      bFlags  = flags,
                                      fpwatch = self.fp,
                                      debug   = self.debug)
            return

        #otherwise, get user choice of variables to be plotted..
        var     = []
        varlist = '\n'
        count   = 0
        for var in X.storage.tlist[0].vlist:
            if len(var.mesh) :
                writevar = writeStorageObject().writeThisVar(var)
                if writevar:
                    if var.rank <=2 and var.rank > 0:
                        s        = '%+25s' %(str(var.nickname))
                        s       += ','
                        varlist += s
                        count   += 1
                        if (count%self.varformat == 0):
                            varlist = varlist + '\n'
        varlist += '%+25s, %+25s' %('all','truchas')
        prompt   = 'List the variables to be written. Enclose the list in double quotes with commas between variables (i.e "T,P") %s ' %(varlist)
        varnames = self.input.usetc(prompt, 'truchas')
        seq_no   = 0

        for idx in idx2:
            # get the values for this timestep
            thisid   = X.storage.tlist[idx].id
            
            if (choice == 'gmv') and (len(idx2) > 1): 
                fileout = '%s.%s.%06d'%(fileroot, choice, idx)
            elif (choice == 'vtk') and (len(idx2) > 1): 
                fileout = '%s.%d.%s'%(fileroot, seq_no, choice)
            elif (choice == 'tecplot') and (len(idx2) > 1): 
                fileout = '%s.dat'%(fileroot)
            
            # Get the region information
            for r in X.regions:
                if r.indices.has_key('CELL'):
                    flags[r.name]     = r.indices['CELL'][0]
                    if 'VAR' in r.selector and thisid > 0:
                        r.getIndices('VAR',cycleid=thisid)
                        flags[r.name] = r.indices['CELL'][thisid]
                else:
                    flags[r.name] = []
            
            #now obtain variables
            thisvlist    = [] 
            if varnames == 'all':
                thisvlist = X.storage.tlist[idx].vlist
            elif varnames == 'truchas':
                for k in X.storage.tlist[idx].vlist:
                    if     ('_SENS'      not in k.name) \
                       and ('_sens'      not in k.name) :
                        thisvlist.append(k)
            else:
                newlist = string.split(varnames,',')
                for nic in newlist:
                    if nic[0] == ' ':
                        nic  = nic[1:len(nic)]
                    for k in X.storage.tlist[idx].vlist:
                        if k.nickname == nic:
                            thisvlist.append(k)

            print >> self.fp, 'For cycle number: %i' %(thisid)

            X.storage.getValues(thisvlist)
            if choice == 'ensight':
                ENSIGHTwriteStorageObject(fileout,   outputfmt,
                                          thismesh,  X.storage.tlist[idx].time, 
					  X.storage.tlist[idx].cycle,
                                          thisvlist, seq_no,
                                          dformat = self.fieldfmt,
                                          bFlags  = flags,
                                          fpwatch = self.fp,
                                          debug   = self.debug)
            elif choice == 'gmv':
                GMVwriteStorageObject(fileout,   outputfmt,
                                      thismesh,  X.storage.tlist[idx].time, 
				      X.storage.tlist[idx].cycle,
                                      thisvlist, seq_no,
                                      dformat = self.fieldfmt,
                                      bFlags  = flags,
                                      fpwatch = self.fp,
                                      debug   = self.debug)
            elif choice == 'tecplot':
                TECPLOTwriteStorageObject(fileout,   outputfmt,
                                          thismesh,  X.storage.tlist[idx].time,
					  X.storage.tlist[idx].cycle,
                                          thisvlist, seq_no,
                                          dformat = self.fieldfmt,
                                          bFlags  = flags,
                                          fpwatch = self.fp,
                                          debug   = self.debug)            
            elif choice == 'vtk':
                VTKwriteStorageObject(fileout,   outputfmt,
                                      thismesh,  X.storage.tlist[idx].time, 
				      X.storage.tlist[idx].cycle,
                                      thisvlist, seq_no,
                                      dformat = self.fieldfmt,
                                      bFlags  = flags,
                                      fpwatch = self.fp,
                                      debug   = self.debug)            
            X.storage.tlist[idx].tossData()
            seq_no += 1


    def restart(self,Xs):

        import os,sys,string
        
        prompt         = 'Create a standard or a mapped restart?'
        type           = self.input.usetc(prompt,'standard')
        type.lower()
        if type == 'standard':
            mapping    = 0
            print >> self.fp
            print >> self.fp, '     Creating standard restart file(s)....'
        elif type == 'mapped':
            mapping         = 1
            self.mapcount  += 1
            print >> self.fp
            print >> self.fp, '   Creating a restart file by mapping to a user specified mesh....'
            
        meshchs      = []
        ignoredfiles = []
        for X in Xs:
            m      = []   
            meshes = X.storage.mlist
            if (len(meshes)):
                meshlist = ''
                for mesh in meshes:
                    m        += [str(mesh.name)]
                    meshlist += str(mesh.name) + ' '
                print >> self.fp
                print >> self.fp, '     For simulation %s:' %(X.file)
                print >> self.fp, '     Available meshes are %s ' %(meshlist)
                print >> self.fp
                prompt = 'Create restart from which mesh? ' 
                meshchs.append(self.input.usetc(prompt, m[-1]))
            else:
                print >> self.fp
                print >> self.fp, '     No meshes in the simulation: %s' %(X.file)
                print >> self.fp, '     Will not create a restart for this simulation'
                print >> self.fp
                ignoredfiles.append(Xs.index(X))

        T  = []
        if len(ignoredfiles)> 0 and len(Xs) > 0:
            for ind in ignoredfiles: T.append(Xs[ind])
            if len(T) == len(Xs): return
            for val in T:
                if val in Xs: del Xs[Xs.index(val)]

        #now get which cycles you want to output for all restarted simulations
        ids          = {}
        cycles       = {}
        cnt          = 0
        ignoredfiles = []

        for X in Xs:
            ignorefile           = 0
            ids[cnt], ignorefile = self.__chooseRestartSteps(X, type)
            if ignorefile:
                ignoredfiles.append(Xs.index(X))
            else:
                cnt += 1
        T  = []
        if len(ignoredfiles) > 0 and len(Xs) > 0:
            for ind in ignoredfiles: T.append(Xs[ind])
            if len(T) == len(Xs): return
            for val in T:
                if val in Xs: del Xs[Xs.index(val)]

        these_storages = []
        these_files    = []
        these_meshes   = []
        
        for X in Xs:
            these_storages.append(X.storage)
            these_files.append(X.file)
            for mesh in X.storage.mlist:
                for meshch in meshchs:
                    if meshch == mesh.name:
                        these_meshes.append(mesh)
                            
        if mapping:
            try:
                #check for any existing exodus files
                exof    =  self.def_exofile
                for thisfile in os.listdir(os.getcwd()):
                    if thisfile[-4:] == '.exo':
                        exof = thisfile
                print >> self.fp
                print >> self.fp, ' '*10,'Restarting on a different Exodus II mesh. '
                print >> self.fp, ' '*10,'Rules to map variables to the new mesh are provided in:'
                print >> self.fp
                print >> self.fp, ' '*13,'your_truchasdir/tools/PythonPackages/TBrookParser/MODutils/MAPutils/defaultmaps.txt'
                print >> self.fp
                print >> self.fp, ' '*10,'You may change these rules to suit your application.'
                print >> self.fp, ' '*10,'Now specify the Exodus II input file containing the target mesh to restart on.. '
                print >> self.fp
                stepids = []
                for i in range(len(these_files)):
                    stepids.append(ids[i][0])

                Y  = getFileObject(self.input, self.def_exofile,
                                   these_files, self.files_storages,
                                   stepids, self.mapcount,
                                   fpwatch = self.fp,
                                   debug   = self.debug)
                """
                Since we have chosen to remap we must now alter the storage 
                 restart features so that the features:
                      1.temperature
                      2.solid_mechanics
                      3.fluid_flow
                      4.joule_heat
                 are ignored in the restart
                (i.e will not be read in in the newly mapped restart file) 
                """
                feats  = Y.storage.specs.feats

                for feat in feats:
                    if 'temperature' in feat:
                        del feats[feats.index(feat)]
                for feat in feats:
                    if 'fluid_flow' in feat:
                        del feats[feats.index(feat)]
                for feat in feats:
                    if 'joule_heat' in feat:
                        del feats[feats.index(feat)]
                for feat in feats:
                    if 'solid_mechanics' in feat:
                        del feats[feats.index(feat)]
                        
                Y.storage.specs.feats  = feats
                Y.storage.specs.nfeats = len(feats)

                these_files       = []
                these_storages    = []
                these_meshes      = []
                these_files.append(Y.file)
                these_storages.append(Y.storage)
                these_meshes.append(Y.storage.mlist[0])
                # if we have chosen to map to a new mesh,
                # adjust the cycles list and id list
                cnt          = 0
                cycles[cnt]  = Y.storage.timesteps.cyclelist
                ids[cnt]     = []
                ids[cnt].append(Y.storage.timesteps.idlist[-1])
                if len(cycles[cnt]) > 1 :
                    print >> self.fp
                    print >> self.fp, '     In this postprocessing event you have performed '
                    print >> self.fp, '     %i mappings (i.e cycles %s) with these mesh choices.' \
                          %(len(cycles[cnt]), str(cycles[cnt]))
                    print >> self.fp
                    prompt         = 'Which of these mappings would you like to use to create the restart?'
                    cycles[cnt][0] = self.input.useti(prompt,
                                                      cycles[cnt][-1]) 
                    ids[cnt]       = [cycles[cnt][0]]
                    
                #now save this file object in the post-processor's storage structure
                if not (Y.file in self.files):
                    self.files.append(Y.file)
        
                self.files_storages[Y.file] = [Y.storage]
                self.files_regions[Y.file]  =  Y.regions
            except:
                print >> self.fp
                print >> self.fp, 'Error in obtaining new restart file from %s '%(exofile)
                print >> self.fp, 'Aborting restart writer.'
                print >> self.fp
                print >> self.fp, 'Please ensure your original file:'
                print >> self.fp, '1.is an XML file'
                print >> self.fp, 'Please ensure your new mesh: '
                print >> self.fp, '1. Is a non-degenerate HEX Mesh'
                print >> self.fp, '2. Overlaps with the Truchas source mesh'
                print >> self.fp, '3. Is in Exodus II format'
                print >> self.fp
                return

        #get choice of ascii or binary
        print >> self.fp
        prompt = 'Available output formats : [binary, binaryswap]. Which format do you want to output in? '
        out    = ['ascii', 'binaryswap', 'binary']
        co     = self.input.usetc(prompt, out[-1])

        try:
            for output in out:
                if co == output: outputfmt = output
        except:
            print >> self.fp
            print >> self.fp, '     Sorry, invalid entry.  Must pick one of: ',out
            print >> self.fp, '     Aborting restart writer'
            print >> self.fp
            return

        # generate restart files now
        print >> self.fp
        print >> self.fp, '     Creating all restart files now..'
        
        for this_file in these_files:
            cnt     = these_files.index(this_file)
            steps   = these_storages[cnt].tlist            
            print >> self.fp, '     For simulation %s:' %(this_file)
            print >> self.fp
            iend    = string.find(this_file,".TBrook.xml")
            thefile = this_file[0:iend]
            
            if(len(ids[cnt]) == 1):
                fileout = '%s.restart.%d'%(thefile, steps[ids[cnt][0]].cycle)
                prompt  = 'Output file name:'
            else:
                fileout = '%s'%(thefile)
                prompt  = 'Root of restart output file names:'
         
            fileout  = self.input.usetc(prompt, fileout)
            fileroot = fileout

            for idx in ids[cnt]:
                # get the values for this timestep
                these_storages[cnt].getValues(steps[idx].vlist)
                if(len(ids[cnt]) > 1): 
                    fileout = '%s.restart.%d'%(fileroot, steps[idx].cycle)
                from RESTARTwriteStorageObject import RESTARTwriteStorageObject
                RESTARTwriteStorageObject(fileout, outputfmt, these_storages[cnt],
                                          these_meshes[cnt], idx,
                                          steps[idx].time, steps[idx].vlist,
                                          dformat = self.fieldfmt,
                                          fpwatch = self.fp,
                                          debug   = self.debug)

    def createregion(self):

        if len(self.files):
            thisfile = self.input.usetc('From which output file?',
                                        self.files[-1])
        else:
            print >> self.fp
            print >> self.fp, 'Try loading a file first before defining regions'
            print >> self.fp
            return
        
        if self.files_regions.has_key(thisfile):
            curregions = self.files_regions[thisfile]
            name     = 'R_%02d'%len(curregions)
            name     = self.input.usetc('Name for your region:',name)
        else:
            print >> self.fp
            print >> self.fp, 'This output file has not been loaded '
            print >> self.fp, 'Please load the file and then try creating regions from it'
            print >> self.fp, 'Aborting region creator'
            print >> self.fp
            return

        this_storage = self.files_storages[thisfile][0]
        meshes       = this_storage.mlist

        if (len(meshes)):
            meshlist   = ''
            m          = []
            for mesh in meshes:
                meshlist = meshlist + str(mesh.name) + ' '
                m        = m + [str(mesh.name)]
            prompt       = 'Available meshes: '+ meshlist + ' Which mesh do you want to select your region from? '
            meshchoice   = self.input.usetc(prompt, m[-1])
        else:
            print >> self.fp
            print >> self.fp, 'No meshes in file so aborting region creator'
            print >> self.fp
            return

        try:
            for mesh in meshes:
                if meshchoice == mesh.name:
                    thismesh = mesh
            x = thismesh.name
        except:
            print >> self.fp
            print >> self.fp, 'Sorry, invalid entry.  Must pick one of: ',meshlist
            print >> self.fp, 'Aborting this region creator.'
            print >> self.fp
            return

        prompt  = 'Select from a combination of index ID (ID), variable (VAR), spatial (S) and mesh block (MB)'
        selector = self.input.usetc(prompt, 'ID,VAR')
        selector = selector.split(',')
        
        Y        = getRegionObject(self.input,thisfile,this_storage,
                                   thismesh,name,selector)

        if len(Y.name):
            self.files_regions[thisfile].append(Y)

    def deleteregions(self):

        #delete all regions except default regions
        for thisfile in self.files:
            tmplist = []
            for region in self.files_regions[thisfile]:
                if 'allMESH' in region.describe:
                    tmplist.append(region)
            self.files_regions[thisfile] = tmplist

    def help(self,commands):

        print >> self.fp, '  Commands are: \n\n    ',commands.keys()
        print >> self.fp, commands['description']


    def __chooseSteps(self, X):
        "get timesteps to write vis files"

        import string

        steps = X.storage.tlist
        c     = []
        if(len(steps)):
            for step in steps:
                c += [step.cycle]

        if(len(c)):
            print >> self.fp
            print >> self.fp, '     Available cycles: ' + str(c)
            print >> self.fp

            crange = '[%s->%s]' %(str(steps[0].cycle),str(steps[-1].cycle))
            prompt = ' Which cycle(s) do you want to output? Use -1 for all or %s for range. ' %(crange)
            i      = self.input.usetc(prompt, str(c[-1]))

            try:
                if '->' in i:
                    if '[' == i[0]:  i = i[1:]
                    if ']' == i[-1]: i = i[:-1]
                    i    = i.split('->')
                    minc = int(i[0])
                    maxc = int(i[1])
                    idx2 = range(c.index(minc),c.index(maxc)+1)
                else:
                    i    = int(i)
                    if (i < 0): idx2 = range(len(c))
                    else: idx2 = [c.index(i)]            
            except:
                print >> self.fp
                print >> self.fp, '      Sorry, invalid entry.  Must pick one of: ',[-1]+c
                print >> self.fp, '      Viz file will not be written for simulation: '
                print >> self.fp, '      %s' %(X.file)
                print >> self.fp
                return
        else:
            print >> self.fp
            print >> self.fp, '      No timesteps in file %s.' %(X.file)
            print >> self.fp, '      Ensure XML_Output_Dt_Multiplier > 0 in your Truchas input file.'
            print >> self.fp, '      The viz file will only contain the mesh \n'
            print >> self.fp
            idx2 = []

        return idx2



    def __chooseRestartSteps(self, X, type):
        "get timesteps to write restart files"

        ignorefile = 0
        steps      = X.storage.tlist
        c          = []
        if(len(steps)):
            for step in steps:
                c += [step.cycle]

        if(len(c)):
            print >> self.fp
            print >> self.fp, '     Available cycles: ' + str(c)
            print >> self.fp

            if type == 'mapped':
                prompt = 'Which specific cycle do you want to obtain data from?'
                i      =  self.input.useti(prompt, c[-1])
                if (i < 0):
                    print >> self.fp
                    print >> self.fp, '      No! You must choose a specific cycle that houses the data'
                    print >> self.fp, '      which will be mapped to the new mesh. '
                    print >> self.fp, '      For this event, cycle %i will be chosen. ' %(c[-1])
                    print >> self.fp
                    i = -1
                idx2 = [c.index(i)]
            else:
                crange = '[%s->%s]' %(str(steps[0].cycle),str(steps[-1].cycle))
                prompt = ' Which cycle(s) do you want to restart from? Use -1 for all or %s for range. ' %(crange)
                i      = self.input.usetc(prompt, str(c[-1]))
                try:
                    if '->' in i:
                        if '[' == i[0]:  i = i[1:]
                        if ']' == i[-1]: i = i[:-1]
                        i    = i.split('->')
                        minc = int(i[0])
                        maxc = int(i[1])
                        idx2 = range(c.index(minc),c.index(maxc)+1)
                    else:
                        i    = int(i)
                        if (i < 0): idx2 = range(len(c))
                        else: idx2 = [c.index(i)]            
                except:
                    print >> self.fp
                    print >> self.fp, '      Sorry, invalid entry.  Must pick one of: ',[-1]+c
                    print >> self.fp, '      A restart file will not be written for simulation: '
                    print >> self.fp, '      %s' %(X.file)
                    print >> self.fp
                    ignorefile = 1
        else:
            print >> self.fp
            print >> self.fp, '      No timesteps in file %s.' %(X.file)
            print >> self.fp, '      Ensure XML_Output_Dt_Multiplier > 0 in your Truchas input file.'
            print >> self.fp, '      A restart file will not be written for simulation: '
            print >> self.fp, '      %s' %(X.file)
            print >> self.fp
            idx2       = []
            ignorefile = 1

        return idx2, ignorefile
    
         
    def str(self):

        print >> self.fp, 'loaded files'
        print >> self.fp, self.files

        print >> self.fp, 'storages'
        for file in self.files:
            for storage in self.files_storages[file]:
                storage.str()

        print >> self.fp, 'all regions'
        for file in self.files:
            for region in self.files_regions[file]:
                region.str()



if __name__=='__main__':
    
    'for testing this component'
    import sys

    dfltdir = '../../scripts/test_TBrookParse/samples/'
    def_file    = 'map_output/map.TBrook.xml'
    def_exofile = 'mesh-a_f.exo'

    from PYTHONutils import uTestOpts
    opts = uTestOpts('fed',
                     defaults={'f': dfltdir + def_file,    # XML file
                               'e': dfltdir + def_exofile},# EXODUS file
                     dir     = dfltdir)
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    fpwatch    = sys.stdout

    try:
        inBuffer    = ''
        input       = usubs.input(inBuffer)
        ppobject    = getPostProcessorObject(opt.f,opt.e,input,fp=fpwatch,debug=opt.d)
        fileObj     = ppobject.load(opt.f)
        #ppobject.defaults()
        #ppobject.write(fileObj,'gmv')
        #ppobject.write(fileObj,'vtk')
        #ppobject.write(fileObj,'ensight')
        ppobject.write(fileObj,'tecplot')
        ppobject.restart([fileObj])
        #ppobject.probe(fileObj)
        #ppobject.query(fileObj)
        #ppobject.stat(fileObj)
        #ppobject.timeseries(fileObj)
        #ppobject.write(fileObj,'vtk')
        #ppobject.restart([fileObj])
        #ppobject.write(fileObj,'ensight')
        #ppobject.write(fileObj,'tecplot')
        #ppobject.write(fileObj,'vtk')
        #ppobject.defaults()
        #ppobject.createregion()
        #ppobject.query(fileObj)

        if None: # parking lot for other tests
            ppobject.timeseries(fileObj)
            ppobject.probe(fileObj)
            ppobject.write(fileObj,'gmv')
            
            ppobject.createregion()
            ppobject.listregions()
            ppobject.deleteregions()
            ppobject.listregions()

            ppobject.stat(fileObj)
            ppobject.listfiles()

            #ppobject.listvariables()
            #ppobject.str()

    except:
        print >> fpwatch, "---> Test failed in some aspect <---"
        print >> fpwatch,"\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise




