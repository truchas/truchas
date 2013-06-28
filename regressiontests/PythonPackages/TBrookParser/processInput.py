"""
COMMENTS HERE PLEASE!!!!

"""
import os, sys, string
from POSTPROCESSORutils     import usubs,createFileObject
from getPostProcessorObject import getPostProcessorObject

class processInput:

    def __init__(self,inbuffer,
                 fp          = sys.stdout,
                 commands    = {},
                 files       = [],
                 def_file    = os.path.abspath(os.path.dirname(__file__)) + '/../../scripts/test_TBrookParse/samples/phasechange_mixed_output/phasechange_mixed.TBrook.xml',
                 def_exofile = os.path.abspath(os.path.dirname(__file__)) + '/../../scripts/test_TBrookParse/samples/hytec_quarter_coarse_side_1.exo'):

        iExit  = None
        iCount = 0
        X      = createFileObject()
        inFile = def_file
        input  = usubs.input(inbuffer,debug=0)

        #check if there is a TBrook.xml file in the current working directory
        for file in os.listdir(os.getcwd()):
            if file[-11:] == '.TBrook.xml':
                inFile  = file

        if (len(files)):
            def_file   = files[0]
            ppobject   = getPostProcessorObject(def_file,def_exofile,input,fp,debug=0)
            X          = ppobject.load(def_file)
            inFile     = X.file
    
        if (not len(files)):
            if not len(inbuffer):
                print >> fp 
                print >> fp, ' Try "load" or "help"'
                print >> fp
            ppobject   = getPostProcessorObject(def_file,def_exofile,input,fp,debug=0)

        while (iExit != 1):
    
            command = ''
            command = input.usetc('Command:',command)

            if (command == 'quit' or command=='q' or command=='e' or command=='end'):

                print >> fp
                print >> fp, 'Finished'
                print >> fp
        
                iExit = 1
        
            elif (command == 'h' or command == 'help'):
        
                ppobject.help(commands)

            elif (command == 'defaults'):

                print >> fp
                print >> fp, 'In postprocessor..setting defaults...'
                print >> fp

                try:
                    ppobject.defaults()
                except:
                    print >> fp
                    print >> fp, 'Error in setting postprocessor defaults...'
                    print >> fp

            elif (command == 'stat' or command=='statistics'):
        
                print >> fp
                print >> fp, 'Chosen to query statistics..'
                print >> fp

                try:
                    inFile = ppobject.input.usetc('Name of file to be read:', inFile)
                    if inFile not in ppobject.files:
                        X  = ppobject.load(inFile)
                    else:
                        for thisfile in ppobject.files:
                            if inFile == thisfile:
                                X.file    = thisfile
                                X.storage = ppobject.files_storages[X.file][0]
                                X.regions = ppobject.files_regions[X.file]
                    inFile = X.file
                    ppobject.stat(X)
                except:
                    print >> fp
                    print >> fp, 'Error in viewing statistics from this file.'
                    print >> fp, 'Try viewing statistics from a different file.'
                    print >> fp


            elif (command == 'probe'):
        
                print >> fp
                print >> fp, 'Chosen to probe..'
                print >> fp

                try:
                    inFile = ppobject.input.usetc('Name of file to be read:', inFile)
                    if inFile not in ppobject.files:
                        X  = ppobject.load(inFile)
                    else:
                        for thisfile in ppobject.files:
                            if inFile == thisfile:
                                X.file    = thisfile
                                X.storage = ppobject.files_storages[X.file][0]
                                X.regions = ppobject.files_regions[X.file]
                    inFile = X.file
                    ppobject.probe(X)
                except:
                    print >> fp
                    print >> fp, 'Error in viewing probes from this file.'
                    print >> fp, 'Try viewing probes from a different file.'
                    print >> fp

            elif (command == 'query'):
        
                print >> fp
                print >> fp, 'Chosen to query field variables..'
                print >> fp

                try:
                    inFile = ppobject.input.usetc('Name of file to be read:', inFile)
                    if inFile not in ppobject.files:
                        X  = ppobject.load(inFile)
                    else:
                        for thisfile in ppobject.files:
                            if inFile == thisfile:
                                X.file    = thisfile
                                X.storage = ppobject.files_storages[X.file][0]
                                X.regions = ppobject.files_regions[X.file]
                    inFile = X.file
                    ppobject.query(X)
                except:
                    print >> fp
                    print >> fp, 'Error in querying variables from this file.'
                    print >> fp, 'Try querying from a different file.'
                    print >> fp
        
            elif (command=='load'):
        
                try:
                    inFile = ppobject.input.usetc('Name of file to be read:', inFile)
                    if inFile not in ppobject.files:
                        X  = ppobject.load(inFile)
                except:
                    print >> fp
                    print >> fp, 'Error in loading file. Try loading a different file'
                    print >> fp

            elif (command=='timeseries'):

                print >> fp
                print >> fp, 'Chosen to create a time series...'
                print >> fp

                try:
                    inFile = ppobject.input.usetc('Name of file to be read:', inFile)
                    if inFile not in ppobject.files:
                        X  = ppobject.load(inFile)
                    else:
                        for thisfile in ppobject.files:
                            if inFile == thisfile:
                                X.file    = thisfile
                                X.storage = ppobject.files_storages[X.file][0]
                                X.regions = ppobject.files_regions[X.file]
                    inFile = X.file
                    ppobject.timeseries(X)
                except:
                    print >> fp
                    print >> fp, 'Error in creating time series.'
                    print >> fp, 'Try viewing a different variable or file.'
                    print >> fp
         
            elif (command=='write'):
        
                print >> fp
                print >> fp, 'Chosen to write graphics...'
                print >> fp

                writeWhat      = input.usetc('write what? gmv/ensight/tecplot/vtk:','gmv')
                choice         = writeWhat.lower()
        
                if choice == 'opendx':
                    print >> fp, '\n Unable to output opendx files. Try ensight or gmv. Will choose gmv by default. \n'
                    choice = 'gmv'

                try:
                    inFile = ppobject.input.usetc('Name of file to be read:', inFile)
                    if inFile not in ppobject.files:
                        X  = ppobject.load(inFile)
                    else:
                        for thisfile in ppobject.files:
                            if inFile == thisfile:
                                X.file    = thisfile
                                X.storage = ppobject.files_storages[X.file][0]
                                X.regions = ppobject.files_regions[X.file]
                    inFile = X.file
                    ppobject.write(X,choice)
                except:
                    print >> fp
                    print >> fp, 'Error in writing graphics file.'
                    print >> fp, 'Try viewing a different file.'
                    print >> fp
            
            elif (command=='restart'):

                print >> fp
                print >> fp, 'Chosen to restart...'
                print >> fp 

                try:
                    rinFiles = []
                    rinFiles.append(inFile)
                    prompt   = 'Name(s) of file(s) to be read (file1,file2): \n'
                    rinFiles = ppobject.input.usetc(prompt, rinFiles[0])
                    if rinFiles[len(rinFiles)-1] == ',':
                        print >> fp
                        print >> fp, 'When providing multiple files please ensure there is only '
                        print >> fp, 'a comma and no spaces seperating the files.'
                        print >> fp, 'Input should be in the form: file1,file2'
                        print >> fp 
                        sys.exit(1)
                    rinFiles = string.replace(rinFiles,' ','')
                    rinFiles = string.split(rinFiles,',')
                    Xs       = []
                    cnt      = 0
                    for inFile in rinFiles:
                        if inFile not in ppobject.files:
                            X   = ppobject.load(inFile)
                            Xs.append(X)
                            cnt += 1
                        else:
                            X   = createFileObject()
                            Xs.append(X)
                            for thisfile in ppobject.files:
                                if inFile == thisfile:
                                    Xs[cnt].file    = thisfile
                                    Xs[cnt].storage = ppobject.files_storages[Xs[cnt].file][0]
                                    Xs[cnt].regions = ppobject.files_regions[Xs[cnt].file]
                                    cnt            += 1
                    inFile = Xs[0].file
                    ppobject.restart(Xs)
                except:
                    print >> fp
                    print >> fp, 'Error in writing restart file(s).'
                    print >> fp, 'Try restarting from a different file'
                    print >> fp 

            elif (command=='list'):
        
                what = input.usetc('List what (files/variables/regions):','regions')
                if   (what == 'regions'):
                    ppobject.listregions()
                elif (what == 'variables'):
                    ppobject.listvariables()
                elif (what == 'files'):
                    ppobject.listfiles()
            
            elif (command=='define' or command=='region'):

                print >> fp
                print >> fp, 'Chosen to define a region...'
                print >> fp 
                ppobject.createregion()
        
            elif (command=='deleteAllRegions'):
        
                ppobject.deleteregions()

            elif (not len(command)):

                iCount    = iCount + 1
                if(iCount%1000 == 0):
                    print >> fp, '\n\n  Too many invalid commands.  Exiting.\n'
                    iExit = 1

            elif (len(command)):

                print >> fp, '\n  Invalid command: %s, see help for valid commands\n'%command
                iCount = 0
