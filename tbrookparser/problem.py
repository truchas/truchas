import os,sys,string,re

class Problem:
    #define a Truchas simulation problem to be compiled, run and plotted 
    
    def compile(self):
       os.chdir(self.truchasdir)
       cmd    = 'make -f GNUmakefile ' + self.make
       proc   = os.popen(cmd)
       output = proc.readlines()
       proc.close()

    def run(self,example):

        exedir   = self.truchasdir + '/bin'
        inptdir  = self.problemsdir+'%s' % (example)
        inptfile = '%s.inp' % (example)
        otptopt  = ' -o:output '
        cmd      = exedir + '/truchas ' + otptopt + inptfile
        os.chdir(inptdir)

        print cmd
        proc     = os.popen(cmd)	
        output   = proc.readlines()
        proc.close

        os.chdir(self.problemsbase)

    def plots(self,example,plotno):

        def gmvtolatexfig(self,gmvfile,gmvattr,xwidth,ywidth,pngfile,outputno,plotmod,max,crop="-trim",cb=0,time=0):

            #produce png file for user manual given gmv and gmv attribute files 

            #first modify gmv attribute file to ensure no colorbar and no time is displayed

            S = open(gmvattr,'r').read()
            if not cb:
                S = re.sub('colorbarflag 1','colorbarflag 0',S)
            if not time:
                S = re.sub('timeflag 1','timeflag 0',S)
            S = re.sub('cellnumflag 1','cellnumflag 0',S)
            S = re.sub('axisflag 1','axisflag 0',S)
            O = open(gmvattr,'w')
            O.write(S)
            O.close()

            #convert gmv output to rgb  

            cmd    = 'gmvbatch -w 0 0 %g %g -a %s -i %s -s tmp.rgb' % (xwidth, ywidth, gmvattr, gmvfile)
            cmd    = self.pltpackage.dir + cmd
            proc   = os.popen(cmd)
            output = proc.readlines()
            proc.close()
            print cmd

            #convert rgb file to png 
          
            cmd     = 'convert tmp.rgb %s %s' % (crop,pngfile)
            proc    = os.popen(cmd)
            output  = proc.readlines()
            proc.close()
            print cmd
            
            #copy chosen png files to latex figures directory 

            if outputno % plotmod == 0 or outputno == 1 or outputno == max:
                cmd     = 'cp %s ' + self.figuredir + '%s'
                cmd     = cmd % (pngfile,pngfile)
                proc    = os.popen(cmd)
                output  = proc.readlines()
                proc.close()     

        if self.pltpackage.name == 'gmv':        

          #now create GMV files using the postprocessor

          #change to output directory 
          outputdir = self.problemsdir + '%s/output/' % (example)
          os.chdir(outputdir)
          print os.getcwd()

          try:
              cmd3   = 'python %s/tools/scripts/TBrookParse.py -f gmv.mac' %(self.truchasdir)
              print cmd3
              proc   = os.popen(cmd3)
              output = proc.readlines()
              proc.close
          except:
              print
              print 'Problems creating .graphics.gmv files from Truchas output'
              print 'will not continue to generate User Manual'
              print
              return

          #find maximum output number for this problem 
          lscmd   = 'ls -al %s.ts.*' % (example)
          proc    = os.popen(lscmd)
          output  = proc.readlines()
          max     = len(output)
          max     = 1
          if (max < plotno) :
              plotno = max
          plotmod = max/plotno

          #for each problem define the cycle numbers to be visualised
          D = {'hc':[0,1,2,4,6,8,10],
               'hc_hs_ipc':[0,2,4],
               'hc_bapc':[0,1,5,10],
               'chem':[0,1,2,3,4,5],
               'flow':[0,1,5,8],
               'hc_mipc_flow':[0,1,4],
               'disp':[0,1],
               'vp':[0,3],
               'ih':[0,2,4,6],
               'rad':[0,5,10],
               'ds':[1,5,15,28]}
              
          for outputno in D[example]:
              
            #interrogate extension on each gmv file
            
            ext = '%06g' % (outputno)

            gmvfile  = '%s.graphics.gmv.%s' % (example, ext)
            gmvattr  = '%s-gmv.attr' % (example)
            pngfile  = '%sgmv%s.png' % (example, ext)
            xpixel   = 600
            ypixel   = 600
            latexdir = self.figuredir
            
            if example == 'ih':

                #special case for ih problem - employ additional joule and title files 

                xpixel  = 450
                ypixel  = 520
                textopt = "-font Helvetica -pointsize 32"
                time    = 5*outputno
                crop    = '''%s -draw "text 20,500 't = %i sec'"''' %(textopt,time)
                pngfile = 'ih-temp-%i.png' %(D[example].index(outputno))
                gmvtolatexfig(self, gmvfile, gmvattr, xpixel, ypixel, pngfile, outputno, plotmod, max, crop)

                for att in ('joule', 'title'):

                    pngfile = 'ih-%s.png' %(att)
                    gmvattr = '%s%s-gmv.attr' % (example, att)                    

                    if att == 'joule':
                        xpixel  = 600
                        ypixel  = 600
                        showcb  = 1 #i.e provide colorbar in plot
                        if outputno == 0:gmvtolatexfig(self, gmvfile, gmvattr, xpixel, ypixel, pngfile, outputno, plotmod, max, cb=showcb)
                    else:
                        xpixel  = 900
                        ypixel  = 1040
                        showcb  = 0 #i.e do not provide colorbar in plot
                        if outputno == 2:gmvtolatexfig(self, gmvfile, gmvattr, xpixel, ypixel, pngfile, outputno, plotmod, max, cb=showcb)
                
            elif example == 'vp':

                #special case for vp problem - employ multiple .attr files

                xpixel   = 600
                ypixel   = 600
                gmvtolatexfig(self, gmvfile, gmvattr, xpixel, ypixel, pngfile, outputno, plotmod, max)

                for att in ('r','dz','sxx','T'):      

                   gmvattr  = '%s%s-gmv.attr' % (example, att)
                   pngfile  = '%s%sgmv%s.png' % (example, att, ext)
                   gmvtolatexfig(self, gmvfile, gmvattr, xpixel, ypixel, pngfile, outputno, plotmod, max)


            elif example == 'disp':

                #special case for disp problem - employ multiple .attr files

                xpixel   = 600
                ypixel   = 600
                gmvtolatexfig(self, gmvfile, gmvattr, xpixel, ypixel, pngfile, outputno, plotmod, max)

                for att in ('dy', 'dz','sxx'):      

                   gmvattr  = '%s%s-gmv.attr' % (example,att)
                   pngfile  = '%s%sgmv%s.png' % (example, att, ext)
                   gmvtolatexfig(self, gmvfile, gmvattr, xpixel, ypixel, pngfile, outputno, plotmod, max)


            elif example == 'hc_hs_ipc':

                #special case for hc_hs_ipc problem - employ problem and sensitivity  .attr files      

                xpixel   = 600
                ypixel   = 600
                gmvtolatexfig(self, gmvfile, gmvattr, xpixel, ypixel, pngfile, outputno, plotmod, max)

                for att in ('_tsens1', '_tsens2'):
                    
                   gmvattr  = '%s%s-gmv.attr' % (example,att)
                   pngfile  = '%s%sgmv%s.png' % (example, att, ext)
                   gmvtolatexfig(self, gmvfile, gmvattr, xpixel, ypixel, pngfile, outputno, plotmod, max)
            
            else:

                #general case...

                xpixel   = 600
                ypixel   = 600
                gmvtolatexfig(self, gmvfile, gmvattr, xpixel, ypixel, pngfile, outputno, plotmod, max)
            
          if example == 'rad':

              #create probe files...
              try:
                  cmd      = 'gnuplot gnuplot.dat'
                  proc     = os.popen(cmd)	
                  output   = proc.readlines()
                  proc.close

                  #copy chosen png files to latex figures directory 
                    
                  prbfile = '%sProbe.png' %(example)
                  cmd     = 'cp %s ' + self.figuredir + '%s'
                  cmd     = cmd % (prbfile,prbfile)
                  proc    = os.popen(cmd)
                  output  = proc.readlines()
                  proc.close()
                    
              except:
                  print "no gnuplot file to read probe data"

        if self.pltpackage.name != 'gmv':
          raise 'unable to handle non gmv files' 


    def colorbar(self,example):
        
       if self.pltpackage.name == 'gmv':

          " change to output directory "
          outputdir = self.problemsdir + '%s/output/' % (example)
	  os.chdir(outputdir)

          if example != "ih":
               #change the problem gmv attribute file to ensure only colorbar shows  
               #if example is "ih" then use existing gmv-colorbar attribute file
               file = outputdir + '%s-gmv.attr' % (example)	
               S    = open(file,'r').read()
               S    = re.sub('magnify (.*e)','magnify 0.00000e',S)  
               S    = re.sub('colorbarflag 0','colorbarflag 1',S)
               S    = re.sub('Timeflag 1','timeflag 0',S)
               S    = re.sub('timeflag 1','timeflag 0',S)
               S    = re.sub('cellnumflag 1','cellnumflag 0',S)
               S    = re.sub('axisflag 1','axisflag 0',S) 
               nfil = '%s-gmv-colorbar.attr' % (example)
               O    = open(nfil,'w')
               O.write(S)
               O.close()


	  " produce the rgb file and convert to png "
          cmd    = 'gmvbatch -w 0 0 400 400 -a %s-gmv-colorbar.attr -i %s.graphics.gmv.000001 -s tmpcb.rgb' % (example,example)
          cmd    = self.pltpackage.dir + cmd
          print cmd
          proc   = os.popen(cmd)
          output = proc.readlines()
          proc.close() 

          if example == 'ih':
              crop     = "-crop 70x310+25+77"
              cbarfile = "ih-temp-cbar.png"
          else:
              crop     = "-trim"
              cbarfile = '%scbgmv.png' %(example)
              
          cmd     = 'convert tmpcb.rgb %s %s' % (crop,cbarfile)
          print cmd
          proc    = os.popen(cmd)
          output  = proc.readlines()
          proc.close()

          #cp this file to latex figures dir 
          cmd    = 'cp %s ' + self.figuredir
          cmd    = cmd % (cbarfile)
          proc   = os.popen(cmd)
          output = proc.readlines()
          proc.close()

          #produce additional rgb files and convert to png for the sensitivity section of the hc_hs_ipc problem
          
          if example == 'hc_hs_ipc':

              for att in ('_tsens1', '_tsens2'):

                  " change gmv attribute file to ensure only colorbar shows " 
                  file = outputdir + '%s%s-gmv.attr' % (example,att)	
                  S    = open(file,'r').read()
                  S    = re.sub('magnify (.*e)','magnify 0.00000e',S)  
                  S    = re.sub('colorbarflag 0','colorbarflag 1',S)
                  S    = re.sub('Timeflag 1','timeflag 0',S)
                  S    = re.sub('timeflag 1','timeflag 0',S)
                  S    = re.sub('cellnumflag 1','cellnumflag 0',S)
                  S    = re.sub('axisflag 1','axisflag 0',S) 
                  nfil = '%s%s-gmv-colorbar.attr' % (example,att)
                  O    = open(nfil,'w')
                  O.write(S)
                  O.close()

                  #produce the rgb file and convert to png 
                  cmd    = 'gmvbatch -a %s%s-gmv-colorbar.attr -i %s.graphics.gmv.000000 -s tmpcb.rgb' % (example,att,example)
                  cmd    = self.pltpackage.dir + cmd
                  proc   = os.popen(cmd)
                  output = proc.readlines()
                  proc.close()
                  print cmd

                  cmd     = 'convert -trim tmpcb.rgb %s%scbgmv.png' % (example,att)
                  proc    = os.popen(cmd)
                  output  = proc.readlines()
                  proc.close()
                  print cmd

                  #cp this file to latex dir 
                  cmd    = 'cp %s%scbgmv.png ' + self.figuredir
                  cmd    = cmd % (example,att)
                  proc   = os.popen(cmd)
                  output = proc.readlines()
                  proc.close()

       if self.pltpackage.name != 'gmv':
	  raise 'unable to handle non gmv files'


          
    def input(self,example):

        #create latex formatted Truchas input file 

	#change to output directory
	thisdir = self.problemsdir + '%s/' % (example)
	os.chdir(thisdir)

        #run input file throught InputConvert.pl
	cmd = self.truchasdir + '/tools/scripts/InputConvert.pl -i %s.inp -o %s_conv.inp' % (example,example)
        proc    = os.popen(cmd)
        output  = proc.readlines()
        proc.close()

        #add formats for latexing this section
	file = '%s_conv.inp' % (example)
	nfil = '%s.inp.tex' % (example)
        I    = open(file,'r')
        O    = open(nfil,'w')
        max  = len(I.readlines())
        I.close()
        T    = "\\inputfilelisting { "
        I    = open(file,'r')
        for line in range(max):
           S = I.readline()
           if line > 1 :
  	     S = S + '\\' + '\\'
             S = re.sub('&','\&',S)
             S = re.sub('_','\protect\_',S)
             if re.search("[^']=[^']", S):
                S = "\\>" + re.sub("[^']=[^']","\\>=",S)
             if string.count(S,",") == 2 and string.count(S,"=") == 0:   
                 S = "\\>" + "\\>" + S
             T = T + S    
        T = T + " } "
	O.write(T)
	O.close()
        #cp the formatted input file to the latex inputfiles dir
        cmd     = 'cp ' + nfil + ' ' + self.inputdir 
        proc    = os.popen(cmd)
        output  = proc.readlines()
        proc.close()

    def clean(self,example):

        #clean up files produced by the script

        #change to problem output directory 
        outputdir = self.problemsdir + '%s/output/' % (example)
        os.chdir(outputdir)
        print os.getcwd()

        #rm bin, gmv, xml, rgb files produced by running this script

        for ext in ['.bin', 'gmv*.png', '.xml', '.rgb', '.tex']:
            try:
                cmd    = 'rm *%s' %(ext)
                proc   = os.popen(cmd)
                output = proc.readlines()
                proc.close
            except:
                print '*.%s files do not exist in this directory' %(ext)

class UMProblem(Problem):

      def __init__(self,truchasdir):

	gmv         = PlottingPackage()
        gmv.setdir('/usr/local/bin/')
        gmv.setname('gmv')

	self.pltpackage   = gmv
        self.truchasdir   = truchasdir
	self.make         = 'all-serial-opt'
	self.figuredir     = self.truchasdir + '/documentation/src/UserManual/figures/'
        self.inputdir     = self.truchasdir + '/documentation/src/UserManual/inputfiles/'
	self.problemsbase = self.truchasdir + '/documentation/src/UserManual/UMProblems/' 
	self.problemsdir  = self.problemsbase + 'Chapter_' 

class PlottingPackage:

      def setdir(self,dir):
         self.dir = dir  

      def setname(self,name):
         self.name = name



