#!/usr/bin/env python
"""
CheckPPVTK

-----------------------------------------------------------------------------
   Purpose:
  
      TestCase to check validity of VTK file (vtk.bin, vtk.ascii) created from
      vtk.mac file 
  
   Public Interface:
  
      T = CheckPPVTK(unittest.TestCase,TruchasBaseTest)
      T.setUp()
      T.shortDescription()
      T.testVTKFile()
      T.tearDown()
      T.str()
  
   Contains:  
      class CheckPPVTK
            __init__(dir,RunSpecs)
            setUp()
            shortDescription()
            testVTKFile()
            tearDown()
            str()
  
   Unit Test Block
  
   Author(s): Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""

import unittest
import os, sys, string, re, fnmatch, shutil, platform

if __name__ == '__main__':
    print "\n for component test in %s \n" %(__file__)
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    testingdir  = os.path.abspath(thisdir + '/../')
    sys.path.append(testingdir)
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

from Runners            import TestRunner, RunTimeSpecs
from TruchasBaseTest    import TruchasBaseTest

class CheckPPVTK(unittest.TestCase,TruchasBaseTest):

    def __init__(self,dir,RunSpecs):

	unittest.TestCase.__init__(self,methodName='testVTKFile')
	self.dir          = dir
	#specifications from a previous Truchas run
	self.runspecs     = RunSpecs        
	#for formatting log file
	self.column1      = 0
	self.column2      = 5	
	self.logdir       = self.runspecs.logdir
        self.logname      = 'CheckPPVTK.log'
        self.basenamelog  = string.split(self.runspecs.logfile,'/')[1]
	self.logfile      = os.path.join(self.logdir,self.basenamelog,self.logname)
	self.outfile      = 'checkppvtk.out'
        self.vtkbinfile   = '*bin*.vtk'
        self.vtkasciifile = '*ascii.vtk'
        self.vtkfiles     = []
	self.batchfile    = 'PVvtktest'
        self.curfile      = None

    def setUp(self):
	"defines vtk file location" 
	
	os.chdir(self.dir)

        L      = string.split(self.logdir,'/')
        logdir = L[-1]
        
        if logdir not in os.listdir(self.dir):
            os.mkdir(logdir)

        wdir = os.path.join(self.dir,self.logdir)
	os.chdir(wdir)

        if self.basenamelog not in os.listdir(wdir):
            os.mkdir(self.basenamelog)
        
        logfl          = os.path.join(self.basenamelog,self.logname)

	self.watcher   = open(logfl,'w')
        
	self.watcher.write('In CheckPPVTK TestCase')
 	self.watcher.write('\n')
 	self.watcher.write('\n')

	tmp     = '*Setting up'
	str     = self.rjustln(self.column1, tmp)
	self.watcher.write(str)
 	self.watcher.write('\n')

    def shortDescription(self):
	"modified unittest short description to allow class name and input file name"

        classname = self.__class__.__name__
        classname = re.sub('Test[ers]*$', '', classname)
        docname   = '['+classname +']' + ' ' 
        return  docname.ljust(40) + ': ' + str(unittest.TestCase.shortDescription(self))

    def testVTKFile(self):
	"tests validity of VTK file created from a postprocessor run or from self.plotter"

	# write to watcher
	tmp     = '*testVTKFile'
	str     = self.rjustln(self.column1, tmp)
	self.watcher.write(str)
 	self.watcher.write('\n')

	# setup error status and message
        success = 0
        tmp     = 'ERROR!! No VTK file created by Truchas postprocessor run or from self.plotter'
        errstr  = '\n \n' + tmp + '\n'
	err     = ''

        for file in os.listdir(os.getcwd()):
            if fnmatch.fnmatch(file,self.vtkbinfile) or fnmatch.fnmatch(file,self.vtkasciifile):
                
		self.batchfile = file
                self.curfile = file
                self.vtkfiles.append(file)
	        self.watcher.write(self.str(self.column2))

		# get the paraview/bin path
		# write the pvbatch batch script
		self.writePVBatchFile(file)
		cmd = '%s/bin/pvbatch %s.pvb 1>stdout 2>stderr' % (os.getenv('PARAVIEW_HOME'),self.batchfile)
		os.system(cmd)

		# test for success
		jpgout = self.batchfile + '.jpg'
		for jpgfile in os.listdir(os.getcwd()):
		    if jpgout in jpgfile:
                        success = 1
		        if jpgfile not in self.vtkfiles:
                            self.vtkfiles.append(jpgfile) 
		        tmp  = 'Output File     : '
			tmp += jpgfile
		        str  = self.rjustln(self.column2, tmp)
                        self.watcher.write(str)
                        
		# write result to watcher
		if success: 
		    tmp = 'Status          : PASSED'
		    str = self.rjustln(self.column2, tmp)
                    self.watcher.write(str)
 	            self.watcher.write('\n')
		else:       
		    tmp    = 'Status          : FAILED -> STDERR '
		    str    = self.rjustln(self.column2, tmp)
                    self.watcher.write(str)
		    # if failed, append stderr contents
		    stderr = open('stderr','r')
		    tmp    = stderr.readlines()
                    self.watcher.writelines(tmp)
		    stderr.close()
		    err   += 'ERROR!!  Failed to generate RGB file from %s \n' % file

        if not success:
            if len(err) > 0:
                errstr  = '\n \n' + err + '\n'
	    self.watcher.write(errstr)
            self.watcher.write('\n')
            raise ValueError(errstr)


    def tearDown(self):
	"ensures working directory returned to original directory before test started" 

        tmp     = '*Tearing down'
        str     = self.rjustln(self.column1, tmp)
	self.watcher.write(str)

	self.watcher.write('\n') 
        self.watcher.close()

        os.chdir(self.dir)

    def str(self,column=15):
        "provides formatted info about vtk file created"

        info  = '\n'
	tmp   = 'VTK file(s)     : ' + str(self.vtkfiles)
        info += self.rjustln(column, tmp)
	tmp   = 'Log dir         : ' + self.runspecs.logdir
        info += self.rjustln(column, tmp)

        return info

    def writePVBatchFile(self,input):
        "writes out a PVbatch script with the current input filename"

	batchfile = self.batchfile + '.pvb'
	fp = open(batchfile,'w')
        print >> fp,'# ParaView Version 2.4'
        print >> fp,''
        print >> fp,''
        print >> fp,'#Initialization'
        print >> fp,''
        print >> fp,'vtkSMObject foo'
        print >> fp,'set proxyManager [foo GetProxyManager]'
        print >> fp,''
        print >> fp,'set smApplication [foo GetApplication]'
        print >> fp,'$smApplication ParseConfigurationFiles'
        print >> fp,'foo Delete'
        print >> fp,''
        print >> fp,'vtkSMProperty foo'
        print >> fp,'foo SetCheckDomains 0'
        print >> fp,'foo Delete'
        print >> fp,''
        print >> fp,''
        print >> fp,'set pvTemp506 [$proxyManager NewProxy sources legacyreader]'
        print >> fp,'  $proxyManager RegisterProxy sources pvTemp506 $pvTemp506'
        print >> fp,'  $pvTemp506 UnRegister {}'
        print >> fp,'  [$pvTemp506 GetProperty FileName] SetElement 0 {%s}' % input
        print >> fp,'  $pvTemp506 UpdateVTKObjects'
        print >> fp,''
        print >> fp,'set pvTemp604 [$proxyManager NewProxy lookup_tables LookupTable]'
        print >> fp,'  $proxyManager RegisterProxy lookup_tables pvTemp604 $pvTemp604'
        print >> fp,'  $pvTemp604 UnRegister {}'
        print >> fp,'  [$pvTemp604 GetProperty ArrayName] SetElement 0 {Enthalpy}'
        print >> fp,'  [$pvTemp604 GetProperty NumberOfTableValues] SetElements1 256'
        print >> fp,'  [$pvTemp604 GetProperty HueRange] SetElements2 0.6667 0'
        print >> fp,'  [$pvTemp604 GetProperty SaturationRange] SetElements2 1 1'
        print >> fp,'  [$pvTemp604 GetProperty ValueRange] SetElements2 1 1'
        print >> fp,'  [$pvTemp604 GetProperty ScalarRange] SetElements2 5580 5580'
        print >> fp,'  [$pvTemp604 GetProperty VectorComponent] SetElements1 0'
        print >> fp,'  [$pvTemp604 GetProperty VectorMode] SetElements1 1'
        print >> fp,'  $pvTemp604 UpdateVTKObjects'
        print >> fp,''
        print >> fp,''
        print >> fp,'set pvTemp611 [$proxyManager NewProxy properties TextProperty]'
        print >> fp,'  $proxyManager RegisterProxy properties pvTemp611 $pvTemp611'
        print >> fp,'  $pvTemp611 UnRegister {}'
        print >> fp,'  [$pvTemp611 GetProperty Bold] SetElement 0 1'
        print >> fp,'  [$pvTemp611 GetProperty Color] SetElement 0 1'
        print >> fp,'  [$pvTemp611 GetProperty Color] SetElement 1 1'
        print >> fp,'  [$pvTemp611 GetProperty Color] SetElement 2 1'
        print >> fp,'  [$pvTemp611 GetProperty FontFamily] SetElement 0 0'
        print >> fp,'  [$pvTemp611 GetProperty FontSize] SetElement 0 12'
        print >> fp,'  [$pvTemp611 GetProperty Italic] SetElement 0 1'
        print >> fp,'  [$pvTemp611 GetProperty Opacity] SetElement 0 1'
        print >> fp,'  [$pvTemp611 GetProperty Shadow] SetElement 0 1'
        print >> fp,'  $pvTemp611 UpdateVTKObjects'
        print >> fp,'set pvTemp612 [$proxyManager NewProxy properties TextProperty]'
        print >> fp,'  $proxyManager RegisterProxy properties pvTemp612 $pvTemp612'
        print >> fp,'  $pvTemp612 UnRegister {}'
        print >> fp,'  [$pvTemp612 GetProperty Bold] SetElement 0 1'
        print >> fp,'  [$pvTemp612 GetProperty Color] SetElement 0 1'
        print >> fp,'  [$pvTemp612 GetProperty Color] SetElement 1 1'
        print >> fp,'  [$pvTemp612 GetProperty Color] SetElement 2 1'
        print >> fp,'  [$pvTemp612 GetProperty FontFamily] SetElement 0 0'
        print >> fp,'  [$pvTemp612 GetProperty FontSize] SetElement 0 12'
        print >> fp,'  [$pvTemp612 GetProperty Italic] SetElement 0 1'
        print >> fp,'  [$pvTemp612 GetProperty Opacity] SetElement 0 1'
        print >> fp,'  [$pvTemp612 GetProperty Shadow] SetElement 0 1'
        print >> fp,'  $pvTemp612 UpdateVTKObjects'
        print >> fp,''
        print >> fp,'set pvTemp609 [$proxyManager NewProxy displays ScalarBarWidget]'
        print >> fp,'  $proxyManager RegisterProxy displays pvTemp609 $pvTemp609'
        print >> fp,'  $pvTemp609 UnRegister {}'
        print >> fp,'# Input to Display Proxy not set properly or takes no Input.'
        print >> fp,'  [$pvTemp609 GetProperty Visibility] SetElement 0 1'
        print >> fp,'  [$pvTemp609 GetProperty LabelFormat] SetElement 0 {%-#6.3g}'
        print >> fp,'  [$pvTemp609 GetProperty LabelTextProperty] AddProxy $pvTemp611'
        print >> fp,'  [$pvTemp609 GetProperty LookupTable] AddProxy $pvTemp604'
        print >> fp,'  [$pvTemp609 GetProperty Orientation] SetElement 0 1'
        print >> fp,'  [$pvTemp609 GetProperty Position] SetElement 0 0.87'
        print >> fp,'  [$pvTemp609 GetProperty Position] SetElement 1 0.25'
        print >> fp,'  [$pvTemp609 GetProperty Position2] SetElement 0 0.13'
        print >> fp,'  [$pvTemp609 GetProperty Position2] SetElement 1 0.5'
        print >> fp,'  [$pvTemp609 GetProperty Title] SetElement 0 {Enthalpy}'
        print >> fp,'  [$pvTemp609 GetProperty TitleTextProperty] AddProxy $pvTemp612'
        print >> fp,'  $pvTemp609 UpdateVTKObjects'
        print >> fp,'#Display Proxy'
        print >> fp,''
        print >> fp,'set pvTemp512 [$proxyManager NewProxy displays LODDisplay]'
        print >> fp,'  $proxyManager RegisterProxy displays pvTemp512 $pvTemp512'
        print >> fp,'  $pvTemp512 UnRegister {}'
        print >> fp,'  [$pvTemp512 GetProperty Input]  AddProxy $pvTemp506'
        print >> fp,'  # skipping not-saveable property CacheUpdate'
        print >> fp,'  # skipping not-saveable property InvalidateGeometry'
        print >> fp,'  [$pvTemp512 GetProperty LODResolution] SetElement 0 50'
        print >> fp,'  [$pvTemp512 GetProperty Representation] SetElement 0 2'
        print >> fp,'  # skipping not-saveable property ResetTransferFunctions'
        print >> fp,'  # skipping not-saveable property Update'
        print >> fp,'  [$pvTemp512 GetProperty Visibility] SetElement 0 1'
        print >> fp,'  [$pvTemp512 GetProperty ColorSpace] SetElement 0 0'
        print >> fp,'  [$pvTemp512 GetProperty HSVWrap] SetElement 0 0'
        print >> fp,'  [$pvTemp512 GetProperty UseStrips] SetElement 0 0'
        print >> fp,'  [$pvTemp512 GetProperty ColorArray] SetElement 0 {Enthalpy}'
        print >> fp,'  [$pvTemp512 GetProperty ColorMode] SetElement 0 0'
        print >> fp,'  [$pvTemp512 GetProperty ImmediateModeRendering] SetElement 0 1'
        print >> fp,'  [$pvTemp512 GetProperty InterpolateScalarsBeforeMapping] SetElement 0 1'
        print >> fp,'  [$pvTemp512 GetProperty LookupTable] AddProxy $pvTemp604'
        print >> fp,'  [$pvTemp512 GetProperty ScalarMode] SetElement 0 4'
        print >> fp,'  [$pvTemp512 GetProperty ScalarVisibility] SetElement 0 1'
        print >> fp,'  [$pvTemp512 GetProperty UseLookupTableScalarRange] SetElement 0 1'
        print >> fp,'  [$pvTemp512 GetProperty Orientation] SetElement 0 0'
        print >> fp,'  [$pvTemp512 GetProperty Orientation] SetElement 1 0'
        print >> fp,'  [$pvTemp512 GetProperty Orientation] SetElement 2 0'
        print >> fp,'  [$pvTemp512 GetProperty Origin] SetElement 0 0'
        print >> fp,'  [$pvTemp512 GetProperty Origin] SetElement 1 0'
        print >> fp,'  [$pvTemp512 GetProperty Origin] SetElement 2 0'
        print >> fp,'  [$pvTemp512 GetProperty Position] SetElement 0 0'
        print >> fp,'  [$pvTemp512 GetProperty Position] SetElement 1 0'
        print >> fp,'  [$pvTemp512 GetProperty Position] SetElement 2 0'
        print >> fp,'  [$pvTemp512 GetProperty Scale] SetElement 0 1'
        print >> fp,'  [$pvTemp512 GetProperty Scale] SetElement 1 1'
        print >> fp,'  [$pvTemp512 GetProperty Scale] SetElement 2 1'
        print >> fp,'  [$pvTemp512 GetProperty Color] SetElement 0 1'
        print >> fp,'  [$pvTemp512 GetProperty Color] SetElement 1 1'
        print >> fp,'  [$pvTemp512 GetProperty Color] SetElement 2 1'
        print >> fp,'  [$pvTemp512 GetProperty Interpolation] SetElement 0 1'
        print >> fp,'  [$pvTemp512 GetProperty LineWidth] SetElement 0 1'
        print >> fp,'  [$pvTemp512 GetProperty Opacity] SetElement 0 1'
        print >> fp,'  [$pvTemp512 GetProperty PointSize] SetElement 0 1'
        print >> fp,'  [$pvTemp512 GetProperty Specular] SetElement 0 0'
        print >> fp,'  [$pvTemp512 GetProperty ScalarMode] SetElement 0 4'
        print >> fp,'  [$pvTemp512 GetProperty SelectScalarArray] SetElement 0 {}'
        print >> fp,'  [$pvTemp512 GetProperty ScalarMode] SetElement 0 4'
        print >> fp,'  [$pvTemp512 GetProperty SelectScalarArray] SetElement 0 {}'
        print >> fp,'  [$pvTemp512 GetProperty ScalarOpacityUnitDistance] SetElement 0 1'
        print >> fp,'  [$pvTemp512 GetProperty ScalarMode] SetElement 0 4'
        print >> fp,'  [$pvTemp512 GetProperty SelectScalarArray] SetElement 0 {}'
        print >> fp,'  $pvTemp512 UpdateVTKObjects'
        print >> fp,''
        print >> fp,'set pvTemp110 [$proxyManager NewProxy axes Axes]'
        print >> fp,'  $proxyManager RegisterProxy axes pvTemp110 $pvTemp110'
        print >> fp,'  $pvTemp110 UnRegister {}'
        print >> fp,'# Input to Display Proxy not set properly or takes no Input.'
        print >> fp,'  [$pvTemp110 GetProperty Origin] SetElement 0 0'
        print >> fp,'  [$pvTemp110 GetProperty Origin] SetElement 1 0'
        print >> fp,'  [$pvTemp110 GetProperty Origin] SetElement 2 0'
        print >> fp,'  [$pvTemp110 GetProperty Position] SetElement 0 250'
        print >> fp,'  [$pvTemp110 GetProperty Position] SetElement 1 5'
        print >> fp,'  [$pvTemp110 GetProperty Position] SetElement 2 5'
        print >> fp,'  [$pvTemp110 GetProperty Scale] SetElement 0 125'
        print >> fp,'  [$pvTemp110 GetProperty Scale] SetElement 1 2.5'
        print >> fp,'  [$pvTemp110 GetProperty Scale] SetElement 2 2.5'
        print >> fp,'  [$pvTemp110 GetProperty Visibility] SetElement 0 1'
        print >> fp,'  $pvTemp110 UpdateVTKObjects'
        print >> fp,'# RenderModule Proxy ---------- '
        print >> fp,'set Ren1 [$proxyManager NewProxy rendermodules LODRenderModule]'
        print >> fp,'  $proxyManager RegisterProxy rendermodules Ren1 $Ren1'
        print >> fp,'  $Ren1 UnRegister {}'
        print >> fp,'  # skipping proxy property CacheUpdate'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp110 } ;#--- Axes part 0'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp122 } ;#--- LODDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp162 } ;#--- CubeAxesDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp165 } ;#--- PointLabelDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp186 } ;#--- LODDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp226 } ;#--- CubeAxesDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp229 } ;#--- PointLabelDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp250 } ;#--- LODDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp290 } ;#--- CubeAxesDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp293 } ;#--- PointLabelDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp313 } ;#--- LineWidgetProxy part 0'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp316 } ;#--- LODDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp356 } ;#--- CubeAxesDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp359 } ;#--- PointLabelDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp380 } ;#--- LODDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp420 } ;#--- CubeAxesDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp423 } ;#--- PointLabelDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp444 } ;#--- LODDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp484 } ;#--- CubeAxesDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp487 } ;#--- PointLabelDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp512 } ;#--- LODDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp570 } ;#--- CubeAxesDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp573 } ;#--- PointLabelDisplay'
        print >> fp,'  catch { [$Ren1 GetProperty {Displays}] AddProxy $pvTemp609 } ;#--- ScalarBarWidget part 0'
        print >> fp,'  # skipping proxy property InteractiveRender'
        print >> fp,'  # skipping proxy property InvalidateGeometry'
        print >> fp,'  [$Ren1 GetProperty {LODResolution}] SetElement 0 50'
        print >> fp,'  [$Ren1 GetProperty {LODThreshold}] SetElement 0 5'
        print >> fp,'  [$Ren1 GetProperty {RenderInterruptsEnabled}] SetElement 0 1'
        print >> fp,'  # skipping proxy property StillRender'
        print >> fp,'  [$Ren1 GetProperty {UseImmediateMode}] SetElement 0 0'
        print >> fp,'  [$Ren1 GetProperty {UseLight}] SetElement 0 0'
        print >> fp,'  [$Ren1 GetProperty {UseTriangleStrips}] SetElement 0 0'
        print >> fp,'  [$Ren1 GetProperty {CameraClippingRange}] SetElement 0 946.699'
        print >> fp,'  [$Ren1 GetProperty {CameraClippingRange}] SetElement 1 993.332'
        print >> fp,'  # skipping proxy property CameraClippingRangeInfo'
        print >> fp,'  [$Ren1 GetProperty {CameraFocalPoint}] SetElement 0 250'
        print >> fp,'  [$Ren1 GetProperty {CameraFocalPoint}] SetElement 1 5'
        print >> fp,'  [$Ren1 GetProperty {CameraFocalPoint}] SetElement 2 5'
        print >> fp,'  # skipping proxy property CameraFocalPointInfo'
        print >> fp,'  [$Ren1 GetProperty {CameraParallelProjection}] SetElement 0 0'
        print >> fp,'  [$Ren1 GetProperty {CameraParallelScale}] SetElement 0 1'
        print >> fp,'  # skipping proxy property CameraParallelScaleInfo'
        print >> fp,'  [$Ren1 GetProperty {CameraPosition}] SetElement 0 250'
        print >> fp,'  [$Ren1 GetProperty {CameraPosition}] SetElement 1 5'
        print >> fp,'  [$Ren1 GetProperty {CameraPosition}] SetElement 2 971.312'
        print >> fp,'  # skipping proxy property CameraPositionInfo'
        print >> fp,'  [$Ren1 GetProperty {CameraViewAngle}] SetElement 0 30'
        print >> fp,'  # skipping proxy property CameraViewAngleInfo'
        print >> fp,'  [$Ren1 GetProperty {CameraViewUp}] SetElement 0 0'
        print >> fp,'  [$Ren1 GetProperty {CameraViewUp}] SetElement 1 1'
        print >> fp,'  [$Ren1 GetProperty {CameraViewUp}] SetElement 2 0'
        print >> fp,'  # skipping proxy property CameraViewUpInfo'
        print >> fp,'  [$Ren1 GetProperty {LightAmbientColor}] SetElement 0 1'
        print >> fp,'  [$Ren1 GetProperty {LightAmbientColor}] SetElement 1 1'
        print >> fp,'  [$Ren1 GetProperty {LightAmbientColor}] SetElement 2 1'
        print >> fp,'  [$Ren1 GetProperty {LightDiffuseColor}] SetElement 0 1'
        print >> fp,'  [$Ren1 GetProperty {LightDiffuseColor}] SetElement 1 1'
        print >> fp,'  [$Ren1 GetProperty {LightDiffuseColor}] SetElement 2 1'
        print >> fp,'  [$Ren1 GetProperty {LightIntensity}] SetElement 0 1'
        print >> fp,'  [$Ren1 GetProperty {LightSpecularColor}] SetElement 0 1'
        print >> fp,'  [$Ren1 GetProperty {LightSpecularColor}] SetElement 1 1'
        print >> fp,'  [$Ren1 GetProperty {LightSpecularColor}] SetElement 2 1'
        print >> fp,'  [$Ren1 GetProperty {LightSwitch}] SetElement 0 1'
        print >> fp,'  [$Ren1 GetProperty {BackLightAzimuth}] SetElement 0 110'
        print >> fp,'  [$Ren1 GetProperty {BackLightElevation}] SetElement 0 0'
        print >> fp,'  [$Ren1 GetProperty {BackLightK:B Ratio}] SetElement 0 3.5'
        print >> fp,'  [$Ren1 GetProperty {BackLightWarmth}] SetElement 0 0.5'
        print >> fp,'  [$Ren1 GetProperty {FillLightAzimuth}] SetElement 0 -10'
        print >> fp,'  [$Ren1 GetProperty {FillLightElevation}] SetElement 0 -75'
        print >> fp,'  [$Ren1 GetProperty {FillLightK:F Ratio}] SetElement 0 3'
        print >> fp,'  [$Ren1 GetProperty {FillLightWarmth}] SetElement 0 0.4'
        print >> fp,'  [$Ren1 GetProperty {HeadLightK:H Ratio}] SetElement 0 3'
        print >> fp,'  [$Ren1 GetProperty {HeadLightWarmth}] SetElement 0 0.5'
        print >> fp,'  [$Ren1 GetProperty {KeyLightAzimuth}] SetElement 0 10'
        print >> fp,'  [$Ren1 GetProperty {KeyLightElevation}] SetElement 0 50'
        print >> fp,'  [$Ren1 GetProperty {KeyLightIntensity}] SetElement 0 0.75'
        print >> fp,'  [$Ren1 GetProperty {KeyLightWarmth}] SetElement 0 0.6'
        print >> fp,'  [$Ren1 GetProperty {MaintainLuminance}] SetElement 0 0'
        print >> fp,'  [$Ren1 GetProperty {FullScreen}] SetElement 0 0'
        print >> fp,'  [$Ren1 GetProperty {OffScreenRendering}] SetElement 0 0'
        print >> fp,'  [$Ren1 GetProperty {RenderWindowSize}] SetElement 0 400'
        print >> fp,'  [$Ren1 GetProperty {RenderWindowSize}] SetElement 1 400'
        print >> fp,'  # skipping proxy property RenderWindowSizeInfo'
        print >> fp,'  [$Ren1 GetProperty {Background}] SetElement 0 0.329412'
        print >> fp,'  [$Ren1 GetProperty {Background}] SetElement 1 0.34902'
        print >> fp,'  [$Ren1 GetProperty {Background}] SetElement 2 0.427451'
        print >> fp,'  # skipping proxy property Viewport'
        print >> fp,'# End of RenderModuleProxy ---- '
        print >> fp,'  [$Ren1 GetProperty OffScreenRendering] SetElement 0 1'
        print >> fp,''
        print >> fp,'set pvTemp114 [$proxyManager NewProxy animation AnimationScene]'
        print >> fp,'$proxyManager RegisterProxy animation pvTemp114 $pvTemp114'
        print >> fp,'[$pvTemp114 GetProperty TimeMode] SetElements1 1'
        print >> fp,'[$pvTemp114 GetProperty StartTime] SetElements1 0'
        print >> fp,'[$pvTemp114 GetProperty EndTime] SetElements1 9'
        print >> fp,'[$pvTemp114 GetProperty AnimatedElement] SetElements1 0'
        print >> fp,'$pvTemp114 UpdateVTKObjects'
        print >> fp,''
        print >> fp,'$pvTemp114 UnRegister {}'
        print >> fp,''
        print >> fp,'  [$pvTemp114 GetProperty Loop] SetElements1 0'
        print >> fp,'  [$pvTemp114 GetProperty FrameRate] SetElements1 1'
        print >> fp,'  [$pvTemp114 GetProperty PlayMode] SetElements1 0'
        print >> fp,'  $pvTemp114 SetRenderModuleProxy $Ren1'
        print >> fp,'  $pvTemp114 UpdateVTKObjects'
        print >> fp,''
        print >> fp,''
        print >> fp,'set saveState 0'
        print >> fp,'for {set i  0} {$i < [expr $argc - 1]} {incr i} {'
        print >> fp,'  if {[lindex $argv $i] == "-XML"} {'
        print >> fp,'    set saveState 1'
        print >> fp,'    set stateName [lindex $argv [expr $i + 1]]'
        print >> fp,'  }'
        print >> fp,'}'
        print >> fp,'if { $saveState } {'
        print >> fp,'   $Ren1 UpdateVTKObjects'
        print >> fp,'   $proxyManager SaveState $stateName'
        print >> fp,'} else {'
        print >> fp,''
        print >> fp,'$Ren1 UpdateVTKObjects'
        print >> fp,'set inBatch 0'
        print >> fp,'for {set i  1} {$i < [expr $argc]} {incr i} {'
        print >> fp,'  if {[lindex $argv $i] == "-BT"} {'
        print >> fp,'    set inBatch 1'
        print >> fp,'  }'
        print >> fp,'}'
        print >> fp,'if { $inBatch } {'
        print >> fp,'  [$Ren1 GetProperty RenderWindowSize]  SetElement 0 300'
        print >> fp,'  [$Ren1 GetProperty RenderWindowSize]  SetElement 1 300'
        print >> fp,'  $Ren1 UpdateVTKObjects'
        print >> fp,'}'
        print >> fp,'$Ren1 StillRender'
        print >> fp,'$Ren1 WriteImage {%s.jpg} vtkJPEGWriter' % self.batchfile
        print >> fp,'}'
        print >> fp,''
        print >> fp,'$proxyManager UnRegisterProxies'
        print >> fp,''
	fp.close()

if __name__=='__main__':
    currdir     = os.getcwd()
    TestRunner  = unittest.TextTestRunner    

    #for testing purposes specify dummy runtime specs 
    for dir in os.listdir(currdir):
	if os.path.isdir(dir) and 'outputs' in dir:
            logdir  = dir
            logfile = logdir+'/static_drop_logs/BasicRun.log'

	    # parse directory name for run time specs
	    vals = string.splitfields(dir,'.')
	    if vals[1][2:] == string.lower(platform.system()):
	        parallel_env = vals[4]
	        np           = 1
	        compiler     = vals[3]
	        compile_mode = vals[5][:-2]
                LastRunSpecs = RunTimeSpecs(logdir,logfile,parallel_env,np,compiler,compile_mode)

                suite        = unittest.TestSuite()
                T            = CheckPPVTK(currdir,LastRunSpecs) 
                suite.addTest(T)
                runner       = TestRunner(verbosity=2)
                result       = runner.run( suite )    
                print T.str()




 
