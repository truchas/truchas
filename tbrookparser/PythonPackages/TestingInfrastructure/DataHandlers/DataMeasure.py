#!/usr/bin/env python
"""
 DataMeasure

-----------------------------------------------------------------------------
   Purpose:
  
      Provides a variety of metrics given a Numeric data array
  
   Public Interface:
  
      T = DataMeasure
      T.meanValue(field)
      T.minValue(field)
      T.maxValue(field)
      T.totalValue(field)
      T.absError(fielda,fieldb)
      T.percentError(first,second)
      T.l1Error(fielda,fieldb)
      T.l2Error(fielda,fieldb)
      T.linfError(fielda,fieldb)
  
   Unit Test Block
  
   Author: Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""

import os, sys
try:
   import numpy.oldnumeric as Numeric
except ImportError:
   import Numeric
except:
   raise

if __name__ == '__main__':
    print "\n for component test in %s \n" %(__file__)
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

from PYTHONutils import getAbortString

class DataMeasure:

    def __init__(self):

        self.abortstring = getAbortString()

    def meanValue(self,field):
	"mean value of a field which must be a Numeric array"

        errstring  = "\n\n                 In meanValue field is not a Numeric array \n"
        errstring += self.abortstring
        
	assert type(field) == type(Numeric.array(1)), errstring

	size = len(field)  

        errstring  = "\n\n                 In meanValue field size <= 0 \n"
        errstring += self.abortstring

	assert size > 0, errstring

	MV   = 1.0*Numeric.sum(field)/size

	return MV

    def minValue(self,field):
	"min value of a field which must be a Numeric array of rank 1"

        errstring  = "\n\n                 In minValue field is not a Numeric array \n"
        errstring += self.abortstring

	assert type(field) == type(Numeric.array(1)), errstring

        errstring  = "\n\n                 In minValue field rank > 1, require field rank = 1 \n"
        errstring += self.abortstring

	assert len(field.shape) == 1, errstring

        minV   = min(field)

	return minV

    def maxValue(self,field):
	"max value of a field which must be a Numeric array or rank 1"

        errstring  = "\n\n                 In maxValue field is not a Numeric array \n"
        errstring += self.abortstring

	assert type(field) == type(Numeric.array(1)), errstring

        errstring  = "\n\n                 In maxValue field rank > 1, require field rank = 1 \n"
        errstring += self.abortstring

	assert len(field.shape) == 1, errstring

        maxV   = max(field)

	return maxV

    def totalValue(self,field):
	"summation of elements in a field which must be a Numeric array"

        errstring  = "\n\n                 In totalValue first field is not a Numeric array \n" 
        errstring += self.abortstring

	assert type(field) == type(Numeric.array(1)), errstring
	
	size = len(field)

        errstring  = "\n\n                 In totalValue field size <= 0 \n"
        errstring += self.abortstring

	assert size > 0, errstring

	T = Numeric.sum(field)
	
	return T

    def absError(self,fielda,fieldb=None):
	"abs error |fielda-fieldb|"

        errstring  = "\n\n                  In absError first field is not a Numeric array \n"
        errstring += self.abortstring

	assert type(fielda) == type(Numeric.array(1)), errstring

	size = len(fielda)

        errstring  = "\n\n                  In absError first field size <= 0 \n" 
        errstring += self.abortstring

        assert size > 0, errstring

        errstring  = "\n\n                  In absError size of first field != size of second field \n" 
        errstring += self.abortstring

	if fieldb == None:
	    fieldb = 0.0*fielda
	else:
	    assert len(fieldb) == len(fielda), errstring

	array = fielda-fieldb

	if len(Numeric.shape(array)) > 1:
	    AB    = Numeric.diagonal(Numeric.sqrt(Numeric.dot(array,Numeric.transpose(array))))
	else:
	    AB    = Numeric.absolute(array)
	
	return AB 

    def percentError(self, first, second):
	"percentage error = |first-second|/|second| where first, second are scalars"

        errstring  = "\n\n  In percentError abs(denominator)=0 \n"
        errstring += self.abortstring
        
	assert abs(second) > 0, errstring 
	
	PE = abs(first-second)/abs(second)
	
	return PE

    def l1Error(self,fielda,fieldb=None):
	"L1 calculation of error fielda-fieldb"

        errstring  = "\n\n                  In l1Error first field is not a Numeric array \n"
        errstring += self.abortstring

	assert type(fielda) == type(Numeric.array(1)), errstring

	size = len(fielda)

        errstring  = "\n\n                  In l1Error first field size <= 0 \n" 
        errstring += self.abortstring

	assert size > 0, errstring

        errstring  = "\n\n                  In l1Error size of first field != size of second field \n" 
        errstring += self.abortstring
        
	if fieldb == None:
	    fieldb = 0.0*fielda
	else:
	    assert len(fieldb) == len(fielda), errstring

	array = fielda-fieldb


	if len(Numeric.shape(array)) > 1:
	    L1    = Numeric.sum(Numeric.sqrt(Numeric.diagonal(Numeric.dot(array,Numeric.transpose(array)))))/size
	else:
	    L1    = Numeric.sum(Numeric.absolute(array))/size

	return L1


    def l2Error(self,fielda,fieldb=None):
	"L2 measure of error fielda-fieldb"

        errstring  = "\n\n                  In l2Error first field is not a Numeric array \n"
        errstring += self.abortstring

	assert type(fielda) == type(Numeric.array(1)), errstring

	size = len(fielda)

        errstring  = "\n\n                  In l2Error first field size <= 0 \n" 
        errstring += self.abortstring

	assert size > 0, errstring

        errstring  = "\n\n                  In l2Error size of first field != size of second field \n" 
        errstring += self.abortstring
        
	if fieldb == None:
	    fieldb = 0.0*fielda
	else:
	    assert len(fieldb) == len(fielda), errstring

	array = fielda-fieldb

	if len(Numeric.shape(array)) > 1:
	    L2    = Numeric.trace((Numeric.dot(array,Numeric.transpose(array))))/size
            L2    = Numeric.sqrt(L2)
	else:
	    L2    = Numeric.sum(array**2)/size
            L2    = L2**(0.5)

	return L2

    def linfError(self,fielda,fieldb=None):
	"Linf measure of error fielda-fieldb"

        errstring  = "\n\n                  In linfError first field is not a Numeric array \n"
        errstring += self.abortstring

	assert type(fielda) == type(Numeric.array(1)), errstring

	size = len(fielda)

        errstring  = "\n\n                  In linfError first field size <= 0 \n" 
        errstring += self.abortstring

	assert size > 0, errstring

        errstring  = "\n\n                  In linfError size of first field != size of second field \n" 
        errstring += self.abortstring
        
	if fieldb == None:
	    fieldb = 0.0*fielda
	else:
	    assert len(fieldb) == len(fielda), errstring

	array = fielda-fieldb

	if len(Numeric.shape(array)) > 1:
	    Linf    = max(Numeric.sqrt(Numeric.diagonal(Numeric.dot(array,Numeric.transpose(array)))))
	else:
	    Linf    = max(Numeric.absolute(array))

	return Linf


if __name__== '__main__':
    " DataMeasure component test"

    measure = DataMeasure()
    field   = Numeric.array(range(10))

    MV      = measure.meanValue(field)
    print 'MV : %5.2f' %(MV) 

    field   = Numeric.array([[1,2,3],[4,5,6]])
    #print field
    fieldb  = Numeric.array([[1,2,3],[4,4,4]])
    #print fieldb
    fieldc  = Numeric.array([1,2,3])
    MV      = measure.meanValue(field)
    print 'MV:'
    print MV
    minv    = measure.minValue(fieldc)
    print minv
    maxv    = measure.maxValue(fieldc)
    print maxv
    totv    = measure.totalValue(field)
    print 'totv:'
    print totv
    L1      = measure.l1Error(field,fieldb)
    print L1
    Linf    = measure.linfError(field,fieldb)
    print Linf
    abs     = measure.absError(field,fieldb)
    print abs

			   
