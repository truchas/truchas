#!/usr/bin/env python

import h5py

fp = h5py.File('fattrs.h5','r')
print fp.keys()
datasetName = 'Mineral names'
mydata = fp.get(datasetName)


print fp.attrs.keys()

for i in fp.attrs.iteritems():
	print i
	print 'Attr name is ' + i[0]
	print 'Attr value is ' + str(i[1])
# (u'GlobalIntAttr', 10)
# (u'GlobalStringAttr', 'important data')

j =	fp.attrs.get('GlobalIntAttr')
# print j
	
print mydata.value
# Output 
# ['Calcium' 'Magnesium' 'Uranium' 'Unobtainium'
#  'My Favorite Mineral in the Whole World']

# print mydata[2] 
# 'Uranium'

fp.close()

