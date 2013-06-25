#!/usr/bin/env python

# ############################################################################ #
#
# Simple Danu Python script
#   
# 

import os

import Danu

# Create an output file with APPEND mode
fh = Danu.Output('dummy-file.h5', 'a')


# Create some attributes for this file
fh.set_attribute('An double value', 1.0)
fh.set_attribute('An integer value', 1)
fh.set_attribute('Created by', 'Me')

# Print out the attributes for the file
print 'Attributes in this file'
print fh.attributes()

# Now let's add a TET mesh
mesh = fh.add_unstruct_mesh('Mesh 1', Danu.TET_ELEM)
mesh.set_attribute('Notes', 'A TET mesh created in the Python iface')

# All done! Close the file
fh.close()

exit()
