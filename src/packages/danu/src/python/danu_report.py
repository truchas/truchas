#!/usr/bin/env python
#===============================================================================
#
#  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
#  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
#  in the LICENSE file found in the top-level directory of this distribution.
#
#===============================================================================

# ---------------------------------------------------------------------------- #
# Danu File Reporter
# ---------------------------------------------------------------------------- #
# Import
import os, sys
from optparse import OptionParser

try:
  import Danu
except ImportError:
  danu_mod_install='@Danu_Python_INSTALL_DIR@'
  sys.path.append(danu_mod_install)
  try:
    import Danu
  except ImportError:
    print 'Faile to locate Danu module'
    raise


# Read in the input
parser=OptionParser(usage="usage: %prog [options] FILE", version="1.0")

parser.add_option("-l", "--list", action="store_true", dest="quick_list",
                  help="print out the meshes and simulations found in the file")

parser.add_option("-m", "--mesh", dest="mesh",
                  help="print out information stored in the MESH group", metavar="MESH")

parser.add_option("-s", "--simulation", dest="sim",
                  help="print out information stored in the SIM group", metavar="SIM")

(options,args) = parser.parse_args()
try:
  output_file=args[-1]
except IndexError:
  parser.error("Failed to define output file")


# Open the file
try:
  fh=Danu.Output(output_file,'r')
except:
  print 'Failed to open %s' % output_file
  raise






