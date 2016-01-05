#! /usr/bin/env python
#===============================================================================
#
#  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
#  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
#  in the LICENSE file found in the top-level directory of this distribution.
#
#===============================================================================
#
#  Usage:
#
#    xdmf-parser.py DANU_FILE.h5
#
#    Creates light data model described in file DANU_FILE.xmf of the mesh
#    associated data in DANU_FILE.h5.
#
#

import sys

danu_py_install='@Danu_Python_INSTALL_DIR@'
sys.path.append(danu_py_install)
try:
    import xdmf_conversion
except ImportError:
    print "Attempted to add %s to Python to import Danu. Failed" % (danu_py_install)
    raise

if __name__ == '__main__':
    xdmf_conversion.start()
