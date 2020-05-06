#!/usr/bin/env python

import os

tbin="t-linux.x86_64.intel.opt-3.1.0-alpha"
dist="truchas"
os.system("rm -rf %s" % dist)
os.system("mkdir -p %s/bin" % dist)
os.system("mkdir -p %s/lib" % dist)

os.system("cp %s %s/bin" % (tbin, dist))
os.system("patchelf --set-rpath '$ORIGIN/../lib' %s/bin/%s" % (dist, tbin))

# Copy openmpi dependencies:

# Disable for now
#os.system("cp -r /home/swuser/ext/lib/* %s/lib/" % dist)

# Copy other dependencies

os.system("""\
ldd %s \
    | grep "=>" \
    | sed -e '/^[^\t]/ d' \
    | sed -e 's/\t//' \
    | sed -e 's/.*=..//' \
    | sed -e 's/ (0.*)//' \
    > ldd_output
""" % tbin)

for l in open("ldd_output").readlines():
    lib = l.strip()
    filename = os.path.basename(lib)
    print("Direct dependency %s" % filename)
    os.system("cp %s %s/lib" % (lib, dist))
    lib = "%s/lib/%s" % (dist, filename)
    os.system("patchelf --set-rpath '$ORIGIN/.' %s" % lib)
    os.system("""\
    ldd %s \
        | grep "=>" \
        | sed -e '/^[^\t]/ d' \
        | sed -e 's/\t//' \
        | sed -e 's/.*=..//' \
        | sed -e 's/ (0.*)//' \
        > ldd_output
    """ % lib)
    for m in open("ldd_output").readlines():
        lib2 = m.strip()
        filename2 = os.path.basename(lib2)
        print("    Indirect dependency %s" % filename2)
        os.system("cp %s %s/lib" % (lib2, dist))
        lib2 = "%s/lib/%s" % (dist, filename2)
        os.system("patchelf --set-rpath '$ORIGIN/.' %s" % lib2)

# Remove libraries that we need to use from the system
os.system("rm %s/lib/libc.so.*" % dist)
os.system("rm %s/lib/libm.so.*" % dist)
os.system("rm %s/lib/libpthread.so.*" % dist)
os.system("rm %s/lib/libdl.so.*" % dist)
os.system("rm %s/lib/librt.so.*" % dist)
# We will distribute this one for now:
#os.system("rm %s/lib/libgcc_s.so.*" % dist)
os.system("tar cjf %s.tar.bz2 %s" % (dist, dist))
