#!/usr/bin/env python

import os

tbin="t-linux.x86_64.intel.opt-3.1.0-alpha"
dist="truchas"
os.system("rm -rf %s" % dist)
os.system("mkdir -p %s/bin" % dist)
os.system("mkdir -p %s/lib" % dist)

os.system("cp %s %s/bin" % (tbin, dist))
os.system("patchelf --set-rpath '$ORIGIN/../lib' %s/bin/%s" % (dist, tbin))

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

os.system("tar cjf truchas.tar.bz2 truchas")
