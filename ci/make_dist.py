#!/usr/bin/env python

import os

# Name of the main executable
tbin="t-linux.x86_64.intel"

# mpich root directory
mpich_root="/home/swuser/ext"

# Determine Truchas version
print("Determining version")
r = os.system("git describe --tags --dirty > version")
if r != 0:
    raise Exception("`git describe` failed")
version = open("version").read().strip()
print("Version:", version)

# Name of the distribution directory / tarball
dist="truchas-%s-Linux" % version


def copy_deps(binary_executable):
    os.system("""\
    ldd %s \
        | grep "=>" \
        | sed -e '/^[^\t]/ d' \
        | sed -e 's/\t//' \
        | sed -e 's/.*=..//' \
        | sed -e 's/ (0.*)//' \
        > ldd_output
    """ % binary_executable)

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



os.system("rm -rf %s" % dist)
os.system("mkdir -p %s/bin" % dist)
os.system("mkdir -p %s/lib" % dist)

# Copy the main binary and mpich executables
os.system("cp %s %s/bin" % (tbin, dist))
os.system("cp genre %s/bin" % (dist))
os.system("cp vizre %s/bin" % (dist))
os.system("cp %s/bin/mpiexec.hydra %s/bin/mpiexec" % (mpich_root, dist))
os.system("cp %s/bin/hydra_pmi_proxy %s/bin/" % (mpich_root, dist))

# Copy all dependencies and set rpath properly
for b in [tbin, "genre", "vizre", "mpiexec", "hydra_pmi_proxy"]:
    copy_deps("%s/bin/%s" % (dist, b))
    os.system("patchelf --set-rpath '$ORIGIN/../lib' %s/bin/%s" % (dist, b))

# Remove libraries that we need to use from the system
for lib in ["libc", "libm", "libpthread", "libdl", "librt"]:
    os.system("rm %s/lib/%s.so.*" % (dist, lib))

# Create a distribution tarball
os.system("tar cjf %s.tar.bz2 %s" % (dist, dist))

# Print the contents of the tarball

os.system("tar tjf %s.tar.bz2" % dist)
