#!/usr/bin/env python

import os

# Name of the main executable
tbin="t-linux.x86_64.intel"

# mpich root directory
mpich_root="/home/swuser/ext"

# Determine Truchas version
this_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.realpath(this_dir + "/..")
version_file = root_dir + "/version"
print("Reading version from %s" % version_file)
version = open(version_file).read().strip()
print("Version:", version)

# Name of the distribution directory / tarball
dist="truchas-%s-Linux" % version

def run(x):
    r = os.system(x)
    if r != 0:
        raise Exception("The command '%s' failed." % x)

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



run("rm -rf %s" % dist)
run("mkdir -p %s/bin" % dist)
run("mkdir -p %s/lib" % dist)

# Copy the main binary and mpich executables
run("cp %s %s/bin" % (tbin, dist))
run("cp genre %s/bin" % (dist))
run("cp vizre %s/bin" % (dist))
run("cp %s/bin/mpiexec.hydra %s/bin/mpiexec" % (mpich_root, dist))
run("cp %s/bin/hydra_pmi_proxy %s/bin/" % (mpich_root, dist))

# Copy the truchas launch script
run("cp %s/ci/truchas-binary-Linux.sh %s/bin/truchas" % (root_dir, dist))

# Copy an example
run("mkdir -p %s/examples/broken-dam" % dist)
run("cp %s/test/broken-dam/broken-dam.inp %s/examples/broken-dam/" % \
        (root_dir, dist))

# Copy a README
run("cp %s/ci/README-binary-Linux.md %s/README.md" % \
        (root_dir, dist))

# Copy Python
run("cp /home/swuser/ext/python-install/bin/python %s/bin/" % (dist))
run("cp -r /home/swuser/ext/python-install/lib/python3.9 %s/lib/" % (dist))
run("cp write-restart.py %s/bin" % (dist))
run("cp -r ../lib/python3.9/site-packages/truchas %s/lib/python3.9/site-packages/" % (dist))
run("cp %s/ci/TruchasConfigInstall-binary-Linux.py %s/lib/python3.9/site-packages/truchas/TruchasConfigInstall.py" % (root_dir, dist))
run("cp -r ../lib/python3.9/site-packages/fortran_write %s/lib/python3.9/site-packages/" % (dist))
run("cp -r ../lib/python3.9/site-packages/grid_mapping %s/lib/python3.9/site-packages/" % (dist))
run("cp -r ../lib/libfwrite.so %s/lib/" % (dist))
run("cp -r ../lib/libgridmap.so %s/lib/" % (dist))
for lib in ["libfwrite.so", "libgridmap.so"]:
    run("patchelf --set-rpath '$ORIGIN/.' %s/lib/%s" % (dist, lib))

# Copy all dependencies and set rpath properly
for b in [tbin, "genre", "vizre", "mpiexec", "hydra_pmi_proxy"]:
    copy_deps("%s/bin/%s" % (dist, b))
    os.system("patchelf --set-rpath '$ORIGIN/../lib' %s/bin/%s" % (dist, b))

# Remove libraries that we need to use from the system
for lib in ["libc", "libm", "libpthread", "libdl", "librt"]:
    run("rm %s/lib/%s.so.*" % (dist, lib))

# Create a distribution tarball
run("tar cjf %s.tar.bz2 %s" % (dist, dist))

# Print the contents of the tarball

run("tar tjf %s.tar.bz2" % dist)
