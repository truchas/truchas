Truchas Third-Party Library Superbuild
------------------------------------------------------------------------------
This directory contains a cmake system for building most of the third party
packages needed by Truchas.  While we endeavor to make this as robust as
possible, especially for the platforms we use and test on, it is inevitable
that this will not work for some people.  In that circumstance remember that
you can build them manually (or obtain them otherwise); you don't have to
make this system work.  All that matters is that the Truchas cmake step can
find them.

### Quick Start Guide
The packages that can be built are HDF5, NetCDF, Exodus, HYPRE, Petaca, YAJL,
and SWIG.  Compressed tarfiles of their source distributions can be found in
the tarfiles subdirectory.

The basic procedure is simple (when it works). You create a build directory,
run cmake from that directory, and then run make. What you choose for a build
directory is irrelevant (other than it can't be this directory).  You will
also provide cmake command line arguments that customize the build in several
ways.  Here's an example:

    $ mkdir build
    $ cd build
    $ cmake -C ../config/linux-intel.cmake ..
    $ make

The `-C` argument preloads the cmake cache with settings from the following
file.  The `config` subdirectory contains many examples.  If none of those
are right for your situation, create your own, or simply define the various
variables directly on the cmake command line (using the `-D` flag).  By
default the packages are installed into the `install` subdirectory of the
build directory. Add `-D CMAKE_INSTALL_PREFIX=<truchas_tpl_dir>` to the cmake
command line to specify a different directory (replace `<truchas_tpl_dir>`
with the directory path; it must be an absolute path.)  Note that no `make
install` command is necessary; the packages are automatically installed by the
`make` command.  Also note that moving the libraries after they are installed
will often break things, so it is better to decide where you want them before
starting.

By default cmake will search for an existing installation of each package
and only configure a package build if it cannot find a suitable version.
You can disable this search and just build the package by setting
`-D SEARCH_FOR_<pkg>=no` on the cmake command line.  Here `<pkg>` can be
`HDF5`, `NETCDF`, `EXODUS`, `HYPRE`, `YAJL`, `PETACA`, or `SWIG`.  This is
sometimes necessary when cmake finds a package that you don't want it to use.
Also note that that when searching, cmake always looks first in the
installation directory.

Another command line variables is `ENABLE_SHARED`.  The default is `yes`.
Setting `ENABLE_SHARED` to `no`, meaning find and build static libraries,
generally doesn't work currently except on platforms that provide a "full"
set of static system libraries (Cray, for example).

### Package Configuration Notes
If you find you have to build some packages manually (or obtain pre-built
ones), here are some points to keep in mind.  Also look at the files in the
`cmake` directory to see how we are configuring them.

#### HDF5
* Only the C interface is needed.
* The high-level library (HL) is needed.
* Be sure to enable parallel HDF5.

#### NetCDF
* Use `--with-netcdf-4`
* Only need the C interface (`--disable-fortran --disable-cxx`)

#### Exodus
* We use and test with a relative old version 5.14.  There are reported
  incompatibilities with current 6.x versions.

#### HYPRE
* Only the C interface is needed (`--disable-fortran`)
* Use `--with-MPI`.
* Use `--without-fei`; it is not needed and has compilation problems.
* We "require" an older 2.6.0b version, but newer version really should work
  too, but small numerical differences are producing enough variation in the
  output to cause some regression tests to report failures, and we haven't
  had the time to evaluate them yet.

#### SWIG
* We require an old version 2
