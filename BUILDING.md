Building Truchas
===============================================================================
## Quick Start Guide

### Requirements
The Truchas build system assumes a UNIX-like environment. Current development
is done primarily on 64-bit Linux and MacOS.
* Fortran and C/C++ compilers.  The compiler executables must be in your path.
  Truchas requires a modern Fortran 2018 compiler. The following suggested
  compilers are recommended for new builds and have worked in recent ad hoc
  testing, but this is not an exhaustive compatibility list.
  - Intel oneAPI Fortran/C/C++ 2025 or later.
  - LLVM flang/Clang/C++ 22.1 or later.
  - NAG Fortran 7.2 or later on Linux (with GNU C/C++; most any version should
  be okay).
  GitLab CI currently tests with Intel oneAPI 2025.2. Some internal CI jobs use
  an additional compiler configuration to match internal computing resources;
  that configuration is not recommended for new external builds.
* CMake version 3.20.2 or later
* Standard software development tools: make, patch, perl
* Zlib development library and header files
* LAPACK library or equivalent
* Python, version 3.12 or later, along with the packages h5py (version 3.16.0
  or later) and numpy (version 2.4.4 or later). Scipy is needed for 1 test (wg).
* MPI.  The C compiler wrapper (`mpicc`, for example) must be in your path.

Truchas requires some additional libraries, but these can be built by the
third party library build step described below.

### Compiling
Compiling Truchas for the first time is usually a two-stage process.  The
first stage involves building and installing additional third party libraries
(TPL) that Truchas requires and which are not present on your system.  This
only needs to be done once.  A CMake superbuild project for this stage can be
found in the [truchas-tpl](https://gitlab.com/truchas/truchas-tpl) repository
on GitLab. This version of Truchas uses the "v26" bundle of TPLs; do a
`git checkout v26` after cloning the TPL repository. See its README file
for further instructions.

Once the required TPLs are installed, the procedure for building Truchas is
straightforward. You create a build directory, run `cmake` from that directory,
and then run `make`. What you choose for a build directory is irrelevant (other
than it must be different than the current directory). Here is an example:

    $ mkdir build
    $ cd build
    $ cmake .. -C ../config/linux-intel.cmake \
            -D CMAKE_BUILD_TYPE=Release \
            -D CMAKE_PREFIX_PATH=<truchas_tpl_dir>
    $ make
    $ make install

* The `-C` argument pre-loads the CMake cache with settings from the following
  file. The `config` subdirectory contains some examples. If none of those are
  suitable, create your own, or simply define the various variables directly
  on the cmake command line (using the `-D` flag).
* `Release` directs CMake to configure an optimized build of Truchas. Another
  option is `Debug` for an unoptimized build with lots of additional runtime
  checking.
* Set the `CMAKE_PREFIX_PATH` variable to the TPL installation directory you
  used in the first stage. It must be an absolute path.
* By default Truchas will be installed into the `install` subdirectory of the
  top-level source directory. Use the `-D CMAKE_INSTALL_PREFIX=<truchas_dir>`
  cmake argument to specify a different directory.

#### LAPACK
A system-installed version of the LAPACK library is required. This can be a
generic reference version provided by a Linux distribution, or a specialized
one such as the Intel Math Kernel Library (MKL). CMake is generally able to
locate the library, though it may require some help through the setting of
additional CMake variables on the `cmake` command line. See this CMake
[page](https://cmake.org/cmake/help/latest/module/FindLAPACK.html)
for a description of the relevant variables, and `BLA_VENDOR` in particular.

* Nothing additional should need to be specified when using a generic reference
  LAPACK library installed in a standard system location.

* To use the Intel oneAPI MKL library with the Intel oneAPI compiler, use the
  additional argument `-D BLA_VENDOR=Intel10_64lp`. The environment variable
  `MKLROOT` must also be set to the root of the MKL installation, but this is
  likely to have been done as part of configuring the environment to use the
  oneAPI compiler.

* The Intel oneAPI MKL library can also be used when compiling with NAG. In
  this case the environment variable `MKLROOT` will likely need to be set
  manually, and `-D BLA_VENDOR=Intel10_64lp_seq` must be used instead to get
  the basic sequential version.

Note that Truchas currently makes limited direct use of LAPACK and then only
for very small systems, so that use of a highly optimized library is unlikely
to yield any noticeable performance improvement.

#### Optional Portage Data Mapping Component
Truchas provides optional support for using the Portage library to do solution
field mapping in induction heating simulations. To enable support, which is
not included by default, add `-D USE_PORTAGE=ON` to the cmake command line.
This requires that the portage library has been compiled and installed. See
the TPL superbuild project referenced above.


### Compiling on a MacOS
To build Truchas on MacOS, use gfortran provided by [Homebrew](https://brew.sh/)
alongside either gcc/g++ (also from Homebrew) or apple-clang. First, install
GCC:

```sh
$ brew install gcc
```

Set up environment variables for your compiler choice, using one of the
following (assuming GCC 15 here).

```sh
$ export CC=gcc-15 CXX=g++-15 FC=gfortran-15
$ export CC=clang CXX=clang++ FC=gfortran-15
```

> Note: MacOS provides the binary `gcc`, however this is Apple Clang and not the
> GNU compiler. Brew installs the GNU GCC with the version number baked into the
> binaries, e.g. `gcc-15`, `g++-15`, and `gfortran-15`.

Brew's OpenMPI formula is built over Apple Clang rather than true GCC, so
OpenMPI must be built manually. Download and unpack the latest supported OpenMPI
tarball, then configure and build using the following:

```sh
$ mkdir build
$ cd build
$ ../configure --prefix=<mpi_install_dir>
$ make all
$ make install
```

Now build the Truchas TPLs and Truchas with this newly-built OpenMPI in your
`PATH`.

For the Truchas TPLs, use either the `macos-gcc.cmake` configuration or the
`macos-gcc-clang.cmake` configuration, depending on whether you are using GCC or
Apple Clang for the C/C++ compilers, respectively. In either case, the
configuration will grab the environment variables set earlier. For example:

```sh
$ mkdir build
$ cd build
$ cmake -C ../config/macos-gcc-clang.cmake \
        -D CMAKE_INSTALL_PREFIX=<truchas_tpl_dir> \
        ..
$ make
```

Then build Truchas:

```sh
$ mkdir build
$ cd build
$ cmake -C ../config/mac-gcc.cmake \
        -D TRUCHAS_TPL_DIR=<truchas_tpl_dir> \
        -D CMAKE_BUILD_TYPE=Release \
        ..
$ make
```

> Note: For the Truchas build, `mac-gcc.cmake` is used for both GCC and Apple
> Clang configurations.

### Testing
From the build directory run the command

    $ ctest

to run the regression test suite. On multi-core systems use the `-j<n>` option
to tell ctest how many processes it can run simultaneously; `-j8`, for example.
All tests should pass.
