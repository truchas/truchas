Building Truchas
===============================================================================
## Quick Start Guide

### Requirements
The Truchas build system assumes a UNIX-like environment. Current development
and testing is done on 64-bit Linux with both Intel, NAG, and GNU compilers,
and macOS with GNU and NAG compilers.
* Fortran and C/C++ compilers.  The compiler executables must be in your path.
  We use and test with the following compilers.
    - Intel Fortran and C/C++:
        - versions 18.0.x, x >= 5; 19.1.0; oneAPI 2021.1 and 2021.2
        - Many versions will **not work** because of compiler bugs.
          These include 19.0.x (any x); 19.1.x, x > 0; and oneAPI 2021.0
    - NAG Fortran (with GNU C/C++):
        - version 6.2 (build 6252 or later preferred)
        - version 7.0 (build 7026 or later required)
        - most any version of GNU C/C++ should be okay
    - NAG Fortran on MacOS (with Apple Clang C/C++):
	    - NAG version 6.2 and 7.0 (build 7028 or later required)
		- Likely any version of Apple Clang should be fine
    - GNU Fortran and C/C++. Versions 9.x and 10.x appear to be working now.
      Version 11.1.0 is not working, but 11.1.1 versions dated after 5/27/2021
      are working.
* Cmake version 3.16 or later
* Standard software development tools: make, patch, perl
* Zlib development library and header files
* Python, version 3.5 or later, along with the packages h5py (version 2.6.0 or
  later) and numpy (version 1.12.0 or later)
* MPI.  The C compiler wrapper (`mpicc`, for example) must be in your path.

Truchas requires some additional libraries, but these can be built by the
third party library build step described below.

### Compiling
Compiling Truchas for the first time is usually a two-stage process.  The
first stage involves building and installing additional third party libraries
(TPL) that Truchas requires and which are not present on your system.  This
only needs to be done once.  A cmake superbuild project for this stage can be
found in the [truchas-tpl](https://gitlab.com/truchas/truchas-tpl) repository
on GitLab. This version of Truchas is tested against the "v18" bundle of TPLs;
do a `git checkout v18` after cloning the TPL repository. See its README file
for further instructions.

Once the required TPLs are installed, the procedure for building Truchas is
straightforward. You create a build directory, run cmake from that directory,
and then run make. What you choose for a build directory is irrelevant (other
than it must be different than the current directory). Here is an example:

    $ mkdir build
    $ cd build
    $ cmake -C ../config/linux-intel.cmake \
            -D CMAKE_BUILD_TYPE=Release \
            -D TRUCHAS_TPL_DIR=<truchas_tpl_dir> ..
    $ make
    $ make install

* *Don't overlook the final `..` argument on the cmake command line!*
* The `-C` argument pre-loads the cmake cache with settings from the following
  file. The `config` subdirectory contains some examples. If none of those are
  suitable, create your own, or simply define the various variables directly
  on the cmake command line (using the `-D` flag).
* `Release` directs CMake to configure an optimized build of Truchas. Another
  option is `Debug` for an unoptimized build with lots of additional runtime
  checking.
* Set the `TRUCHAS_TPL_DIR` variable to the TPL installation directory you
  used in the first stage. It must be an absolute path.
* By default Truchas will be installed into the `install` subdirectory of the
  top-level source directory. Use the `-D CMAKE_INSTALL_PREFIX=<truchas_dir>`
  cmake argument to specify a different directory.

#### Optional Portage data mapping component
Truchas provides optional support for using the Portage library to do solution
field mapping in induction heating simulations. To enable support, which is
not included by default, add `-D USE_PORTAGE=ON` to the cmake command line.
This requires that the portage library has been compiled and installed. See
the TPL superbuild project referenced above.


### Compiling on a Mac with GNU Compilers
To build Truchas on MacOS via GCC provided by [Homebrew](https://brew.sh/),
first install GCC:

```sh
$ brew install gcc
```

MacOS provides the binary `gcc`, however this is Apple Clang and not the GNU
compiler. Brew installs the GNU GCC with the version number baked into the
binaries, e.g. `gcc-10`, `g++-10`, and `gfortran-10`. These are what we will
use.

Brew's OpenMPI formula is built over Apple Clang rather than true GCC, so
OpenMPI must be built manually. Download and unpack the latest supported OpenMPI
tarball, then configure and build using the following:

```sh
$ mkdir build
$ cd build
$ ../configure CC=gcc-10 CXX=g++-10 FC=gfortran-10 --prefix=<mpi_install_dir>
$ make all
$ make install
```

Now build the Truchas TPLs and Truchas with this newly-built OpenMPI in your
`PATH`. Note CMake must be configured to use the GCC compilers rather than Apple
Clang. This is done with the `CMAKE_*_COMPILER` variables, shown below.

For the Truchas TPLs with GCC, the `linux-gcc.cmake` configuration will do just
fine:

```sh
$ mkdir build
$ cd build
$ cmake -C ../config/linux-gcc.cmake \
        -D CMAKE_C_COMPILER=gcc-10 \
        -D CMAKE_CXX_COMPILER=g++-10 \
        -D CMAKE_Fortran_COMPILER=gfortran-10 \
        -D CMAKE_INSTALL_PREFIX=<truchas_tpl_dir> \
        ..
$ make
```

Then build Truchas:

```sh
$ mkdir build
$ cd build
$ cmake -C ../config/mac-gcc.cmake \
        -D CMAKE_C_COMPILER=gcc-10 \
        -D CMAKE_CXX_COMPILER=g++-10 \
        -D CMAKE_Fortran_COMPILER=gfortran-10 \
        -D TRUCHAS_TPL_DIR=<truchas_tpl_dir> \
        -D CMAKE_BUILD_TYPE=Release \
        ..
$ make
```


### Compiling on a Mac with NAG Fortran and Apple Clang
This is a 3 step process

1. OpenMPI
2. TPL
3. Truchas

##### Notes for OpenMPI
There is an issue with failing to trigger the NAG preprocessor on the
case insensitive Mac filesystem.  This should be fixed for OpenMPI
versions > 4.0.3.  If compiling for version <= 4.0.3, the following
configuration step should suffice (from the top level source
directory):

	$ mkdir build
	$ cd build
	$ ../configure FC=nagfor FCFLAGS=-fpp --prefix=<my_openmpi_install_dir>

The `-fpp` option can be removed for later OpenMPI releases

##### Notes for TPL
Configuring cmake using the `config/linux-nag.cmake` is sufficent.  No other
special flags have been needed

##### Notes for Truchas
`PGSLib` is a subpackage of Truchas that depends heavily on the linux
shared library linking model.  The linux linking behavior can be
emulated with the `CMAKE_SHARED_LINKER_FLAGS` that are set in
`mac-nag-chk.cmake`.  With this in mind, the following allowed for
building truchas on a mac

	$ mkdir build
    $ cd build
    $ cmake -C ../config/mac-nag.cmake \
          -D CMAKE_BUILD_TYPE=Release \
	        -D TRUCHAS_TPL_DIR=<truchas_tpl_dir> ..
    $ make
    $ make install

### Testing
From the build directory run the command

    $ ctest

to run the regression test suite. On multi-core systems use the `-j<n>` option
to tell ctest how many processes it can run simultaneously; `-j8`, for example.
All tests should pass.
