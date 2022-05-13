Building Truchas with spack
======================================================================

A spack recipe for Truchas is available beginning with 22.04.1.  The
recipe will be included in the upcoming spack v0.17.3 release but
currently requires that one is using an up-to-date version of the
spack `develop` branch.  The [spack
docs](https://spack.readthedocs.io/en/latest/) are a great resource
for getting started and learning about spack.  Some helpful hints
about environments and compiler setup will be given below but they are
meant to supplement the official spack docs rather than replace.

## Advanced Spack Users
Following the compiler restrictions in [BUILDING.md](BUILDING.md) a
truchas binary can be installed via:

```shell
$ spack install truchas
```

Instead of installing the binary, truchas developers can install all the dependencies via:

```shell
$ spack install --only dependencies truchas
```
and then continue to develop/build/test truchas from the git repo.

## Newer Spack Users

### Initial spack setup

Clone the spack repo and initialize some environment variables
```bash
> git clone git@github.com:spack/spack.git
> source spack/share/spack/setup-env.sh
```
The `setup-env.sh` file is specific to bash.  If you use another shell, source the appropriate file. I use fish on my mac: `source spack/share/setup-env.fish`.
Note that this script will need to be sourced in each terminal making use of spack so adding it to a startup file (i.e. `.bashrc`) will make life easier.

>>>
## Note for internal LANL users
_If_ building on a yellow workstation, you can cut down on the download time by having spack use the LANL mirror
```bash
> spack mirror add pe-serve https://pe-serve.lanl.gov/spack-mirror
```
You will also need to add `pe-serve.lanl.gov` to the `no_proxy` list
```bash
export no_proxy="pe-serve.lanl.gov"
export NO_PROXY="pe-serve.lanl.gov"
```
>>>

Checkout the develop branch
```bash
> cd /path/to/cloned/spack
> git checkout develop
```

Spack requires the python library `clingo`.  It will try to build it on a new checkout and bootstrap itself but I've had issues with this so I usually just install it with `pip install clingo` (On my system `pip` and `python` default to python3 rather than python2).  If your pip/python can't install clingo, you need to upgrade your pip (or pip3 if pip is python2) via:
```bash
> pip3 install --upgrade pip --user
> python3 -m pip install clingo
```

At this point, the `spack spec` command should work and produce something like the following:
```bash
> spack spec zlib
Input spec
--------------------------------
zlib

Concretized
--------------------------------
zlib@1.2.12%gcc@11.2.1+optimize+pic+shared patches=0d38234 arch=linux-fedora35-broadwell

```

### Compilers
Spack can find compilers in system paths via:
```bash
$ spack compiler find
```

General compiler restrictions can be found in
[BUILDING.md](BUILDING.md).  If building truchas using intel (on
linux) it is _highly_ recommended to setup a module file (either lmod
or tcl) that loads the compiler and its dependencies into your path.
To teach spack about the intel compiler, load the module and then run
`spack compiler find`.  If all went well, there should be a file,
`~/.spack/linux/compilers.yaml`, that contains (among other things)
something like the following:

```yaml
- compiler:
    spec: intel@2021.5.0
    paths:
      cc: /opt/intel/oneapi/compiler/2022.0.2/linux/bin/intel64/icc
      cxx: /opt/intel/oneapi/compiler/2022.0.2/linux/bin/intel64/icpc
      f77: /opt/intel/oneapi/compiler/2022.0.2/linux/bin/intel64/ifort
      fc: /opt/intel/oneapi/compiler/2022.0.2/linux/bin/intel64/ifort
    flags: {}
    operating_system: fedora35
    target: x86_64
    modules: []
    environment: {}
    extra_rpaths: []
```

It is very important to change the `modules` entry to use the module
file for the intel compiler.  For example, if one loads the intel
environment via `ml oneapi`, the modules entry should be changed to
`modules: [oneapi]` for reasons given in the [spack docs](https://spack.readthedocs.io/en/latest/getting_started.html#vendor-specific-compiler-configuration).

On a mac, if gfortran doesn't work for Truchas, then one is really only
left with nag since intel on mac has lots of issues.  If one is fine
with using nag as the fortran compiler for _all_ spack packages, then
an entry like the following in `.spack/darwin/compilers.yaml` will
suffice:

```yaml
- compiler:
    spec: apple-clang@13.0.0
    paths:
      cc: /usr/bin/clang
      cxx: /usr/bin/clang++
      f77: /usr/local/bin/nagfor
      fc: /usr/local/bin/nagfor
    flags: {}
    operating_system: bigsur
    target: x86_64
    modules: []
    environment: {}
    extra_rpaths: []
```

Where the the location of `nagfor` should point to where ever it was
installed.  However, if it is desired to use gfortran by default and
only use `nagfor` for Truchas, then an environment will need to be
created and customized which will be discussed later

### External Packages

If the only thing spack knows about your system is the compilers, then spack will build _everything_ else from scratch.  This is usually the _safest_ route but can take quite a while.  To cut down on build times we can tell spack to use external packages that it finds on our system.  This does expose us to potential problems with ABI so if an error is encountered at some point due to one of our external packages, we will need to remove it from our `packages.yaml` file and rebuild (or reconcretize in spack speak).  The spack folks are aware of ABI issues when reusing packages, so, by default, only certain external packages are supported.  We can add all these packages (if found) with:
```bash
spack external find --not-buildable
```
The `--not-buildable` option ensures that spack will use our external package even if the spec is not an exact match.   At this point, spack should have generated a `~/.spack/packages.yaml` file with a bunch of entries that look like:
```yaml
packages:
  bzip2:
    externals:
    - spec: bzip2@1.0.8
      prefix: /usr
    buildable: false
  diffutils:
    externals:
    - spec: diffutils@3.7
      prefix: /usr
    buildable: false
  findutils:
    externals:
    - spec: findutils@4.8.0
      prefix: /usr
    buildable: false
  gawk:
    externals:
    - spec: gawk@5.1.0
      prefix: /usr
    buildable: false
  gmake:
    externals:
    - spec: gmake@4.3
      prefix: /usr
    buildable: false
  ...
```

### Spack Environments

[Environments](https://spack.readthedocs.io/en/latest/environments.html)
are a powerful spack feature and can be used to customize the
compilers and packages needed for Truchas.  They also work with cmake,
simplifying development.  The focus here will be on using spack to
build the truchas dependencies and using environments to tell cmake
where those dependencies are.

We begin by creating a new spack environment
```bash
> spack env create truchas-devel
==> Updating view at /home/ptb/gitrepos/spack/var/spack/environments/truchas-devel/.spack-env/view
==> Created environment 'truchas-devel' in /home/ptb/gitrepos/spack/var/spack/environments/truchas-devel
==> You can activate this environment with:
==>   spack env activate truchas-devel
```

If building with gfortran on mac or linux, most things will just work
and we can activate the environment and install the truchas
dependencies:

```bash
> spack env activate truchas-devel
> spack install --only dependencies truchas %gcc
```

Once completed, note that the spack environment has set
`CMAKE_PREFIX_PATH` to environment directory so that cmake will
automatically find all the dependencies spack just installed.  The
environment can be deactivated via `spack env deactivate`.
Deactivating the environment will remove the entry from
`CMAKE_PREFIX_PATH`.  If building the truchas binary simply remove the
`--only dependencies` qualifier from the install command.  Note the
that environment also sets `PATH` giving direct access to the truchas
executable.

Things can get tricky using intel or nag and are somewhat setup dependent.

#### Linux + Intel

If your intel module file is relatively sparse and doesn't load things
like `mkl` or `mpi` into your paths then the instructions above should
work fine.  However, by default, the intel module file loads up such
things into your path.  Spack will happily ignore these externals
unless told otherwise but this does complicate the development
environment.  One particularly bad outcome is when spack builds all
the truchas dependencies on its own `lapack` and `mpi` but
then cmake picks up on the external intel packages when one goes to
build truchas from source.  The conflicting versions seem to compile
fine but then trigger all sorts of runtime errors due to mismatching
ABI's.  To fix this, one could customize the spack environment to also
use the intel externals rather than building duplicate versions.

The environment can be customized by editing its `spack.yaml` file.  This file is pretty bare for a new environment
```bash
> cat $(spack location -e truchas-devel)/spack.yaml
# This is a Spack Environment file.
#
# It describes a set of packages to be installed, along with
# configuration settings.
spack:
  # add package specs to the `specs` list
  specs: []
  view: true
```

Having run into the issue just described above, the following customization told spack to use the intel mpi and mkl that were installed as part of the oneapi toolkits on my system:

```yaml
spack:
  # add package specs to the `specs` list
  specs: []
  view: true
  packages:
    all:
      compiler: [intel]
      providers:
        mpi: [intel-oneapi-mpi]
        blas: [intel-oneapi-mkl]
        lapack: [intel-oneapi-mkl]
    intel-oneapi-mpi:
      externals:
      - spec: intel-oneapi-mpi@2021.5.1
        prefix: /opt/intel/oneapi
    intel-oneapi-mkl:
      externals:
      - spec: intel-oneapi-mkl@2022.0.2
        prefix: /opt/intel/oneapi
```

Note that since this is in an environment's `spack.yaml` file, it
_only_ impacts this environment.  With this customization, the
linux+intel build should proceed without issue.

#### Mac + Nag

The spack recipe for openmpi doesn't work with nag.  The recipes for
various lapack/blas also do not work.  To get around this, we need to
build these externally and then modify our spack environment to use
these external packages.  In addition, due to complications building
various packages with nag, we probably want to limit our use of nag to just
the truchas spack environment.  This can be accomplished with
something like the following environment `spack.yaml`

```yaml
spack:
  compilers:
  - compiler:
      spec: apple-clang@13.0.1
      paths:
        cc: /usr/bin/clang
        cxx: /usr/bin/clang++
        f77: /usr/local/bin/nagfor
        fc: /usr/local/bin/nagfor
      flags: {}
      operating_system: bigsur
      target: x86_64
      modules: []
      environment: {}
      extra_rpaths: []  # add package specs to the `specs` list
  packages:
    all:
      compiler: [apple-clang%13.0.1]
      providers:
        mpi: [openmpi]
        lapack: [intel-oneapi-mkl]
    openmpi:
      externals:
      - spec: openmpi@4.1.2
        prefix: /usr/local/clang/13.0/nag/7.1/openmpi/4.1.2
    intel-oneapi-mkl:
      externals:
      - spec: intel-oneapi-mkl@2022.0.2
        prefix: /opt/intel/oneapi
  specs: []
  view: true
```

There are a few things to note about this.  The first is the
`apple-clang` version.  On this system it is actually version `13.0.0`
but bumping the version in this config file helps spack to pick this
clang/nag combination.  It's hacky but there doesn't seem to be a
better solution at this time.  The openmpi version was built and
installed in the location indicated.  This should be changed to suit
the machine.  The mkl being used here was the one installed via the
oneapi toolkits.
