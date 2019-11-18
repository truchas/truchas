This is a 2-stage build. The first stage builds Truchas and its TPLs using
the Intel compilers mounted from the build host machine. The second stage
installs the freely distributable Intel run-time libraries and does the
final configuration of the image.

### Stage 1 (`fedora29-intel18`)

The Intel compilers need to be installed on the build host. The docker file
assumes this is `/opt/intel`. The container will be able to run the compilers
as long as they were set up to use a network license server on the build host.
If compilers are using a node-locked license the container will not be able
run the compilers, due to the container having a different MAC address than
the host. It is possible to set the MAC address that the container will use
but I have not had any luck making this work (albeit with the NAG compiler
which uses a custom license server) and it may mess up network access.

The command is

> docker build -v /opt/intel:/opt/intel:ro --file=fedora29-intel18 [options] .

### Stage 2 (`fedora29-intel18-stage2`)

You need to capture the hash of the stage 1 image, edit the docker file
and use this hash in the initial `FROM` command.

This time we do not mount the build host's Intel directory, but in its
place install the freely-distributable run time libraries into the image.
The libraries can be obtained from
https://software.intel.com/en-us/articles/redistributable-libraries-for-intel-c-and-fortran-2018-compilers-for-linux
You'll want both the C and Fortran libraries for the version that matches
the build host's compiler version.

The docker file assumes that these have been untarred into the `intel`
subdirectory of the docker file directory, with the top-level component
stripped away (i.e., untar in the docker file directory and rename the
created directory `intel`); ensure that the hierarchy of the `intel`
directory looks like that of the host's `/opt/intel` directory. The docker
file will copy this directory to `/opt/intel` in the image. In the future
it would be better to just copy the tar files and have the docker file
unpack them.

The build command is just

> docker build --file=fedora29-intel18-stage2 [options] .
