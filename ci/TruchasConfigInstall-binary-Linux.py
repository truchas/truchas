import os

this_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.realpath(this_dir + "/../../../..")

class TruchasConfig:
    mpiexec = root_dir + "/bin/mpiexec"

    bin_dir = root_dir + "/bin"
    lib_dir = root_dir + "/lib"
    lib_suffix = ".so"
    python_executable = root_dir + "/bin/python"

    truchas_executable = os.path.join(bin_dir, "t-linux.x86_64.intel")
    write_restart_executable = os.path.join(bin_dir, "write-restart.py")
    libfwrite = os.path.join(lib_dir, "libfwrite" + lib_suffix)
    libgridmap = os.path.join(lib_dir, "libgridmap" + lib_suffix)
