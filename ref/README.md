# The Thompson-Micro-Physics reference kernel

NOTE: If you are reading this with a plain text editor, please note that this document is
formatted with Markdown syntax elements.  See https://www.markdownguide.org/cheat-sheet/
for more information.

This is the reference implementation of the `thompson` kernel.  It contains a base CPU implementation
parallelized with OpenMP.  There are also two GPU implementations that were created by adding
OpenACC directives to the CPU version.  The number of OpenMP threads is determined by the
`OMP_NUM_THREADS` environment variable, which defaults to the number of hyperthreads on the
machine if it is not set.

## Dependencies
The following packages are required for building and running this kernel:

* C compiler
* Fortran compiler (PGI is required for GPU support)
* [cmake](https://cmake.org/download/) (version >= 3.10)
* git

## Building the kernel

This kernel uses an out-of-source cmake build, meaning that the build must be done in 
directory that is not in the source tree.

### Basic build procedure (from the directory containing this file)

The basic build steps are as follows (**NOTE**: Make sure not to omit the two dots at the end
of the `cmake` step.):

```bash
$ rm -rf build ; mkdir build ; cd build
$ export CC=<name of C compiler>
$ export FC=<name of fortran compiler> 
$ cmake -DCMAKE_BUILD_TYPE=<debug | release> -DENABLE_GPU=<on | off>..
$ make VERBOSE=1
```

On many systems, the above will suffice. However, some systems will require you to help cmake
find dependencies, particularly if software depencencies are not installed in standard locations.
See below for more information.

NOTE: If GPU support is enabled, you must be on a machine that has GPU hardware and drivers as well as
CUDA installed.  Currently, GPU versions are only supported for the PGI/NVHPC compilers.

### Machines that use modules to manage software

Most HPC systems use modules to manage software.  Make sure you have loaded the versions of
the compiler and software you want to use before running the build steps above.  This will allow build
dependencies to be found properly.  For example:

```bash
$ module load pgi/19.10 cuda/10.1 cmake
```

### Machines that do not use modules to manage software

If compilers or other dependencies are not installed in a standard location where cmake can find it,
you may need to add their installation paths to the `CMAKE_PREFIX_PATH` before running the steps
above. For example:

```bash
$ export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/path/to/pgi_compiler
```

### Building on a Mac

By default, gcc points to the clang compiler on Mac.  To use the GNU compiler on Mac, depending
on how the GNU compiler was installed, you may need to specify the C compiler name as gcc-$version.
For example:

```bash
$ export CC=gcc-10
```

## Testing the kernel

First, set the OpenMP variables, including number of threads you want to use for the tests. For example:

```bash
$ export OMP_PLACES=cores
$ export OMP_PROC_BIND=close
$ export OMP_NUM_THREADS=4
```

Then, to run the test suite (from the build directory):

```bash
$ ctest
```

To run a specific test (for example):

```bash
$ ctest -R cpu_kernel
```

To run a specific test with full output to get more information about a failure (for example):

```bash
$ ctest -VV -R cpu_kernel
```

## Installation and running

To (optionally) install the built executable into exe/

```bash
$ make install
```

To run the installed executable (for example):

```bash
$ export OMP_PLACES=cores
$ export OMP_PROC_BIND=close
$ export OMP_NUM_THREADS=4
$ exe/thompson_kernel_cpu
```

## NOTES:

1. The test suite does not measure performance, but reports how long each test takes to run.
2. Detailed performance timings are printed in the stdout when the kernel runs.
3. To view kernel output, either run ctest in verbose mode, or run the kernel manually.

## Here is a list of the files and what they contain.

- `cmake/` contains compiler flag settings and cmake helper scripts
- `src/` contains the kernel source code
- `test/` contains the tests and test output
- `exe/` contains the installed executable


