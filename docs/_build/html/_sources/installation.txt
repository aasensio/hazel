Installation
============

The Hazel distribution is organized as follows:

#. ``src`` contains the Fortran 90 sources and a makefile that can be
   used to build the binary file.

#. ``run`` contains a directory tree structure with the appropriate
   configuration files to run the code in command line mode.

#. ``runMPI`` contains a directory tree structure with the appropriate
   configuration files to run the parallel version of the code.

#. ``synthGUI`` contains all the files that are needed to run the IDL
   front-end for the synthesis problem.

#. ``invGUI`` contains all the files that are needed to run the IDL
   front-end for the inversion problem.

#. ``IDL_routines`` contains some IDL routines that are needed by the
   front-ends.

#. ``pyRoutines`` contains some Python routines that can be useful

#. ``docs`` contains this manual.

Hazel has a source directory that can generate both the serial and parallel
versions of the code. It has been tested on Linux platforms using the Intel Fortran
Compiler (``ifort``) and the free GFortran compiler. The compilation is performed
by means of a makefile. We provide several examples that allow us to compile
the code in several machines. As an example, we focus on ``makefile.Intel'', that
is appropriate for compiling the code in standard Intel machines. It can be run
by providing the version that will be compiled, together with the compiler. By
default, it will compile the serial version with the Intel Fortran compiler:

::
      make -f makefile.Intel version=serial compiler=ifort

After compiling and linking, the executable is copied to appropriate
places. If you want to run the GUI, you first need to compile the serial
version of the code. Finally, the generated object and module files can be cleaned typing:

::

           make -f makefile.Intel clean

There is an additional requirement for compiling the parallel version
of the code. The compilation depends on the precompiled library NetCDF for reading and writing
output files. NetCDF is a standard for platform independent binary files
that you need to have installed in your system.  The location of the NetCDF
libraries have to be provided inside the makefile. The variables ``NETCDF_INCLUDE`` and
``NETCDF_LIB`` have to point to the ``include`` and ``lib`` directories
of the NetCDF distribution.

The code makes use of the MPI package for parallelization, so it has to
be installed on your system. In order to obtain the executable file (for
instance for the Intel compiler), just type:

::

           make -f makefile.Intel version=mpi compiler=ifort

The code admits up to three command line parameters:

-  Filename with the main configuration file.

-  Starting pixel of the inversion. This is used if you want to rerun
   the inversion of some pixels.

-  Final pixel of the inversion. This is used if you want to rerun the
   inversion of some pixels.