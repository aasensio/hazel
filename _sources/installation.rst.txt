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

#. ``pyGUI`` contains a Python front-end for the synthesis problem

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

            make version=serial compiler=ifort

After compiling and linking, the executable is copied to appropriate
places. There is an additional requirement for compiling the parallel version
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

           make version=mpi compiler=ifort supercomputer=local


The current supported options for the makefile are the following:

#. ``version``: serial/mpi/python

#. ``compiler``: ifort/gfortran

#. ``supercomputer``: local/teide (for the Teide Supercomputer)

If you want to add your configuration, appropriately modify the makefile. You would need to define
the location of the NetCDF include files and libraries (``NETCDF_INCLUDE`` and ``NETCDF_LIB``) 
and the correct compilation and linking flags for your compiler (``COMPILER_OPTS`` and ``COMPILER_LINKS``).
If you want to run the GUI, you first need to compile the serial
version of the code. Finally, the generated object and module files can be cleaned typing:

::

           make clean