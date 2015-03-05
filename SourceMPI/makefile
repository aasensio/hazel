# makefile for the MPI version of Hazel (P-Hazel)

#******************
# NET-CDF
#******************
# If nf-config exists in your path
# NETCDF_INCLUDE = $(shell nf-config --fflags)
# NETCDF_LIB = $(shell nf-config --flibs)

# Otherwise, hardwire the path to the include files
NETCDF_INCLUDE = /net/nas4/fis/aasensio/scratch/LIBS/NETCDF/include
NETCDF_LIB = /net/nas4/fis/aasensio/scratch/LIBS/NETCDF/lib

#******************
# COMPILER (MPI)
#******************
# COMPILER = mpif90

# In my desktop computer
COMPILER = /net/nas4/fis/aasensio/scratch/MPICH/bin/mpif90
COMPILER_LINKS = -L/net/nas4/fis/aasensio/scratch/MPICH/lib

# ------------------ CHANGE THIS ACCORDING TO YOUR COMPILER ------------------
# If mpif90 calls gfortran, uncomment the following line
# COMPILER_OPTS = -c -O3 -ffree-line-length-none -I/net/nas4/fis/aasensio/scratch/MPICH/include

# ------------------ CHANGE THIS ACCORDING TO YOUR COMPILER ------------------
#  If mpif90 calls ifort, uncomment the following line
COMPILER_OPTS = -c -O3 -ipo -I/net/nas4/fis/aasensio/scratch/MPICH/include

#******************
# Hazel
#******************
LIBS =  -L$(NETCDF_LIB) -lnetcdf -lcurl 
FINAL_EXECUTABLE = phazel

OBJECTS = vars.o maths.o allen.o svd.o io.o SEE.o rt_coef.o synth.o marquardt.o hazel.o DIRect.o \
	DIRserial.o DIRsubrout.o singleton.o mpi_routines.o inversion.o
	
hazel: $(OBJECTS)
	$(COMPILER) $(COMPILER_LINKS) $(OBJECTS) $(LIBS) -I$(NETCDF_INCLUDE) -o $(FINAL_EXECUTABLE)
	cp $(FINAL_EXECUTABLE) ../RunMPI

all:
	make hazel
		
clean: 
	find . -maxdepth 2 -name "*.o" -delete ; find . -maxdepth 1 -name "*.mod" -delete ;
	find . -maxdepth 1 -name "*.f90~" -delete ; find . -maxdepth 2 -name "*.a" -delete ;
	find . -maxdepth 1 -name "phazel" -delete
	find ../ -maxdepth 2 -name "phazel" -delete
	
vars.o: vars.f90
	$(COMPILER) $(COMPILER_OPTS)  vars.f90

maths.o: maths.f90 vars.o singleton.o
	$(COMPILER) $(COMPILER_OPTS)  maths.f90

allen.o: allen.f90 maths.o vars.o
	$(COMPILER) $(COMPILER_OPTS)  allen.f90
	
svd.o: svd.f90
	$(COMPILER) $(COMPILER_OPTS)  svd.f90

io.o: io.f90 maths.o vars.o
	$(COMPILER) $(COMPILER_OPTS) -I$(NETCDF_INCLUDE) io.f90

SEE.o: SEE.f90 allen.o maths.o vars.o
	$(COMPILER) $(COMPILER_OPTS)  SEE.f90

rt_coef.o: rt_coef.f90 vars.o maths.o
	$(COMPILER) $(COMPILER_OPTS)  rt_coef.f90

synth.o: synth.f90 rt_coef.o SEE.o vars.o
	$(COMPILER) $(COMPILER_OPTS)  synth.f90

marquardt.o: marquardt.f90 maths.o svd.o synth.o rt_coef.o SEE.o vars.o
	$(COMPILER) $(COMPILER_OPTS) -I$(NETCDF_INCLUDE) marquardt.f90
	
DIRect.o: DIRect.f90
	$(COMPILER) $(COMPILER_OPTS)  DIRect.f90
	
DIRserial.o: DIRserial.f90
	$(COMPILER) $(COMPILER_OPTS)  DIRserial.f90
	
DIRsubrout.o: DIRsubrout.f90
	$(COMPILER) $(COMPILER_OPTS)  DIRsubrout.f90

singleton.o: singleton.f90
	$(COMPILER) $(COMPILER_OPTS)  singleton.f90

mpi_routines.o: mpi_routines.f90
	$(COMPILER) $(COMPILER_OPTS) mpi_routines.f90

inversion.o: inversion.f90 marquardt.o
	$(COMPILER) $(COMPILER_OPTS) -I$(NETCDF_INCLUDE) inversion.f90
	
hazel.o: hazel.f90 allen.o synth.o rt_coef.o SEE.o maths.o vars.o mpi_routines.o inversion.o io.o
	$(COMPILER) $(COMPILER_OPTS) -I$(NETCDF_INCLUDE) hazel.f90
	
