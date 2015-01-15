# Compiler 
# ifort -> Intel Fortran
# gfortran -> gfortran

# The following lines are used to choose the compiler
# We have tested gfortran and ifort
# In case you want to use another compiler, use equivalent
# keywords
#  FLAGS -> is used to indicate that the preprocessor has to be invoked
#  OPTIONS -> these are general compilation flags that, for the moment, only indicate
#             that the object file is generated, without linking
# It should be easy to find the equivalent flags in your compiler

# gfortran (comment these lines if you want to use another compiler)
COMPILER = gfortran
# PYCOMPILER = gnu95
# FLAGS =
OPTIONS = -c -ffree-line-length-none -fPIC -O3
 
# ifort (comment these lines if you want to use another compiler)
# COMPILER = ifort
# PYCOMPILER = intel
# LIBS = -lifcore -lsvml -limf
# OPTIONS = -c -g -traceback -fPIC -limf
#-fPIC -O3
# -g -traceback
# 

#-O3
 

all:
	make pyhazel

pyhazel: vars.o maths.o allen.o svd.o io.o SEE.o rt_coef.o synth.o hazel.o singleton.o

clean: 
	find . -maxdepth 2 -name "*.o" -delete ; find . -maxdepth 1 -name "*.mod" -delete ;
	find . -maxdepth 1 -name "*.f90~" -delete ; find . -maxdepth 2 -name "*.a" -delete ;
	find . -maxdepth 1 -name "pyhazel.so" -delete ; find . -maxdepth 1 -name "*.pyf" -delete
	find ../ -maxdepth 2 -name "pyhazel.so" -delete ; find . -maxdepth 1 -name "*.c" -delete
	
vars.o: vars.f90
	$(COMPILER) $(OPTIONS)  vars.f90

maths.o: maths.f90 vars.o singleton.o
	$(COMPILER) $(OPTIONS)  maths.f90

allen.o: allen.f90 maths.o vars.o
	$(COMPILER) $(OPTIONS)  allen.f90
	
svd.o: svd.f90
	$(COMPILER) $(OPTIONS)  svd.f90

io.o: io.f90 maths.o vars.o
	$(COMPILER) $(OPTIONS)  io.f90

SEE.o: SEE.f90 allen.o maths.o vars.o
	$(COMPILER) $(OPTIONS)  SEE.f90

rt_coef.o: rt_coef.f90 vars.o maths.o
	$(COMPILER) $(OPTIONS)  rt_coef.f90

synth.o: synth.f90 rt_coef.o SEE.o vars.o
	$(COMPILER) $(OPTIONS)  synth.f90
	
singleton.o: singleton.f90
	$(COMPILER) $(OPTIONS)  singleton.f90
	
test.o: test.f90
	$(COMPILER) $(OPTIONS)  test.f90

hazel.o: hazel.f90 allen.o synth.o rt_coef.o SEE.o maths.o vars.o io.o
	$(COMPILER) $(OPTIONS)  hazel.f90
