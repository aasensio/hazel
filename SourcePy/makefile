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
# COMPILER = gfortran
# PYCOMPILER = gnu95
# FLAGS =
# OPTIONS = -c -ffree-line-length-none -fPIC
 
# ifort (comment these lines if you want to use another compiler)
COMPILER = ifort
PYCOMPILER = intelem
LIBS = -lifcore -lsvml
OPTIONS = -c -fPIC -O3
# -g -traceback
# 

#-O3
 
FINAL_EXECUTABLE = hazel

all:
# 	make hazel
	make pyhazel
	
hazel: vars.o maths.o allen.o svd.o io.o SEE.o rt_coef.o synth.o hazel.o singleton.o test.o
	$(COMPILER) -l$(LIBS) vars.o maths.o allen.o svd.o io.o SEE.o rt_coef.o synth.o hazel.o singleton.o test.o -o $(FINAL_EXECUTABLE)
	/usr/pkg/python/Python-2.7.4/bin/f2py -lifcore --fcompiler=$(PYCOMPILER) -c hazel.pyf *.o
	cp $(FINAL_EXECUTABLE) ../RunPy
	
pyhazel: vars.o maths.o allen.o svd.o io.o SEE.o rt_coef.o synth.o hazel.o singleton.o	
	/usr/pkg/python/Python-2.7.4/bin/f2py --fcompiler=intelem $(LIBS) -c hazel.pyf vars.o maths.o allen.o svd.o io.o SEE.o rt_coef.o synth.o hazel.o singleton.o	
	cp hazel.so ../RunPy

clean: 
	find . -maxdepth 2 -name "*.o" -delete ; find . -maxdepth 1 -name "*.mod" -delete ;
	find . -maxdepth 1 -name "*.f90~" -delete ; find . -maxdepth 2 -name "*.a" -delete ;
	find . -maxdepth 1 -name "hazel.so" -delete
	find ../ -maxdepth 2 -name "hazel.so" -delete
	
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
	/usr/pkg/python/Python-2.7.4/bin/f2py --overwrite-signature -m hazel -h hazel.pyf hazel.f90
	
