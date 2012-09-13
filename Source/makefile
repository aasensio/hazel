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
#COMPILER = gfortran
#FLAGS =
#OPTIONS = -c -ffree-line-length-none
 
# ifort (comment these lines if you want to use another compiler)
COMPILER = ifort
FLAGS = 
OPTIONS = -c -O3
# -g -traceback
 
FINAL_EXECUTABLE = hazel

hazel: vars.o maths.o allen.o svd.o io.o SEE.o rt_coef.o synth.o pikaia.o marquardt.o hazel.o DIRect.o \
	DIRserial.o DIRsubrout.o singleton.o lbfgsb.o
	$(COMPILER) vars.o maths.o allen.o svd.o io.o SEE.o rt_coef.o synth.o pikaia.o marquardt.o hazel.o \
		DIRect.o DIRserial.o DIRsubrout.o singleton.o lbfgsb.o -o $(FINAL_EXECUTABLE)
	cp $(FINAL_EXECUTABLE) ../Run
	cp $(FINAL_EXECUTABLE) ../Widget_Synth
	cp $(FINAL_EXECUTABLE) ../Widget_Inv

all:
	make

clean: 
	find . -maxdepth 2 -name "*.o" -delete ; find . -maxdepth 1 -name "*.mod" -delete ;
	find . -maxdepth 1 -name "*.f90~" -delete ; find . -maxdepth 2 -name "*.a" -delete ;
	find . -maxdepth 1 -name "hazel" -delete
	find ../ -maxdepth 2 -name "hazel" -delete
	
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

pikaia.o: pikaia.f
	$(COMPILER) $(OPTIONS)  pikaia.f

marquardt.o: marquardt.f90 maths.o svd.o synth.o rt_coef.o SEE.o vars.o lbfgsb.o
	$(COMPILER) $(OPTIONS)  marquardt.f90
	
DIRect.o: DIRect.f
	$(COMPILER) $(OPTIONS)  DIRect.f
	
DIRserial.o: DIRserial.f
	$(COMPILER) $(OPTIONS)  DIRserial.f
	
DIRsubrout.o: DIRsubrout.f
	$(COMPILER) $(OPTIONS)  DIRsubrout.f

lbfgsb.o: lbfgsb.f90
	$(COMPILER) $(OPTIONS) lbfgsb.f90
# fft.o: fft.f90
# 	$(COMPILER) $(OPTIONS)  fft.f90

singleton.o: singleton.f90
	$(COMPILER) $(OPTIONS)  singleton.f90
	
hazel.o: hazel.f90 allen.o synth.o marquardt.o rt_coef.o SEE.o maths.o vars.o
	$(COMPILER) $(OPTIONS)  hazel.f90
	
