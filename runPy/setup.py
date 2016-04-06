from __future__ import print_function
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

# compile the fortran modules without linking
fortran_mod_comp = 'cd ../src ; make -f makefile.Intel version=python'
system(fortran_mod_comp)

ext_modules = [Extension(# module name:
                         "pyhazel",
                         # source file:
                         ['pyhazel.pyx'],
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3'],
                         # other files to link to
                         extra_link_args=['../src/vars.o', '../src/maths.o', '../src/allen.o', '../src/svd.o', '../src/io_py.o', '../src/SEE.o', '../src/rt_coef.o', '../src/synth.o',
														'../src/hazel_py.o','../src/singleton.o', '-lgfortran', '-L/usr/local/gfortran/lib'])]

setup(name = 'pyhazel',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)

system('cp *.so ../pyGUI')
