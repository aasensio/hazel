from __future__ import print_function
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

# compile the fortran modules without linking
fortran_mod_comp = 'make'
system(fortran_mod_comp)

ext_modules = [Extension(# module name:
                         "pyhazel",
                         # source file:
                         ['pyhazel.pyx'],
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3'],
                         # other files to link to
                         extra_link_args=['vars.o', 'maths.o', 'allen.o', 'svd.o', 'io.o', 'SEE.o', 'rt_coef.o', 'synth.o',
														'hazel.o','singleton.o', '-lgfortran'])]

setup(name = 'pyhazel',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)

system('cp *.so ../RunPy')
