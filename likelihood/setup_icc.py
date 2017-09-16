import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

# Set an extra C compiler argument based on using either Intel or gcc compilers

cc_compiler_arg = "-axCORE-AVX2,AVX,SSE4.2"

extensions = [
    Extension("*", ["*.pyx"],
        include_dirs=[numpy.get_include()], libraries=["m"],
              extra_compile_args = ["-O3", "-ffast-math", cc_compiler_arg, "-fopenmp" ],
              extra_link_args=['-fopenmp'])
]
setup(
    ext_modules = cythonize(extensions),
)
