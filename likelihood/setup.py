import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

# Get CC compiler from environment variable, if none then set to gcc
cc_compiler = os.getenv('CC', 'gcc')

# Set an extra C compiler argument based on using either Intel or gcc compilers
if cc_compiler == 'icc':
	cc_compiler_arg = "-axCORE-AVX2,AVX,SSE4.2"
else:
	cc_compiler_arg = "-march=native"

extensions = [
    Extension("*", ["*.pyx"],
        include_dirs=[numpy.get_include()], libraries=["m"],
              extra_compile_args = ["-O3", "-ffast-math", cc_compiler_arg, "-fopenmp" ],
              extra_link_args=['-fopenmp'])
]
setup(
    ext_modules = cythonize(extensions),
)
