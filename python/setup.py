from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension("*", ["*.pyx"],
        include_dirs=[numpy.get_include()], libraries=["m"],
              extra_compile_args = ["-O3", "-ffast-math", "-march=native", "-fopenmp" ],
              extra_link_args=['-fopenmp'])
        #extra_compile_args=["-ffast-math",'-O3']) 
]
setup(
    #name = "My hello app",
    ext_modules = cythonize(extensions),
)