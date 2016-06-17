from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("convolve_dump_161.pyx")
)
