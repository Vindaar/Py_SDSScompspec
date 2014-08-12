from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("dust_extinction_array", ["dust_extinction_array.pyx"])]

setup(
    name = 'dust extinction array',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
