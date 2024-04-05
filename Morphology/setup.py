#
# Standard imports.
#
import os
import codecs
#
from distutils.command.sdist import sdist as DistutilsSdist
from setuptools import setup, dist, find_packages
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext as _build_ext


#
# Build extensions wrapper to handle numpy includes.
#
class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        print(numpy.get_include())
        self.include_dirs.append(numpy.get_include())


#
# Setup function.
#
setup(
    name='morphology',
    description='galaxy morphology type package',
    long_description="",
    author='Julia Largett, University of Rochester',
    author_email='jlargett@u.rochester.edu',
    license='BSD 3-clause License',
    url='https://github.com/DESI-UR/DESI_SGA/tree/Julias-updates',
    version='1.0.0',

    #packages=find_packages(),

    # Requirements.
    requires=['Python (>3.7.0)'],
    #install_requires=open(path_prefix + 'requirements.txt', 'r').read().split('\n'),
    zip_safe=False,
    use_2to3=False,

    # Unit tests.
    #test_suite='tests',
    #tests_require='pytest',

    # Set up cython modules.
    setup_requires=['Cython', 'numpy'],
    ext_modules = [
          Extension('Cythonformult',
                    ['Cythonformult.pyx'], 
                    library_dirs=['m']), 
    ],

    cmdclass={'build_ext':build_ext}
)