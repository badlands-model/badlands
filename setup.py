#!/usr/bin/env python

"""
setup.py for pyBadlands
"""
from numpy.distutils.core import setup, Extension

ext_modules = []
'''
# NOTE: compilation of f90 source is not working yet. You must run `make` from the libUtils directory.
# You can install the package using 'pip install -e <checkout dir>', but a packaged installation is not likely to work.
PDstack_mod = Extension(name = "PDstack",
                         sources=['libUtils/PDstack.f90'],
                      #extra_f90_compile_args=["-ffixed-form"],
                      )
ext_modules = [
    Extension(
        name="SFDalgo",
        sources=['libUtils/SFDalgo.f90'],
        #extra_f90_compile_args=["-ffixed-form"],
    ),
    Extension(
        name="FLOWalgo",
        sources=['libUtils/FLOWalgo.f90'],
        #extra_f90_compile_args=["-ffixed-form"],
    ),
    Extension(
        name="FLWnetwork",
        sources=['libUtils/FLOWstack.f90', 'libUtils/FLWnetwork.f90'],
        #extra_f90_compile_args=["-ffixed-form"],
    ),
    Extension(
        name="PDalgo",
        # sources=['libUtils/PDstack.f90', 'libUtils/PDalgo.f90'],
        sources=['libUtils/PDalgo.f90'],
        #extra_f90_compile_args=["-ffixed-form"],
    ),
    #Extension(
        #name="FVclass",
        #sources=['libUtils/FVclass.f90'],
        ##extra_f90_compile_args=["-ffixed-form"],
    #),
    #Extension(
        #name="FVframe",
        ##sources=['libUtils/FVframe.f90'],#, 'libUtils/FVclass.f90'],
        #sources=['libUtils/FVclass.f90', 'libUtils/FVframe.f90'],
        ##extra_f90_compile_args=["-ffixed-form"],
    #),
]
'''

setup(
    name="pyBadlands",
    version="0.1",
    author="Tristan Salles",
    author_email="",
    description=("Python interface to Badlands"),
    long_description=open('README.md').read(),
    classifiers=[
        "Development Status :: 3 - Alpha",
    ],
    packages=['pyBadlands', 'pyBadlands.libUtils'],#, 'pyBadlands.raster2TIN', 'pyBadlands.FVmethod'],
    # py_modules = modules,
    ext_package='pyBadlands',
    ext_modules=ext_modules,
    scripts=[],
)

