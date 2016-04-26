#!/usr/bin/env python

"""
setup.py for pyBadlands

"""
from numpy.distutils.core import setup, Extension

ext_modules = [
    Extension(
        name="SFDalgo",
        sources=['pyBadlands/libUtils/SFDalgo.f90'],
    ),
    Extension(
        name="FLOWalgo",
        sources=['pyBadlands/libUtils/FLOWalgo.f90'],
    ),
    Extension(
        name="FLWnetwork",
        sources=['pyBadlands/libUtils/FLOWstack.f90', 'pyBadlands/libUtils/FLWnetwork.f90'],
    ),
    Extension(
        name="PDalgo",
        sources=['pyBadlands/libUtils/PDalgo.f90'],
    ),
    Extension(
        name="FVframe",
        sources=['pyBadlands/libUtils/FVclass.f90', 'pyBadlands/libUtils/FVframe.f90'],
    ),
    Extension(
        name="FASTloop",
        sources=['pyBadlands/libUtils/FASTloop.f90'],
    )
]

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
    packages=['pyBadlands', 'pyBadlands.flow', 'pyBadlands.forcing', 'pyBadlands.hillslope', 'pyBadlands.libUtils', 'pyBadlands.surface'],
    ext_package='pyBadlands.libUtils',
    ext_modules=ext_modules,
    scripts=[],
)
