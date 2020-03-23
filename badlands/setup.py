from setuptools import setup, find_packages
from numpy.distutils.core import setup, Extension

import glob
import subprocess
from os import path
import io

this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

sys_includes = []
subprocess.call(["make", "-C", "utils"])
sys_includes += glob.glob("utils/*.o")
sys_includes += glob.glob("utils/*.so")
sys_includes += glob.glob("utils/*.mod")

# interface for fortran code
ext1 = Extension(
    name="badlands.flowalgo",
    sources=["utils/flowalgo.pyf", "utils/flowalgo.f90"],
    extra_link_args=["utils/classfv.o"],
)

ext2 = Extension(
    name="badlands.fvframe",
    sources=["utils/fvframe.pyf", "utils/fvframe.f90"],
    extra_link_args=["utils/classfv.o"],
)

ext3 = Extension(
    name="badlands.pdalgo",
    sources=["utils/pdalgo.pyf", "utils/pdalgo.f90"],
    extra_link_args=["utils/classpd.o"],
)

ext4 = Extension(
    name="badlands.ormodel",
    sources=["utils/ormodel.pyf", "utils/ormodel.f90"],
    extra_link_args=["utils/classoro.o"],
)

ext5 = Extension(
    name="badlands.waveseds", sources=["utils/waveseds.pyf", "utils/waveseds.f90"]
)

ext6 = Extension(name="badlands.sfd", sources=["utils/sfd.pyf", "utils/sfd.c"])

if __name__ == "__main__":
    setup(
        name="badlands",
        author="Tristan Salles",
        author_email="tristan.salles@sydney.edu.au",
        url="https://github.com/badlands-model",
        version="2.0.15",
        description="Basin and Landscape Dynamics (Badlands) is a TIN-based landscape evolution model",
        long_description=long_description,
        long_description_content_type="text/markdown",
        ext_modules=[ext1, ext2, ext3, ext4, ext5, ext6],
        packages=[
            "badlands",
            "badlands.flow",
            "badlands.forcing",
            "badlands.hillslope",
            "badlands.simulation",
            "badlands.surface",
            "badlands.underland",
        ],
        package_data={"badlands": [], "badlands": sys_includes},
        data_files=[("badlands", sys_includes)],
        include_package_data=True,
        install_requires=[
            "tribad",
            "numpy>=1.15.0",
            "six>=1.11.0",
            "setuptools>=38.4.0",
            "gFlex>=1.1.0",
            "scikit-image>=0.15",
            "pandas>=0.24",
            "scipy>=1.2",
            "h5py>=2.8.0",
            "matplotlib>=3.0",
        ],
        python_requires=">=3.5",
        classifiers=[
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
        ],
    )
