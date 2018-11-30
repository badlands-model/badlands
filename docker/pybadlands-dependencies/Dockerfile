# Pull base image.
FROM ubuntu:trusty

MAINTAINER Ben Mather

RUN apt-get update -y

RUN apt-get install -y git python3-pip python3-dev libzmq3 libzmq3-dev pkg-config libfreetype6-dev libpng3 libopenmpi-dev openmpi-bin libhdf5-dev liblapack-dev llvm-3.6 libedit-dev gfortran

RUN python3 -m pip -U setuptools
RUN python3 -m pip -U pip  # fixes AssertionError in Ubuntu pip
RUN python3 -m pip enum34
RUN LLVM_CONFIG=llvm-config-3.6 python3 -m pip llvmlite==0.8.0
RUN python3 -m pip jupyter markupsafe zmq singledispatch backports_abc certifi jsonschema ipyparallel path.py matplotlib mpi4py==1.3.1 git+https://github.com/badlands-model/triangle pandas plotly
RUN apt-get install -y libnetcdf-dev python-mpltoolkits.basemap
RUN python3 -m pip Cython==0.20
RUN python3 -m pip h5py
RUN python3 -m pip scipy
RUN python3 -m pip numpy
RUN python3 -m pip numba==0.23.1 ez_setup
RUN python3 -m pip gFlex
RUN python3 -m pip netcdf4==1.3.1
RUN python3 -m pip colorlover
RUN python3 -m pip cmocean
RUN python3 -m pip scikit-fuzzy
RUN python3 -m pip pyevtk