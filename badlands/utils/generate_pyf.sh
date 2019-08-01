#!/bin/bash

f2py --overwrite-signature -m flowalgo -h flowalgo.pyf flowalgo.f90
f2py --overwrite-signature -m fvframe -h fvframe.pyf fvframe.f90
f2py --overwrite-signature -m ormodel -h ormodel.pyf ormodel.f90
f2py --overwrite-signature -m pdalgo -h pdalgo.pyf pdalgo.f90
f2py --overwrite-signature -m waveseds -h waveseds.pyf waveseds.f90
