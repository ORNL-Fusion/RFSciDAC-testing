#!/bin/bash
set -x # echo on
gfortran -c kind_mod.f90
gfortran -c bspline90_22.f90
gfortran -c numerics_module.f90
gfortran -c g3d_module.f90
ar rcs libutil.a *.o
