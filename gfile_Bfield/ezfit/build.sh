#!/bin/bash
set -x # echo on
#DEBUG="-DDEBUG"
gfortran -c $DEBUG ezfit.F90 -I../PSPLINE/Pspline -L../PSPLINE/Pspline -lpspline

