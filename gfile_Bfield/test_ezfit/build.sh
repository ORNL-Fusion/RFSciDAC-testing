#!/bin/bash
set -x # echo on
EZFIT_DIR=../ezfit
PSPLINE_DIR=../PSPLINE/Pspline
gfortran -o test_ezfit test_ezfit.f90 $EZFIT_DIR/ezfit.o -I$EZFIT_DIR -I$PSPLINE_DIR -L$PSPLINE_DIR -lpspline

