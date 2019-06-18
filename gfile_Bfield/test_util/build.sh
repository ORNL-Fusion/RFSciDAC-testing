#!/bin/bash
set -x # echo on
UTIL_DIR=../util-library
LAPACK_DIR=/usr/lib
gfortran -o test_util test_util.f90 -I$UTIL_DIR -L$UTIL_DIR -lutil -L$LAPACK_DIR -llapack

