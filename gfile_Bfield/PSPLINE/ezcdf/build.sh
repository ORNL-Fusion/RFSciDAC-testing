#!/bin/bash
FC=gfortran
NETCDF_INC_DIR=/usr/include
set -x
$FC -c ezcdf_opncls.f90 -I$NETCDF_INC_DIR
$FC -c ezcdf_inqvar.f90 -I$NETCDF_INC_DIR
$FC -c ezcdf_GenPut.f90 -I$NETCDF_INC_DIR
$FC -c ezcdf_GenGet.f90 -I$NETCDF_INC_DIR
$FC -c ezcdf_attrib.f90 -I$NETCDF_INC_DIR
$FC -c ezcdf.f90 -I$NETCDF_INC_DIR
ar rcs libezcdf.a *.o
