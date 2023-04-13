#!/bin/bash

tlim=180

${SIFDECODE}/bin/sifdecoder -A ${MYARCH} ${1}.SIF

# COMPILE ALL SOURCE FILES
gfortran -O3 -c ELFUN.f -o ELFUN.o
gfortran -O3 -c EXTER.f -o EXTER.o
gfortran -O3 -c GROUP.f -o GROUP.o
gfortran -O3 -c RANGE.f -o RANGE.o
gfortran -O3 -c nitsolma-forcutest.f90 -o nitsolma-forcutest.o

# LINK FILES
gfortran \
    -O3 -L ${CUTEST}/objects/${MYARCH}/double -L${NITSOL}/Nitsol \
    -o nitsolma-forcutest \
    ELFUN.o EXTER.o GROUP.o RANGE.o nitsolma-forcutest.o \
    -lcutest -lnitsol -llapack -lblas

# RUN ACCELERATED DF-SANE
ulimit -St $tlim
./nitsolma-forcutest
ulimit -St unlimited

cat tabline.txt >> table-nitsol-cutest.txt

# CLEAN FILES
rm -f ELFUN.f EXTER.f GROUP.f RANGE.f \
   ELFUN.o EXTER.o GROUP.o RANGE.o nitsolma-forcutest.o \
   problem.mod AUTOMAT.d OUTSDIF.d tabline.txt nitsolma-forcutest
