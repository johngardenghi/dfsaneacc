#!/bin/bash

tlim=180

make -C ../dfsaneacc-in-Fortran cutest PROBNAME=$1

# RUN ACCELERATED DF-SANE
mv ../dfsaneacc-in-Fortran/bin/cutest/OUTSDIF.d .
ulimit -St $tlim
../dfsaneacc-in-Fortran/bin/cutest/dfsaneaccma-forcutest
ulimit -St unlimited

cat tabline.txt >> table-dfsaneacc-fortran-cutest.txt
rm -f tabline.txt OUTSDIF.d

# CLEAN FILES
make -C ../dfsaneacc-in-Fortran cutest-distclean
