# Paper examples and experiments

This folder includes files to reproduce the paper example, given in
the beginning of the Section 3, and to reproduce the experiments
present in the Section 4.

## Pre-requisites

All the following instructions were tested in an Ubuntu 22.04.1 LTS
operating system. You must have installed `GCC` compilation suite
(including `gfortran`).

## Reproducing the paper example

In order to reproduce the first example cited in the beginning of the
Section 3 of the manuscript, just run the R script
`dfsaneacc-example.R`. It may install `dfsaneacc` CRAN package using
`install.packages` R command.

## Reproducing the paper experiments

### Running a single CUTEst problem

In order to reproduce the paper experiments regarding CUTEst problems
as well as the second example presented in the Section 3 of the paper
manuscript:

1. Install the CUTEst problem testset. Instructions are available at
[https://github.com/ralna/CUTEst/wiki/](https://github.com/ralna/CUTEst/wiki/). Be
sure to properly configure the environment varibles as described in
the installation instructions.

2. Run `dfsaneacc-forcutest.R` script, giving, as input, the name of
the CUTEst problem to run. Suitable CUTEst problems that fits in
nonlinear systems of equations follows.
```
BOOTH, CLUSTER, CUBENE, DENSCHNCNE, DENSCHNFNE,
FREURONE, GOTTFR, HIMMELBA, HIMMELBC, HIMMELBD,
HS8, HYPCIR, POWELLBS, POWELLSQ, PRICE3NE, PRICE4NE,
RSNBRNE, SINVALNE, WAYSEA1NE, WAYSEA2NE, DENSCHNDNE,
DENSCHNENE, HATFLDF, HATFLDFLNE, HELIXNE, HIMMELBE,
RECIPE, ZANGWIL3, POWELLSE, POWERSUMNE, HEART6,
HEART8, COOLHANS, MOREBVNE, OSCIPANE, TRIGON1NE,
INTEQNE, HATFLDG, HYDCAR6, METHANB8, METHANL8,
HYDCAR20, LUKSAN21, MANCINONE, QINGNE, ARGTRIG,
BROWNALE, CHANDHEU, 10FOLDTR, KSS, MSQRTA, MSQRTB,
EIGENAU, EIGENB, EIGENC, NONMSQRTNE, BROYDN3D,
BROYDNBD, BRYBNDNE, NONDIANE, SBRYBNDNE, SROSENBRNE,
SSBRYBNDNE, TQUARTICNE, OSCIGRNE, CYCLIC3, YATP1CNE,
YATP1NE, YATP2CNE, YATP2SQ
```

*IMPORTANT*: To run `dfsaneacc-forcutest.R` script, please maintain
the folder structure from this repository or change source code to
point to the correct place of `cutest.R` file, located in this
repository at `cutest-R-interface/cutest.R` path.

### Running `dfsaneacc` code implemented in Fortran

Run script
```
run-dfsaneaccma-forcutest.sh PROBLEM
```
where `PROBLEM` must be replaced by one of the CUTEst problem
described above.

This is a `bash` script, requires `gfortran` compiler installed, as
well as CUTEst installed.

### Running `NITSOL` solver

1. Download NITSOL solver at
[https://users.wpi.edu/~walker/NITSOL/](https://users.wpi.edu/~walker/NITSOL/).

2. Uncompress and compile using `make`. User may face problem with a
version of `gfortran` newer than 9 to compile NITSOL. In this case,
please install version 9 of gfortran.

3. Set the environment variable `NITSOL` to point to the NITSOL root
folder.

4. Run script
```
./run-nitsolma-forcutest.sh PROBLEM
```
where `PROBLEM` must be replaced by one of the CUTEst problem
described above.

### Running all experiments at once

1. Run script
```
runall-cutest.sh [algorithm]
```
where `[algorithm]` must be 1, 2 or 3:
* 1 to run Accelerated DF-SANE in R
* 2 to run Accelerated DF-SANE in Fortran
* 3 to run NITSOL
