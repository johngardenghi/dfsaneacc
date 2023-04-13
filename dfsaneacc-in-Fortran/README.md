# Accelerated DF-SANE (Fortran implementation)

This is the Fortran implementation of Accelerated DF-SANE published by
Birgin and Mart&iacute;nez (2021) with minor modifications.

## How to compile and use

### Using the Makefile

You can use the `Makefile` available in this folder in three ways:

1. *To generate the library with compiled code to Accelerated DF-SANE*. Just enter
```
make
```

2. *To generate a program to run a
[CUTEst](https://github.com/ralna/CUTEst/wiki) problem*. Just enter
```
make cutest PROBNAME=[problem]
```
where `problem` is a SIF file available in your `${MASTSIF}`
folder. For instructions in how to install and configure `CUTEst`,
follow this [link](https://github.com/ralna/CUTEst/wiki).

3. *To generate a program to run the example contained in
`src/examples/dfsaneaccma.f90`*. Just enter
```
make example
```

### Instructions for manual compilation

In order to clarify the compilation procedure to whom it may concern,
instructions on how to manually compile the code follows. We consider
`gfortran` compiler and `ar` for library generation.

#### To build the library `libdfsaneacc.a`

The first step consists in generating the library `libdfsaneacc.a` to
call the function `dfsaneacc`. To do this, move to `src/dfsaneacc`
folder and compile sources as follows:
```
gfortran -O3 -c -o qr.o qr.f90
gfortran -O3 -c -o dfsaneacc.o dfsaneacc.f90
ar rcs ./lib/libdfsaneacc.a qr.o dfsaneacc.o
```

The library `libdfsaneacc.a` will be generated in `src/dfsaneacc/lib`
folder.

#### How to compile the example file (and any other you may need)

Consider that you cloned this repository in your computer, and the
complete path to `dfsaneacc-in-Fortran` folder is `${DFSANEACC}`. To
compile the example, for instance (or any other source you written by
your own), move into `src/examples` folder, compile the `main`
source by invoking
```
gfortran -O3 -c -o dfsaneaccma.o dfsaneaccma.f90
```
and link with the library generated in the previous section by entering:
```
gfortran -O3 -L${DFSANEACC}/src/dfsaneacc/lib -o dfsaneaccma-example dfsaneaccma.o -ldfsaneacc
```

With this, the executable file `dfsaneaccma-example` will be generated
in `src/examples` folder.

#### How to compile the cutest interface

Don't remember to follow the installation and environment variables
configuration in [CUTEst](https://github.com/ralna/CUTEst/wiki) page.

The procedure is very similar to the one shown in previous
section. You just need to generate source files for a specific problem
and include the CUTEst library when generating the executable
file. So, suppose you need to run `BOOTH` problem. You need to move
into `src/interfaces/cutest` directory and

1. Generate SIF files 
```
sifdecoder -A ${MYARCH} BOOTH
```
This will generate four Fortran source files (`ELFUN.f`, `EXTER.f`,
`GROUP.f`, and `RANGE.f`) and two data files (`AUTOMAT.d` and
`OUTSDIF.d`).

2. Compile the generated files and the main program
```
gfortran -O3 -c ELFUN.o ELFUN.f
gfortran -O3 -c EXTER.o EXTER.f
gfortran -O3 -c GROUP.o GROUP.f
gfortran -O3 -c RANGE.o RANGE.f
gfortran -O3 -c dfsaneaccma.o dfsaneaccma.c
```

3. Link everything with `libdfsaneacc.a` and `libcutest.a`.
```
gfortran -O3                                                                      \
	 -L${DFSANEACCSRC}/src/dfsaneacc/lib -L${CUTEST}/objects/${MYARCH}/double \
	 -o dfsaneaccma-forcutest $^ -lcutest -ldfsaneacc
```

## References

E. G. Birgin, J. M. Mart&iacute;nez (2021), Secant acceleration of
sequential residual methods for solving large-scale nonlinear systems
of equations. *Technical Report
[arXiv:2012.13251v1](https://arxiv.org/abs/2012.13251)*.
