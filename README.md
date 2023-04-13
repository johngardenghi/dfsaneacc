# Accelerated Sequential Residual Method for Nonlinear Systems

This is `dfsaneacc`, an [R](https://www.r-project.org/) package for
solving large-scale nonlinear systems of equations using
derivative-free accelerated sequential residual method.

This repository includes all source files related to the publication
Birgin, Gardenghi, Marcondes and Mart&iacute;nez (2021), namely:

* the folder `cutest-R-interface` includes an interface to use
  [CUTEst](https://github.com/ralna/CUTEst/wiki) nonlinear systems
  without derivatives in R.

* the folder `dfsaneacc-R-package` includes the R package
  `dfsaneacc`. For more details on installation in R, see the
  following section.

* the folder `dfsaneacc-in-Fortran` includes the source of the Fortran
  implementation of Accelerated DF-SANE, made by Birgin and
  Mart&iacute;nez (2021). See the README inside this folder for more
  details.

* the folder `paper-experiments` includes all the source files
  necessary to reproduce the results in the publication (Birgin,
  Gardenghi, Marcondes, Mart&iacute;nez 2021).

## How to install and use the package in R

### How to install

Download `dfsaneacc_1.0.tar.gz` inside the `dfsaneacc-R-package`
folder. Then, unpack this and type in your terminal
```
R CMD INSTALL dfsaneacc
```

### How to use dfsaneacc

You must call the routine
```
dfsaneacc (x, evalr, nhlim=6, epsf=1.e-06, maxit=Inf, iprint=-1, ...)
```
where:

* `x` is an initial estimate for the solution of the nonlinear system.

* `evalr` is a function that computes the nonlinear system evaluated
  at a given point as parameter and return the evaluated value (i.e, a
  vector of the same size of the input). See details below.

* `nhlim` is an integer that determines how many previous iterates
  must be considered in the sequential secant acceleration step. The
  default is 6.

* `epsf` is a real value determining the absolute convergence
  tolerance. The default is `1.0e-6`. See details below.

* `maxit` is an integer determining the maximum number of
  iterations. The default is `Inf`.

* `iprint` is the output level. The default is `-1`. See details
  below.

* `...` represents the additional arguments that must be passed to
  `evalr`.

This routine returns:

* `x` is the final estimate to the solution.

* `res` is the final nonlinear system value.

* `normF` is the final nonlinear system value squared L2-norm.

* `iter` is the total number of iterations.

* `fcnt` is the total number of functional evaluations.

* `istop` is an integer indicating the convergence type. Possible
  values are `0` for successful convergence (L2-norm of the
  residual) and `1` for maximum number of iterations exceeded.

### Details

Convergence of the algorithm is declared when the L2-norm is less or
equal than `epsf`. The default value for `epsf` is `1.0e-6`.

The algorithm employs the function `evalr` to compute the value of the
nonlinear system at a given point `x`. The function `evalr` must have
the form `evalr (x, ...)`.

The function has four output levels, based on the value of the input
parameter `iprint`: `iprint=-1` no output is generated, `iprint=0`
means basic information at every iteration, `iprint=1` adds additional
information related to the backtracking strategy, and `iprint=2` adds
information related to the computation of the acceleration step. Its
default value is `iprint=-1`.

### Example

To solve the Exponential Function 2 (see La Cruz and Raydan, 2003):

```
n <- 3
x0 <- rep(1/n^2, n)

expfun2 <- function(x) {
    n <- length(x)
    f <- rep(NA, n)
    f[1] <- exp(x[1]) - 1.0
    f[2:n] <- (2:n)/10.0 * (exp(x[2:n]) + x[1:n-1] - 1)
    f
}

ret <- dfsaneacc(x=x0, evalr=expfun2, nhlim=6, epsf=1.0e-6*sqrt(n), iprint=0)
ret
```

## The CUTEst interface

This interface consists in a set of R routines to evaluate the
objective function of the nonlinear system problems from CUTEst
testset.

It has an initialization routine that aims to generate a dynamic
library `R-problem-forcutest.so` to be loaded in R using `dyn.load`
routine. This initialization routine tries to generate `ELFUN.f`,
`EXTER.f`, `GROUP.f`, and `RANGE.f` Fortran source files using
`sifdecoder` from CUTEst and then to compile them together with
`cutest.c` interface to generate the dynamic library.

All this machinery was tested in an Ubuntu Linux with gfortran and
gcc. If your environment is different, you can decode your problem
manually using `sifdecoder`, compile the `.f` files with `cutest.c`
interface using your compiler and generate a dynamic library called
`R-problem-forcutest.so`. Then, you can call `cutest_init` with
optional logical variables `sifdecode` and `compile` set to
`FALSE`.

See further details in the next section.

### How to use the interface

The interface consists in five routines:

* `cutest_init(problem, interfaceDir=".", sifdecode=TRUE,
  compile=TRUE)` receives an `string` with the problem name, generates
  a dynamic library and loads it into R. This routine depends on a
  Fortran and C compiler. This was tested with GCC 9 on a Ubuntu Linux
  20.10 operating system. Optional variables are:

    * `interfaceDir` to point out the path where `cutest.R` and
      `cutest.c` files are located. If you are using these files in
      the same directory of your script, the default value is ok.

    * `sifdecode` is a logical variable indicating whether
      `sifdecoder` from CUTEst must be called to decode the SIF
      problem file.

    * `compile` is a logical variable indicating whether source file
      (`cutest.c` and files generated by `sifdecoder`) need to be
      compiled to generate the dynamic library.

* `cutest_end()` finalizes all global variables and allocated structures.

* `cutest_getn()` gets the dimension of the system.

* `cutest_getx0()` gets the initial estimate to the solution.

* `cutest_evalr(x)` computes the residual of the nonlinear system at a
  given point `x`.

This interface works only with the 70 problems from the set of square
nonlinear systems of equations.

For example, to solve `BOOTH` problem using this interface:
```
cutest_init("BOOTH")

n <- cutest_getn()
x0 <- cutest_getx0()
ret <- dfsaneacc(x=x0, evalr=cutest_evalr, nhlim=6,
                 epsf=1.0e-6*sqrt(n),iprint=0)

cutest_end()
```

You can find this example in `cutest/example.R` file.

## About the tests from the publication

### Requirements

To reproduce the results, one must have:

* a `bash` terminal.

* updated versions `gfortran`, `gcc`, `lapack`, `blas`, and `R` installed.

* [NITSOL](https://github.com/cmacmackin/nitsol) cloned and
  installed. After installation, an environment variable `NITSOL` must
  be set with the folder where the NITSOL was cloned.

* [CUTEST](https://github.com/ralna/CUTEst/wiki) installed and all the
  environment variables pointed out in the package installation
  instructions properly set.

* `dfsaneacc` package installed in R (see details above).

### Source included

The folder `paper-experiments` includes all necessary code to
reproduce the results from (Birgin, Gardenghi, Marcondes,
Mart&iacute;nez 2021), namely:

* to run **Accelerated DF-SANE (Fortran)**

    * the script `run-dfsaneaccma-forcutest.sh` compiles and run
      Accelerated DF-SANE (Fortran) for a given nonlinear system from
      CUTEst testset. To use it, enter:
      ```
      ./run-dfsaneaccma-forcutest.sh PROBLEM-NAME
      ```

* to run **Accelerated DF-SANE (R)**

    * `dfsaneacc-forcutest.R` is an R script to run a problem from
      CUTEst using the CUTEst interface developed (see details
      above). All codes are added in this script considering the
      directory organization in this repository. To use it, enter:
      ```
      Rscript dfsaneacc-forcutest.R
      ```
      It expects to enter, after this call, the CUTEst `PROBLEM-NAME`.

* to run **NITSOL (Fortran)**

    * `nitsolma-forcutest.f90` is the source to run NITSOL with CUTEst
      nonlinear systems of equations, without derivatives. The script
      `run-nitsolma-forcutest.sh` compiles and runs NITSOL by entering:
      ```
      ./run-nitsolma-forcutest.sh PROBLEM-NAME
      ```

The script `runall-cutest.sh` compiles and runs one of the three above
in all the 70 nonlinear systems from CUTEst testset. Its usage is:
```
./runall-cutest.sh [algorithm]
[algorithm] must be 1, 2, or 3:
   1  to run Accelerated DF-SANE in R
   2  to run Accelerated DF-SANE in Fortran
   3  to run NITSOL
```

## References

E. G. Birgin, J. L. Gardenghi, D. S. Marcondes, J. M. Mart&iacute;nez (2021),
Accelerated derivative-free spectral residual method for nonlinear
systems of equations, *Technical Report 
[arXiv:2104.13447](https://arxiv.org/abs/2104.13447)*.

E. G. Birgin, J. M. Mart&iacute;nez (2021), Secant acceleration of sequential
residual methods for solving large-scale nonlinear systems of
equations. *Technical Report
[arXiv:2012.13251v1](https://arxiv.org/abs/2012.13251)*.

W. La Cruz, and M. Raydan (2003), Nonmonotone spectral methods for
large-scale nonlinear systems, *Optimization Methods and Software*,
18, 583-599. [DOI
10.1080/10556780310001610493](https://doi.org/10.1080/10556780310001610493)

W. La Cruz, J. M. Mart&iacute;nez, and M. Raydan (2006), Spectral residual
method without gradient information for solving large-scale nonlinear
systems of equations, *Mathematics of Computation*, 75,
1429-1448. [DOI
10.1090/S0025-5718-06-01840-0](https://doi.org/10.1090/S0025-5718-06-01840-0)
