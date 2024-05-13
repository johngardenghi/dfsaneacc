\name{dfsaneacc}
\alias{dfsaneacc}
\title{Accelerated derivative-free spectral residual method for
  nonlinear systems of equations}
\description{Accelerated derivative-free algorithm to solve nonlinear
  systems of equations.}
\usage{
dfsaneacc(x, evalr, nhlim = 6, epsf = 1e-06, maxit = Inf, iprint = -1, ...)
}
\arguments{
  \item{x}{initial estimate for the solution of the nonlinear system.}
  \item{evalr}{a function that computes the nonlinear system evaluated
    at a given point as parameter and return the evaluated value. See
    details below.}
  \item{nhlim}{an integer that determines how many previous iterates
    must be considered in the sequential secant acceleration step. The
    default is 6.}
  \item{epsf}{a real value determining the absolute convergence
    tolerance. The default is \code{1.0e-6}. See details below.}
  \item{maxit}{an integer determining the maximum number of
    iterations. The default is \code{Inf}.}
  \item{iprint}{the output level. The default is \code{-1}. See details below.}
  \item{\dots}{represents additional arguments that must be passed to \code{evalr}.}
}
\details{
  The function \code{dfsaneacc} implements sequential residual methods
  (La Cruz and Raydan 2003; La Cruz, Mart\'inez, and Raydan 2006) with
  sequential secant acceleration approach proposed by Birgin and
  Mart\'inez (2022).

  Convergence of the algorithm is declared when \eqn{\|F(x)\|_2^2 \leq
  \mbox{epsf}^2}. The default value for \code{epsf} is \code{1.0e-6}.

  The algorithm employ the function \code{evalr} to compute the value of
  the nonlinear system at a given point \code{x}. The function
  \code{evalr} must have the form \code{evalr (x, ...)}.

  The function has four output levels, based on the value of the input
  parameter \code{iprint}: \code{iprint=-1} no output is generated,
  \code{iprint=0} means basic information at every iteration,
  \code{iprint=1} adds additional information related to the
  backtracking strategy, and \code{iprint=2} adds information related to
  the computation of the acceleration step. Its default value is
  \code{iprint=-1}.
}
\value{A list with
\item{x}{the final estimate to the solution.}
\item{res}{the final nonlinear system value.}
\item{normF}{the final nonlinear system value squared L2-norm.}
\item{iter}{the total number of iterations.}
\item{fcnt}{the total number of functional evaluations.}
\item{istop}{an integer indicating the convergence type. Possible values
  are \code{0} for successful convergence (squared L2-norm of the
  residual) and \code{1} for maximum number of iterations exceeded.}
}
\author{
  \packageAuthor{dfsaneacc}

  Maintainer: \packageMaintainer{dfsaneacc}
}
\examples{
n <- 3
x0 <- rep(1/n^2, n)

expfun2 <- function(x) {
    n <- length(x)
    f <- rep(NA, n)
    f[1] <- exp(x[1]) - 1.0
    f[2:n] <- (2:n)/10.0 * (exp(x[2:n]) + x[1:n-1] - 1)
    f
}

ret <- dfsaneacc(x=x0, evalr=expfun2, nhlim=6, epsf=1.0e-6*sqrt(n),
iprint=0)
ret
}
\references{
  J. Barzilai, and J. M. Borwein (1988), Two-point step size gradient
  methods, \emph{IMA J Numerical Analysis}, 8, 141-148.

  E. G. Birgin, J. M. Mart\'inez (2022), Secant acceleration of
  sequential residual methods for solving large-scale nonlinear systems
  of equations, \emph{SIAM Journal on Numerical Analysis}, 60(6),
  3145-3180.

  W. LaCruz, and M. Raydan (2003), Nonmonotone spectral methods for
  large-scale nonlinear systems, \emph{Optimization Methods and
  Software}, 18, 583-599.

  W. LaCruz, J. M. Mart\'inez, and M. Raydan (2006), Spectral residual
  method without gradient information for solving large-scale nonlinear
  systems of equations, \emph{Mathematics of Computation}, 75,
  1429-1448.

  M. Raydan (1997), Barzilai-Borwein gradient method for large-scale
  unconstrained minimization problem, \emph{SIAM Journal on
    Optimization}, 7, 26-33.
}
\keyword{ large-scale }
