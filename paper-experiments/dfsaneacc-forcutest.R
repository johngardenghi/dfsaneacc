library("dfsaneacc")
source("../cutest-R-interface/cutest.R")

lerstd <- file( "stdin" )
input <- readLines( lerstd, n=1 )
close( lerstd )

probname <- input[1]
problem_current <- formatC(toString(probname), format="s", width=10,
                           flag="-")

cutest_init(probname, interfaceDir="../cutest-R-interface/")

n <- cutest_getn()
x0 <- cutest_getx0()

## Print initial information to file
F0 <- cutest_evalr( x=x0 )
normF0 <- formatC( sqrt( sum( F0 * F0 ) ), format="E", digits=1, width=7, mode="double" )
istop <- formatC( 9, format="d", digits=2 )
iter <- formatC( 0, format="d", digits=9 )
fcnt <- formatC( 0, format="d", digits=9 )
time.taken <- formatC( 0.0, format="f", digits=6, width=12, mode="double" )

textfile <- file.path( "tabline.txt" )
printer <- file( textfile, "w" )
write( c( problem_current, n, istop, normF0, iter, fcnt, time.taken ),
      textfile, sep=" ", append=TRUE, ncolumns=7 )
close(printer)

## Parameters to call Accelerated DF-SANE
eps <- 1.e-06 * sqrt(n)
maxit = 1.0e+20

## Calls Accelerated DF-SANE
start.time <- Sys.time()
ret <- dfsaneacc(x=x0, evalr=cutest_evalr, nhlim=6, epsf=eps,
                 maxit=Inf, iprint=-1)
end.time <- Sys.time()

## Print final results to file
Ffinal <- cutest_evalr( x=ret$x )
normFfinal <- formatC( sqrt( sum( Ffinal * Ffinal ) ),
                      format="E", digits=1, width=7, mode="double" )
n <- formatC( n, format="d", digits=6 )
istop <- formatC( ret$istop, format="d", digits=2 )
iter <- formatC( ret$iter, format="d", digits=9 )
fcnt <- formatC( ret$fcnt, format="d", digits=9 )
time.taken <- formatC( as.numeric( difftime(end.time, start.time, units="secs") ),
                      format="f", digits=6, width=12, mode="double" )

cat( "Satisfied stopping criterion = ",istop,"\n" )
cat( "Residual Euclidian norm = ",normFfinal, "\n" )
cat( "Number of iterations = ",ret$iter,"\n" )
cat( "Number of iterations safeguarded = ",ret$itersafe, "\n" )
cat( "Number of functional evaluations = ",ret$fcnt,"\n" )
cat( "CPU time in seconds = ",time.taken,"\n" )

textfile <- file.path( "tabline.txt" )
printer <- file( textfile, "w" )
write( c( problem_current, n, istop, normFfinal, iter, fcnt, time.taken ),
      textfile, sep=" ", append=TRUE, ncolumns=7 )
close(printer)

cutest_end()
