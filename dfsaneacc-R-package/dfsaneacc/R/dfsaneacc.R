dfsaneacc <- function (x, evalr, nhlim=6, epsf=1.e-06, maxit=Inf,
                       iprint=-1, ...) {
    ## Set parameters
    macheps  <- .Machine$double.eps # Machine precision
    acc      <- TRUE                # Whether employ acceleration
    hsmall   <- 1.0e-01             # Used to prepare extra point to enter in history
    hlarge   <- 1.0e-01             # Used to rebuild history
    gamma    <- 1.e-04              # Constant used in descence criterion of backtracking
    lamsgmin <- sqrt( macheps )     # Lower bound on the spectral step
    lamsgmax <- 1.0 / lamsgmin      # Upper bound on the spectral step
    M        <- 10                  # Length of prior funcional values to consider in backtracking
    ndiis    <- 1                   # Interval of iterations to employ acceleration
    sigma1   <- 0.1                 # Lower bound for step size in backtracking
    sigma2   <- 0.5                 # Upper bound for step size in backtracking
    tabline  <- FALSE               # DEBUG mode: print a file with partial statistics

    n <- as.integer( length(x) )
    istop <- 9

    mlim <- as.integer(nhlim-1)
    rank <- as.integer(0)
    gam <- as.double( numeric(mlim) )
    perm <- integer(mlim)
    Q <- matrix( as.double(numeric(n*mlim)), nrow=n, ncol=mlim )
    R <- matrix( as.double(numeric(mlim*mlim)), nrow=mlim, ncol=mlim )

    xhist <- matrix( numeric(n*nhlim), nrow=n, ncol=nhlim )
    rhist <- matrix( numeric(n*nhlim), nrow=n, ncol=nhlim )

    fargs <- list(...)

    fcnt <- iter <- itersafe <- maxrank <- 0

    nexti <- 1

    lastfv <- array( -Inf, c( M ) )

    creatSY <- maxsizeXY <- FALSE

    res <- try( do.call( evalr, append( list( x ), fargs ) ) )
    fcnt <- fcnt + 1
    if ( inherits( res, "try-error" ) )
        stop( "Failure in initial functional evaluation. \n" )
    else if ( !is.numeric( res ) || !is.vector( res ) )
        stop( "Function must return a vector numeric value." )
    else if ( any( is.nan( res ), is.infinite( res ), is.na( res ) ) )
        stop( "Failure in initial functional evaluation. \n" )

    f <- sum( res * res )

    nh <- 1

    xhist[1:n,1] <- x
    rhist[1:n,1] <- res

    eta <- min( 0.5 * f, sqrt( f ) )

    lamsg <- 1.0

    ## MAIN LOOP
    while ( TRUE ) {
        if ( iprint >= 0 ) {
            if ( iprint > 0 ) message( "\n" )
            message( "Iter: ", iter, " f = ", f, "\n" )
        }

        if ( tabline ) {
            sqrtf <- sqrt(f)
            sqrtf <- formatC(sqrtf, format = "E", digits = 1, mode="double")
            textfile <- file.path ( "tabline-interrupted-tmp.txt" )
            printer <- file( textfile, "w" )
            write(c(
                formatC(n, format="d", digits=6),
		formatC(9, format="d", digits=2),
                formatC(sqrtf, format="E", digits=1, width=7, mode="double"),
                formatC(iter, format="d", digits=9),
                formatC(fcnt, format="d", digits=9),
                formatC(180, format="f", digits=6, width=12, mode="double"),
                formatC(itersafe/iter, format="f", digits=6, width=9, mode="double" )
            ),
            textfile, sep=" ", append=TRUE, ncolumns=7 )
            close( printer )
            file.copy("tabline-interrupted-tmp.txt", "tabline-interrupted.txt", overwrite=TRUE)
        }

        if ( f <= epsf ^ 2 ) {
            if ( iprint >= 0 )
                message ("success!\n")
            istop <- 0
            break
        }

        else if ( iter >= maxit ) {
            if ( iprint >= 0 )
                message("maximum number of iterations reached!\n")
            istop <- 1
            break
        }

        lastfv[ iter%%M + 1 ] <- f
        fmax <- max( lastfv )

        iter <- iter + 1

        backcnt <- 0

        ## Backtracking
        while ( TRUE ) {
            if ( backcnt == 0 )
                lambda1 <- lamsg
            else {
                ltmp <- lambda1 ^ 2 * f / ( ftrial1 + ( 2.0 * lambda1 - 1.0 ) * f )
                if ( is.nan( ltmp ) )
                    lambda1 <- sigma1 * lambda1
                else
                    lambda1 <- max( sigma1 * lambda1, min( ltmp, sigma2 * lambda1 ) )
            }

            xtrial <- x - lambda1 * res

            restrial <- try( do.call( evalr, append( list( xtrial ), fargs ) ) )
            fcnt <- fcnt + 1

            ftrial1 <- sum( restrial * restrial )

            if ( iprint >= 1 )
                message( "lambda = ", lambda1, "  ftrial+ = ", ftrial1, "\n" )

            if ( ! is.nan( ftrial1 ) ) {
                if ( ftrial1 <= fmax + eta - 2.0 * gamma * lambda1 ^ 2 * f ) {
                    ftrial <- ftrial1
                    break
                }
            }

            if ( backcnt == 0 )
                lambda2 <- lamsg
            else {
                ltmp <- lambda2 ^ 2 * f / ( ftrial2 + ( 2.0 * lambda2 - 1.0 ) * f )
                if ( is.nan( ltmp ) )
                    lambda2 <- sigma1 * lambda2
                else
                    lambda2 <- max( sigma1 * lambda2, min( ltmp, sigma2 * lambda2 ) )
            }

            xtrial <- x + lambda2 * res

            restrial <- try( do.call( evalr, append( list( xtrial ), fargs ) ) )
            fcnt <- fcnt + 1

            ftrial2 <- sum( restrial * restrial )

            if ( iprint >= 1 )
                message( "lambda = ", lambda2, "ftrial- = ", ftrial2, "\n" )

            if ( ! is.nan( ftrial2 ) ) {
                if ( ftrial2 <= fmax + eta - 2.0 * gamma * lambda2 ^ 2 * f ) {
                    ftrial <- ftrial2
                    break
                }
            }

            backcnt <- backcnt + 1
        }

        ## ACCELERATION
        if ( acc && iter%%ndiis == 0 ) {
            if ( iprint >=2 )
                message("Including xtrial in history\n")

            if (nh + 1 <= min( nhlim, n + 1 ) ) {
                if ( iprint >= 2 )
                    message("There is space to include xtrial without removing any other point ",
                        "(nh = ",nh," min( nhlim, n + 1 ) = ",min( nhlim, n + 1 ),").\n")

                nh <- nh + 1
                jtbr <- as.integer( 0 )
            }
            else {
                if ( iprint >= 2 )
                    message("History is full. Oldest point is being removed ",
                        "(nh = ",nh," min( nhlim, n + 1 ) = ",min( nhlim, n + 1 ),").\n")

                xhist[1:n,1:(nh-1)] <- xhist[1:n,2:nh]
                rhist[1:n,1:(nh-1)] <- rhist[1:n,2:nh]
                jtbr <- as.integer( 1 )
            }

            xhist[1:n,nh] <- xtrial
            rhist[1:n,nh] <- restrial

            col <- matrix( as.double(rhist[1:n,nh] - rhist[1:n,nh-1]), nrow=n, ncol=1 )

            extra <- FALSE

            if ( nh == 2 ) {
                ncol <- as.integer(1)

                .Call( c_bmqr, n, ncol, mlim, n, col, gam, perm, rank )
                .Call( c_extractQ, n, ncol, mlim, n, col, gam, rank, n, Q )
                .Call( c_extractR, ncol, mlim, n, col, rank, mlim, R )

                maxrank <- max( maxrank, rank )
            }
            else {
                .Call( c_bmqrupdate, n, ncol, rank, mlim, jtbr, as.integer(1), n, Q, mlim, R, perm, col )

                maxrank <- max( maxrank, rank )

                if ( iprint >= 2 )
                    message("Number of columns = ",ncol," rank = ",rank," (largest achieved rank = ",maxrank,")\n")

                if ( ncol != nh - 1)
                    stop("1) ERROR in subroutine dfsane: ncol is not equal to nh - 1! ncol = ",ncol," nh = ",nh)

                if ( rank < maxrank ) {
                    if ( iprint >= 2 )
                        message('Current rank is smaller that largest achieved rank. Thus, an extra point will be considered.\n')

                    extra <- TRUE

                    xtmp <- x
                    xtmp[nexti] <- x[nexti] + hsmall * max( 1.0, abs( x[nexti] ) )

                    nexti <- nexti%%n + 1

                    restmp <- try( do.call( evalr, append( list( xtmp ), fargs ) ) )
                    fcnt <- fcnt + 1

                    if ( iprint >= 2 )
                        message('Including extra point in history\n')

                    if (nh + 1 <= min( nhlim, n + 1 ) ) {
                        if ( iprint >= 2 )
                            message("There is space to include the extra point without removing any other point ",
                                "(nh = ",nh," min( nhlim, n + 1 ) = ",min( nhlim, n + 1 ),").\n")

                        nh <- nh + 1
                        jtbr <- as.integer( 0 )
                    }
                    else {
                        if ( iprint >= 2 )
                            message("History is full. Oldest point is being removed ",
                                "(nh = ",nh," min( nhlim, n + 1 ) = ",min( nhlim, n + 1 ),").\n")

                        xhist[1:n,1:(nh-1)] <- xhist[1:n,2:nh]
                        rhist[1:n,1:(nh-1)] <- rhist[1:n,2:nh]
                        jtbr <- as.integer( 1 )
                    }

                    xhist[1:n,nh] <- xtmp
                    rhist[1:n,nh] <- restmp

                    col <- matrix( as.double( rhist[1:n,nh] - rhist[1:n,nh-1] ), nrow=n, ncol=1 )

                    .Call( c_bmqrupdate, n, ncol, rank, mlim, jtbr, as.integer(1), n, Q, mlim, R, perm, col )

                    maxrank <- max( maxrank, rank )

                    if ( iprint >= 2 )
                        message("Number of columns = ",ncol," rank = ",rank," (largest achieved rank = ",maxrank,")\n")

                    if ( ncol != nh - 1)
                        stop("2) ERROR in subroutine dfsane: ncol is not equal to nh - 1! ncol = ",ncol," nh = ",nh)
                }
            }

            if ( rank > 0 ) {
                omega <- as.double( numeric(ncol) )
                .Call( c_bmqrsolve, n, ncol, rank, n, Q[1:n,1:rank], rank, R[1:rank,1:ncol], perm, restrial, omega )

                xacc <- xtrial - as.matrix( xhist[1:n,2:nh] - xhist[1:n,1:(nh-1)], nrow=n, ncol=ncol ) %*% omega[1:ncol]

                if ( extra ) {
                    if ( iprint >= 2 )
                        message("The extra point is being removed.\n")

                    nh <- nh - 1

                    jtbr <- as.integer( ncol )

                    .Call( c_bmqrupdate, n, ncol, rank, mlim, jtbr, as.integer(0), n, Q, mlim, R, perm, col )

                    maxrank <- max( maxrank, rank )

                    if ( iprint >= 2 )
                        message("Number of columns = ",ncol," rank = ",rank," (largest achieved rank = ",maxrank,")\n")

                    if ( ncol != nh - 1)
                        stop("3) ERROR in subroutine dfsane: ncol is not equal to nh - 1! ncol = ",ncol," nh = ",nh)
                }

                if ( sqrt( sum( xacc * xacc ) ) <= 10.0 * max( 1.0, sqrt( sum( x * x ) ) ) ) {

                    resacc <- try( do.call( evalr, append( list( xacc ), fargs ) ) )
                    fcnt <- fcnt + 1

                    facc <- sum( resacc * resacc )

                    if ( iprint >= 1 )
                        message("    facc = ", facc,"\n")

                    if ( facc < ftrial && sqrt( sum( ( xacc - x ) * ( xacc - x ) ) ) >
                         macheps * max( 1.0, sqrt( sum( x * x ) ) ) ) {

                        if ( iprint >= 2 )
                            message("Accelerated point is being accepted; and history will be updated substituting xtrial by xacc.\n")

                        xhist[1:n,nh] <- xacc
                        rhist[1:n,nh] <- resacc

                        jtbr <- as.integer( ncol )

                        col <- matrix( as.double( rhist[1:n,nh] - rhist[1:n,nh-1] ), nrow=n, ncol=1 )

                        .Call( c_bmqrupdate, n, ncol, rank, mlim, jtbr, as.integer(1), n, Q, mlim, R, perm, col )

                        maxrank <- max( maxrank, rank )

                        if ( iprint >= 2 )
                            message("Number of columns = ",ncol," rank = ",rank," (largest achieved rank = ",maxrank,")\n")

                        if ( ncol != nh - 1)
                            stop("4) ERROR in subroutine dfsane: ncol is not equal to nh - 1! ncol = ",ncol," nh = ",nh)

                        xtrial   <- xacc
                        restrial <- resacc
                        ftrial   <- facc
                    }
                }
            }

            if ( rank == 0 ) {
                if ( iprint >= 0 ) {
                    message("Acceleration matrix has null rank; thus, it will be rebuild from scratch.\n")
                    message("It will be populated with neighbours of the current point taking a small step into coordinate directions.\n")
                }

                nh <- 0

                while ( nh + 1 < min( nhlim, n + 1 ) ) {
                    xtmp <- xtrial
                    xtmp[nexti] <- xtrial[nexti] + hlarge * max( 1.0, abs( xtrial[nexti] ) )

                    nexti <- nexti%%n + 1

                    restmp <- try( do.call( evalr, append( list( xtmp ), fargs ) ) )
                    fcnt <- fcnt + 1

                    nh <- nh + 1

                    xhist[1:n,nh] <- xtmp
                    rhist[1:n,nh] <- restmp


                    if ( nh >= 2 ) {
                        col <- matrix( as.double( rhist[1:n,nh] - rhist[1:n,nh-1] ), nrow=n, ncol=1 )

                        if ( nh == 2 ) {
                            ncol <- as.integer(1)

                            .Call( c_bmqr, n, ncol, mlim, n, col, gam, perm, rank )
                            .Call( c_extractQ, n, ncol, mlim, n, col, gam, rank, n, Q )
                            .Call( c_extractR, ncol, mlim, n, col, rank, mlim, R )

                            maxrank <- max( maxrank, rank )
                        }
                        else {
                            jtbr <- as.integer( 0 )

                            .Call( c_bmqrupdate, n, ncol, rank, mlim, jtbr, as.integer(1), n, Q, mlim, R, perm, col )

                            maxrank <- max( maxrank, rank )

                            if ( iprint >= 2 )
                                message("Number of columns = ",ncol," rank = ",rank," (largest achieved rank = ",maxrank,")\n")

                            if ( ncol != nh - 1)
                                stop("5) ERROR in subroutine dfsane: ncol is not equal to nh - 1! ncol = ",ncol," nh = ",nh)
                        }
                    }
                }

                nh <- nh + 1

                xhist[1:n,nh] <- xtrial
                rhist[1:n,nh] <- restrial

                jtbr <- as.integer( 0 )

                col <- matrix( as.double( rhist[1:n,nh] - rhist[1:n,nh-1] ), nrow=n, ncol=1 )

                .Call( c_bmqrupdate, n, ncol, rank, mlim, jtbr, as.integer(1), n, Q, mlim, R, perm, col )

                maxrank <- max( maxrank, rank )

                omega <- as.double( numeric(ncol) )
                .Call( c_bmqrsolve, n, ncol, rank, n, Q[1:n,1:rank], rank, R[1:rank,1:ncol], perm, restrial, omega[1:ncol] )

                xacc <- xtrial - as.matrix( xhist[1:n,2:nh] - xhist[1:n,1:nh-1], nrow=n, ncol=ncol ) %*% omega[1:ncol]

                resacc<- try( do.call( evalr, append( list( xacc ), fargs ) ) )
                fcnt <- fcnt + 1

                facc <- sum( resacc * resacc )

                if ( iprint >= 1 )
                    message("    facc = ",facc,"\n")

                if ( facc < ftrial && sqrt( sum( ( xacc - x ) * ( xacc - x ) ) ) >
                     macheps * max( 1.0, sqrt( sum( x * x ) ) ) ) {

                    xhist[1:n,nh] <- xacc[1:n]
                    rhist[1:n,nh] <- resacc[1:n]

                    jtbr <- as.integer( ncol )

                    col <- matrix( as.double( rhist[1:n,nh] - rhist[1:n,nh-1] ), nrow=n, ncol=1 )

                    .Call( c_bmqrupdate, n, ncol, rank, mlim, jtbr, as.integer(1), n, Q, mlim, R, perm, col )

                    maxrank <- max( maxrank, rank )

                    xtrial   <- xacc
                    restrial <- resacc
                    ftrial   <- facc

                }
            }
        }

        s <- xtrial - x
        lamsg <- sum( s * s ) / sum( s * ( restrial - res ) )
        if ( !( lamsgmin <= abs( lamsg ) && abs( lamsg ) <= 1.0 ) )
            lamsg = max( lamsgmin, min( sqrt( sum( xtrial * xtrial ) ) / sqrt( ftrial ), lamsgmax ) )

        eta <- 0.5 * eta

        x   <- xtrial
        res <- restrial
        f   <- ftrial
    }

    return ( list( x=x, res=res, normF=f, iter=iter, fcnt=fcnt, istop=istop ) )
}
