cutest_init <- function(probname, interfaceDir=".", sifdecode=TRUE, compile=TRUE) {
    ## Get environment variables
    env <- Sys.getenv( c( "CC", "CFLAGS", "CUTEST", "FC", "FFLAGS", "MYARCH", "SIFDECODE" ), unset=NA )

    if ( is.na( env["CC"] ) ) env["CC"] <- "gcc"
    if ( is.na( env["FC"] ) ) env["FC"] <- "gfortran"

    if ( is.na( env["CFLAGS"] ) ) env["CFLAGS"] <- ""
    if ( is.na( env["FFLAGS"] ) ) env["FFLAGS"] <- ""

    ## Run sifdecoder from CUTEst to generate source code to problem, i.e.:
    ## ELFUN.f, EXTER.f, GROUP.f, RANGE.f, AUTOMAT.d and OUTSDIF.d
    if ( sifdecode ) {
        if ( is.na( env["CUTEST"] ) || is.na( env["MYARCH"] ) || is.na( env["SIFDECODE"] ) )
            stop( "The environment variables CUTEST, MYARCH, and SIFDECODE related to CUTEst are not properly set." )

        cmd <- paste0(env["SIFDECODE"], "/bin/sifdecoder -A ",
                      env["MYARCH"], " ", probname, ".SIF")
        system(cmd)
    }

    ## Compile code with $FC and $CC
    if ( compile ) {
        cmd <- paste0(env["FC"], " ", env["FFLAGS"],
                      " -c ELFUN.f EXTER.f GROUP.f RANGE.f")
        system(cmd)

        cmd <- paste0(env["CC"], " -fPIC -I ",
                      R.home(component="include"), " -I ",
                      env["CUTEST"], "/include/ ", env["CFLAGS"],
                      " -c ", interfaceDir, "/cutest.c")
        system(cmd)

        cmd <- paste0(env["CC"], " -shared -fPIC ", env["CFLAGS"],
                      " -L ", env["CUTEST"], "/objects/",
                      env["MYARCH"], "/double -L ",
                      R.home( component="lib" ),
                      " ELFUN.o EXTER.o GROUP.o RANGE.o cutest.o ",
                      "-o R-problem-forcutest.so -lcutest -lR -lgfortran\n")
        system(cmd)

        ## Remove compilation files
        file.remove("ELFUN.o")
        file.remove("EXTER.o")
        file.remove("GROUP.o")
        file.remove("RANGE.o")
        file.remove("ELFUN.f")
        file.remove("EXTER.f")
        file.remove("GROUP.f")
        file.remove("RANGE.f")
        file.remove("cutest.o")
        file.remove("AUTOMAT.d")
    }

    if ( !file.exists("R-problem-forcutest.so") )
        stop("Something went wrong in generation of shared library problem-forcutest.so.")

    dyn.load("R-problem-forcutest.so")

    n <- .Call( "c_cutest_getn" )
    assign( "Global.n", as.integer(n), envir=.GlobalEnv )
    
    assign( "Global.x0", as.vector( as.double( rep( 0, Global.n ) ) ), envir=.GlobalEnv )
    assign( "Global.crhs", as.vector( as.double( rep( 0, Global.n ) ) ), envir=.GlobalEnv )

    ret <- .Call( "c_cutest_init", Global.x0, Global.crhs )
}

cutest_getn <- function() {
    if ( exists( "Global.n", envir=.GlobalEnv ) ) Global.n
    else NA
}

cutest_getx0 <- function() {
    if ( exists( "Global.x0", envir=.GlobalEnv ) ) Global.x0
    else NA
}

cutest_evalr <- function(x) {
    n <- length( x )
    f <- as.double( rep( 0, n ) )
    .Call( "c_cutest_evalr", n, x, Global.crhs )
}

cutest_end <- function() {
    .Call( "c_cutest_end" )
    if ( is.loaded( "c_cutest_init" ) ) dyn.unload("R-problem-forcutest.so")
    if ( file.exists( "R-problem-forcutest.so" ) ) file.remove("R-problem-forcutest.so")
    if ( exists( "Global.n", envir=.GlobalEnv ) ) remove( Global.n, envir=.GlobalEnv )
    if ( exists( "Global.crhs", envir=.GlobalEnv ) ) remove( Global.crhs, envir=.GlobalEnv )
    if ( exists( "Global.x0", envir=.GlobalEnv ) ) remove( Global.x0, envir=.GlobalEnv )
    if ( file.exists( "OUTSDIF.d" ) ) file.remove("OUTSDIF.d")
}
