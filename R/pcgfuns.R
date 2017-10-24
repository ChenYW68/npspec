
#' @export
pcg_spec <- function( y, spec, obs, precondmethod = "fft", maxit = 500, tol = 1e-6, silent = FALSE, NNarray = NULL ){

    nvec <- dim(spec)


    fmult <- function(x,spec,obs){
        n <- prod(dim(x))
        x[!obs] = 0
        x1 <- fft(x)
        x1 <- spec*x1
        x1 <- fft( x1, inverse = TRUE )
        x1 <- 1/n*Re( x1 )
        #x1 <- Re( 1/n*fft( (spec)*fft(x), inverse = TRUE ) )
    }

    if( precondmethod == "fft"){
        #cat("Using FFT \n")
        pmult <- function(x,spec,obs){
            n <- prod(dim(x))
            x[!obs] = 0
            x1 <- Re( 1/n*fft( fft(x)/(spec), inverse = TRUE ) )
        }
    } else if(precondmethod == "Vecchia"){
        #cat("Using Vecchia \n")
        if(is.null(NNarray)) stop("Must supply nearest neighbor array for Vecchia preconditioner")
        nvec <- dim(y)
        n <- prod(nvec)
        locs <- as.matrix(expand.grid(1:nvec[1],1:nvec[2]))
        locs <- locs[obs,]
        covarray <- 1/n*Re( fft(spec, inverse = TRUE ) )

        LinvEntries <- getLinvEntries(covarray,locs,NNarray)

        pmult <- function(x,spec,obs){
            yvec <- x[obs]

            #xvec <- vecchiaPrecond(covarray,yvec,locs,NNarray)

            xvec <- vecchiaPrecondFromEntries(LinvEntries,yvec,NNarray)

            x1 <- matrix(0,nvec[1],nvec[2])
            x1[obs] <- xvec
            return(x1)
        }
    } else { stop("Invalid Preconditioning Method precondmethod. Use \"fft\" or \"Vecchia\" ") }

    b <- y
    b[!obs] <- 0
    x <- matrix(0,nvec[1],nvec[2])
    r <- fmult(x,spec,obs)-b
    z <- pmult(r,spec,obs)
    z[!obs] <- 0
    p <- -z

    for(k in 1:maxit){
        p[!obs] <- 0
        Ap <- fmult(p,spec,obs)
        a <- sum( r[obs]*z[obs] )/sum( p[obs]*Ap[obs] )
        x <- x + a*p
        r1 <- r + a*Ap
        if( sqrt(sum( r1[obs]^2 )/sum(obs)) < tol ){
            if(silent==FALSE){
                cat( paste0("Converged at Iteration ",k,", RMSE = ", round(sqrt( sum( r1[obs]^2 )/sum(obs)),8),"\n" ))
            }
            break
        }
        z1 <- pmult(r1,spec,obs)
        beta <- sum( z1[obs]*r1[obs] )/sum( z[obs]*r[obs] )
        p <- -z1 + beta*p
        r <- r1
        z <- z1
        if(k==maxit){
            cat(paste0("Did not converge after ",maxit," iterations, RMSE = ", round(sqrt(sum(r1[obs]^2)/sum(obs)),8)))
            cat("\nConsider increasing maxit or trying a different preconditioner\n")
        }
    }

    x[!obs] <- 0
    yhat <- fmult(x,spec,obs)
    return(list(x=x,yhat=yhat))

}
