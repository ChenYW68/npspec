
#' @export
maternFun <- function(distmat, covparms){

    # covariance parameters are...
    # (variance (including nugget), range, smoothness, signal to signal + noise ratio)
    # range of values are...
    # (0,Inf), (0,Inf) , (0,Inf), (0,1)

    n1 <- dim(distmat)[1]
    n2 <- dim(distmat)[2]
    covmat      <-  matrix(0,n1,n2)
    scaledist   <-  sqrt(2*covparms[3])*distmat/covparms[2]

    # to avoid bad spots for smoothness parameter
    covparms[3] <- min(max(covparms[3],1e-12),20)
    normcon     <-  2^(covparms[3]-1)*gamma(covparms[3])

    # special cases for nu = 1/2, 1, and 3/2
    #if( covparms[3] == 1/2 ){ maternpart <- exp(-scaledist) }
    #else if( covparms[3] == 3/2 ){ maternpart <- (1 + scaledist)*exp(-scaledist) }
    #else if( covparms[3] == 1   ){ maternpart <- MaternFun(distmat,c(1,covparms[2:3])) }
    #else {
        besspart <- besselK(scaledist,covparms[3])
        maternpart  <-  1/normcon*(scaledist)^covparms[3]*besspart
    #}
    covmat[distmat != 0] <- covparms[1]*maternpart[distmat!=0]
    covmat[distmat == 0] <- covparms[1]*(1 + covparms[4])
    return(covmat)

}




#' @export
get_exponential2D_spec <- function(nvec,parms){

    # evaluate the 2D exponential spectral density
    # with parameters parms ( parms[1]*exp(-d/parms[2]) )
    # at Fourier frequencies of resolution nvec

    # calculates wrapped covariances on a grid of size
    # nvec until tol is reached, then does fft of wrapped covariances

    tol <- 1e-8

    # to make sure we don't do too much computation
    kmax <- round( sqrt(6000^2/prod(nvec))/2 )

    # distance matrix for (j1,j2) = (0,0)
    # this is used to create translated distances
    dist1 <- matrix( c(0:(nvec[1]-1)), nvec[1], nvec[2] )
    dist2 <- matrix( c(0:(nvec[2]-1)), nvec[1], nvec[2], byrow = TRUE )


    wrapcov <- array(0,nvec)
    for(k in 1:kmax){

        print(k)
        wrapcov_k <- array(0,nvec)
        for( j1 in -k:(k-1) ){
            for( j2 in -k:(k-1) ){
                dist1_trans <- dist1 + j1*nvec[1]
                dist2_trans <- dist2 + j2*nvec[2]

                distmat <- sqrt( dist1_trans^2 + dist2_trans^2 )
                cov_trans <- parms[1]*exp(-distmat/parms[2])
                wrapcov_k <- wrapcov_k + cov_trans
            }
        }

        # compare wrapcov_k to wrapcov_k-1
        if( max( abs( wrapcov_k - wrapcov ) ) < tol ){
            wrapcov <- wrapcov_k
            break
        } else {
            wrapcov <- wrapcov_k
        }

        if( k == kmax ){
            stop("Increase kmax or use different method for calculating spectral density")
        }
    }
    specden <- Re(fft(wrapcov))
    return(specden)
}



#' @export
get_matern2D_spec <- function(nvec,parms){

    # evaluate the 2D matern spectral density
    # at Fourier frequencies of resolution nvec

    # calculates wrapped covariances on a grid of size
    # nvec until tol is reached, then does fft of wrapped covariances

    tol <- 1e-8

    # to make sure we don't do too much computation
    kmax <- round( sqrt(6000^2/prod(nvec))/2 )

    # distance matrix for (j1,j2) = (0,0)
    # this is used to create translated distances
    dist1 <- matrix( c(0:(nvec[1]-1)), nvec[1], nvec[2] )
    dist2 <- matrix( c(0:(nvec[2]-1)), nvec[1], nvec[2], byrow = TRUE )


    wrapcov <- array(0,nvec)
    for(k in 1:kmax){

        print(k)
        wrapcov_k <- array(0,nvec)
        for( j1 in -k:(k-1) ){
            for( j2 in -k:(k-1) ){
                dist1_trans <- dist1 + j1*nvec[1]
                dist2_trans <- dist2 + j2*nvec[2]

                distmat <- sqrt( dist1_trans^2 + dist2_trans^2 )
                cov_trans <- maternFun( distmat, parms )
                wrapcov_k <- wrapcov_k + cov_trans
            }
        }

        # compare wrapcov_k to wrapcov_k-1
        if( max( abs( wrapcov_k - wrapcov ) ) < tol ){
            wrapcov <- wrapcov_k
            break
        } else {
            wrapcov <- wrapcov_k
        }

        if( k == kmax ){
            stop("Increase kmax or use different method for calculating spectral density")
        }
    }
    specden <- Re(fft(wrapcov))
    return(specden)
}




#' @export
uncond_sim_exponential2D <- function(nvec,parms){

    # uses circulant embedding to do exact simulations
    # from exponential covariance model
    nvec0 <- nvec
    posdef <- FALSE
    maxcount <- 5
    counter <- 0
    while( !posdef & counter < maxcount ){
        dist1 <- matrix( c(0:(nvec0[1]-1)), nvec0[1], nvec0[2] )
        dist2 <- matrix( c(0:(nvec0[2]-1)), nvec0[1], nvec0[2], byrow = TRUE )
        distarray <- sqrt( dist1^2 + dist2^2 )
        covarray <- parms[1]*exp(-distarray/parms[2])
        covarray_wide <- cbind( covarray, covarray[1:nvec0[1],nvec0[2]:2] )
        covarray_embed <- rbind( covarray_wide, covarray_wide[ nvec0[1]:2, 1:(2*nvec0[2]-1) ] )
        spec_embed <- Re( fft(covarray_embed) )
        if( min(spec_embed) <= 0 ){
            counter <- counter+1
            nvec0 <- 2*nvec0
        } else {
            posdef <- TRUE
        }
    }
    if( counter == maxcount ){
        print("Did not converge, increase maxcount")
        return(NULL)
    } else {
        y <- uncond_sim( spec_embed )[1:nvec[1],1:nvec[2]]
        return(y)
    }
}





#' @export
uncond_sim_matern2D <- function(nvec,parms){

    # uses circulant embedding to do exact simulations
    # from exponential covariance model
    nvec0 <- nvec
    posdef <- FALSE
    maxcount <- 5
    counter <- 0
    while( !posdef & counter < maxcount ){
        dist1 <- matrix( c(0:(nvec0[1]-1)), nvec0[1], nvec0[2] )
        dist2 <- matrix( c(0:(nvec0[2]-1)), nvec0[1], nvec0[2], byrow = TRUE )
        distarray <- sqrt( dist1^2 + dist2^2 )

        covarray <- maternFun( distarray, parms )

        covarray_wide <- cbind( covarray, covarray[1:nvec0[1],nvec0[2]:2] )
        covarray_embed <- rbind( covarray_wide, covarray_wide[ nvec0[1]:2, 1:(2*nvec0[2]-1) ] )
        spec_embed <- Re( fft(covarray_embed) )
        if( min(spec_embed) <= 0 ){
            counter <- counter+1
            nvec0 <- 2*nvec0
        } else {
            posdef <- TRUE
        }
    }
    if( counter == maxcount ){
        print("Did not converge, increase maxcount")
        return(NULL)
    } else {
        y <- uncond_sim( spec_embed )[1:nvec[1],1:nvec[2]]
        return(y)
    }
}
