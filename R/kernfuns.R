
#' @export
sqexp_kern <- function( r, nvec ){

    # r is a range parameter
    # larger r gives stronger smoothing

    v1 <- as.matrix( c( seq( 0, ceiling((nvec[1]-1)/2) ),
                        seq( floor((nvec[1]-1)/2), 1, by = -1 ) ) )/nvec[1]
    v2 <- as.matrix( c( seq( 0, ceiling((nvec[2]-1)/2) ),
                        seq( floor((nvec[2]-1)/2), 1, by = -1 ) ) )/nvec[2]

    v1arr <-    v1[,rep(1,nvec[2])]
    v2arr <- t( v2[,rep(1,nvec[1])] )

    kern <- exp( -(v1arr^2 + v2arr^2)/r^2 )
    kern <- kern/sum(kern)
    return(kern)

}

#' @export
smooth_pgram <- function(pgram,kern, smoothlog = FALSE){

    n <- prod(dim(pgram))
    if( !smoothlog ){
        smpgram <- Re(1/n*fft( fft(pgram)*fft(kern), inverse = TRUE ))
    } else if( smoothlog ){
        sumpgram <- sum(pgram)
        logpgram <- log(pgram)
        logsmpgram <- Re(1/n*fft( fft(logpgram)*fft(kern), inverse = TRUE ))
        smpgram <- exp(logsmpgram)
        smpgram <- sumpgram*smpgram/sum(smpgram)
    }
    return(smpgram)

}





