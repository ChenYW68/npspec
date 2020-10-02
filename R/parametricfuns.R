# specfun <- function (parm, nvec){
#     n <- prod(nvec)
#     filter0 <- array(0, nvec)
#     filter0[1, 1] = 1
#     filter0[1, 2] = -parm
#     filter0[2, 1] = -parm
#     filter0[nvec[1], 1] = -parm
#     filter0[1, nvec[2]] = -parm
#     fftfilter <- fft(filter0)
#     spec <- 1/(abs(fftfilter)^2)
#     spec <- n * spec/sum(spec)
#     return(spec)
# }

whittle_lik <- function(pgram,parms,specfun){

    nvec <- dim(pgram)
    spec <- specfun(parms,nvec)   # i.e. spec_AR1
    loglik <- -prod(nvec)/2*log(2*pi)-sum(log(spec)) - sum(pgram/spec)
    return(loglik)
}

logit <- function(x){ log(x/(1-x)) }
expit <- function(x){ exp(x)/(1+exp(x)) }




qmatern <- function( nvec, covparms ){

    # this is for 2 dimensions, can we write
    # it generally for nd?

    # covparms contains sigsq,a1,a2,theta,sm,nug

    # compute sine functions
    n1 <- nvec[1]
    n2 <- nvec[2]
    n <- n1*n2
    w1 <- (0:(n1-1))*2*pi/n1
    w2 <- (0:(n2-1))*2*pi/n2
    s1 <- matrix(sin(w1/2))
    s2 <- matrix(sin(w2/2))
    s1mat <-   s1[,rep(1,n2)]
    s2mat <- t(s2[,rep(1,n1)])

    #
    a1 <- covparms[2]
    a2 <- covparms[3]
    th <- covparms[4]
    sm <- covparms[5]
    nug <- covparms[6]
    cth <- cos(th)
    sth <- sin(th)
    con <- 1/( a1^2 + a2^2 )
    con <- 1
    d   <- con +  a1^2*con*s1mat^2 + a2^2*con*s2mat^2
    s <- d^(-sm-1)
    s <- n*covparms[1]*s/sum(s) + 0*nug
    s
}

#' @export
spec_AR1 <- function(parm,nvec){

    n <- prod(nvec)
    filter0 <- array(0,nvec)
    filter0[1,1] = 1
    filter0[1,2] = -parm
    filter0[2,1] = -parm
    filter0[nvec[1],1] = -parm
    filter0[1,nvec[2]] = -parm

    fftfilter <- fft(filter0)
    spec <- 1/(abs(fftfilter)^2)
    spec <- n*spec/sum(spec)
    return(spec)

}


#' @export
spec_qmatern <- function( covparms, nvec ){

    # compute sine functions
    n1 <- nvec[1]
    n2 <- nvec[2]
    n <- n1*n2
    w1 <- (0:(n1-1))*2*pi/n1
    w2 <- (0:(n2-1))*2*pi/n2
    s1 <- matrix(sin(w1/2))
    s2 <- matrix(sin(w2/2))
    s1mat <-   s1[,rep(1,n2)]
    s2mat <- t(s2[,rep(1,n1)])

    #
    a1 <- covparms[1]
    sm <- covparms[2]
    con <- 1
    d   <- con +  a1^2*s1mat^2 + a1^2*s2mat^2
    s <- d^(-sm-1)
    s <- n*s/sum(s)
    return(s)
}
