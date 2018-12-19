

#' @export
uncond_sim <- function(spec){
    nvec <- dim(spec)
    n <- prod(nvec)
    z <- matrix(rnorm(n),nvec[1],nvec[2])
    y <- Re(1/sqrt(n)*fft( sqrt(spec)*1/sqrt(n)*fft(z), inverse = TRUE ))
    return(y)
}





#' @export
condexp_spec <- function(y,spec,obs,silent = FALSE, maxit = 500, precondmethod = "fft", ...){

    pcg_result <- pcg_spec(y,spec,obs, silent = silent, maxit = maxit, precondmethod = precondmethod, ...)
    return(pcg_result$yhat)

}


#' @export
condsim_spec <- function(y,spec,obs, silent = FALSE, maxit = 500, precondmethod = "fft", ...){

    eps <- uncond_sim(spec)
    z1 <- condexp_spec( y + eps, spec = spec, obs = obs, silent = silent, maxit = maxit, precondmethod = precondmethod, ...)
    z2 <- y
    z2[!obs] = z1[!obs] - eps[!obs]
    return(z2)

}






