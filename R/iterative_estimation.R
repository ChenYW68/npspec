
#' @import aldodevel
#' @import FNN
#' @export
iterate_spec_semipar <- function(y, observed, embed_fac = 1.2, n_iter = 100, par_spec_fun = spec_AR1,
                                 kern_parm = 10, n_avg = 30, precondmethod = "fft", m = 10,
                                 silent = TRUE, ncondsim = 0){

    # pass NNarray through "..."

    if(identical(par_spec_fun,FALSE)){
        do_parametric_filter <- FALSE
    } else {
        do_parametric_filter <- TRUE
    }

    # error for bad embedding factor value
    if( embed_fac < 1 ){
        stop("Embedding factor embed_fac cannot be less than 1")
    }

    # create embedded lattice
    nvec_obs <- dim(y)
    nvec <- round( nvec_obs*embed_fac )
    y_embed <- array(NA,nvec)
    y_embed[1:nvec_obs[1],1:nvec_obs[2]] <- y
    observed_embed <- array(FALSE,nvec)
    observed_embed[1:nvec_obs[1],1:nvec_obs[2]] <- observed

    if( precondmethod == "Vecchia"){
        locsfull <- expand.grid( 1:nvec[1], 1:nvec[2] )
        locs <- locsfull[observed_embed,]
        NNarray <- aldodevel::findOrderedNN_kdtree(locs,m)
        NNarray[m+1,] <- (m+1):1 # need to make sure first part goes in the right order
    } else {
        NNarray <- NULL
    }

    # get grid size and define kerenel
    n <- prod(nvec)
    kern <- sqexp_kern(kern_parm, nvec)

    # impute with mean
    y0 <- y_embed
    y0[!observed_embed] <- mean(y, na.rm = TRUE)

    # define likelihood function (written for AR1 filter)
    likfun <- function(x){
        expitx <- expit(x)/4
        return(-whittle_lik(pgram,expitx,par_spec_fun))
    }

    # set an initial value for optimization
    parm <- 1/8
    logitpar <- logit(parm*4)

    # loop over number of iterations
    specs <- array(NA,c(nvec,n_iter))

    for(k in 1:n_iter){

        # periodogram
        pgram <- 1/n*abs( fft(y0) )^2

        # do the optimization
        if(do_parametric_filter){
            res <- optim(logitpar,likfun,method="Brent",lower=-6,upper=6)
            logitpar <- res$par         # update logit parameter
            parm <- expit(logitpar)/4   # update parameter

            # parametric approximation
            param_spec <- par_spec_fun(parm,nvec)
        } else {
            param_spec <- array(1,nvec)
        }

        # smooth ratio of periodogram to parametric approximation
        sm_pgram <- smooth_pgram(pgram/param_spec,kern)
        # estimate is product of smoothed ratio and parametric spec
        specs[,,k] <- sm_pgram*param_spec

        # do a conditional simulation with the new spectrum
        y0 <- condsim_spec(y = y0,spec = specs[,,k],obs = observed_embed, silent = silent, maxit = 500, precondmethod = precondmethod, NNarray = NNarray)
    }

    avg_spec <- apply( specs[,,(n_iter-n_avg+1):n_iter], c(1,2), mean )
    if( ncondsim == 0 ){
        return(avg_spec)
    } else {
        condsims <- array( NA, c( nvec, ncondsim ) )
        for(j in 1:ncondsim ){
            condsims[,,j] <- condsim_spec(y0,avg_spec,observed_embed,silent = silent, maxit = 500, precondmethod = precondmethod, NNarray = NNarray)
        }
        return(list(avg_spec=avg_spec,condsims = condsims))
    }
}


#' @export
iterate_spec <- function(y, observed, embed_fac = 1.2, burn_iters = 100, par_spec_fun = spec_AR1,
                                 kern_parm = 10, precond_method = "fft", m = 10,
                                 silent = TRUE, max_iter = 200, tol = 1e-6){

    if(identical(par_spec_fun,FALSE)){
        do_parametric_filter <- FALSE
    } else {
        do_parametric_filter <- TRUE
    }

    # error for bad embedding factor value
    if( embed_fac < 1 ){
        stop("Embedding factor embed_fac cannot be less than 1")
    }

    # create embedded lattice
    nvec_obs <- dim(y)
    nvec <- round( nvec_obs*embed_fac )
    y_embed <- array(NA,nvec)
    y_embed[1:nvec_obs[1],1:nvec_obs[2]] <- y
    observed_embed <- array(FALSE,nvec)
    observed_embed[1:nvec_obs[1],1:nvec_obs[2]] <- observed

    if( precond_method == "Vecchia"){
        locsfull <- expand.grid( 1:nvec[1], 1:nvec[2] )
        locs <- locsfull[observed_embed,]
        NNarray <- aldodevel::findOrderedNN_kdtree(locs,m)
        NNarray[m+1,] <- (m+1):1 # need to make sure first part goes in the right order
    } else {
        NNarray <- NULL
    }

    # get grid size and define kerenel
    n <- prod(nvec)
    kern <- sqexp_kern(kern_parm, nvec)

    # impute with mean
    y0 <- y_embed
    y0[!observed_embed] <- mean(y, na.rm = TRUE)

    # define likelihood function (written for AR1 filter)
    likfun <- function(x){
        expitx <- expit(x)/4
        return(-whittle_lik(pgram,expitx,par_spec_fun))
    }

    # set an initial value for optimization
    parm <- 1/8
    logitpar <- logit(parm*4)

    # loop over number of iterations
    #specs <- array(NA,c(nvec,n_iter))

    for(k in 1:burn_iters){
        # periodogram
        pgram <- 1/n*abs( fft(y0) )^2
        # do the optimization
        if(do_parametric_filter){
            res <- optim(logitpar,likfun,method="Brent",lower=-6,upper=6)
            logitpar <- res$par         # update logit parameter
            parm <- expit(logitpar)/4   # update parameter
            # parametric approximation
            param_spec <- par_spec_fun(parm,nvec)
        } else {
            param_spec <- array(1,nvec)
        }
        # smooth ratio of periodogram to parametric approximation
        sm_pgram <- smooth_pgram(pgram/param_spec,kern)
        # estimate is product of smoothed ratio and parametric spec
        spec <- sm_pgram*param_spec
        # do a conditional simulation with the new spectrum
        y0 <- condsim_spec(y = y0, spec = spec, obs = observed_embed, silent = silent, maxit = 500, precondmethod = precond_method, NNarray = NNarray, tol = tol)
    }

    spec_old <- spec
    tolval <- 0.05

    for(k in 1:max_iter){
        pgram <- 1/n*abs( fft(y0) )^2
        # do the optimization
        if(do_parametric_filter){
            res <- optim(logitpar,likfun,method="Brent",lower=-6,upper=6)
            logitpar <- res$par         # update logit parameter
            parm <- expit(logitpar)/4   # update parameter
            # parametric approximation
            param_spec <- par_spec_fun(parm,nvec)
        } else {
            param_spec <- array(1,nvec)
        }
        # smooth ratio of periodogram to parametric approximation
        sm_pgram <- smooth_pgram(pgram/param_spec,kern)
        # estimate is product of smoothed ratio and parametric spec
        spec <- sm_pgram*param_spec
        # update the spectrum estimate
        spec_new <- (k-1)/k*spec_old + 1/k*spec
        # compare difference to tolerance
        spec_sd <- sqrt( smooth_pgram( spec_old^2, kern^2 ) )

        criterion <- max( abs(spec_new - spec_old)/spec_sd )
        if( k %% 10 == 0 ) cat(paste("Averaging Iteration",k,"Criterion =",round(criterion,4),"\n"))
        if( criterion < tolval ){
            break
        }

        spec_old <- spec_new
        y0 <- condsim_spec(y = y0, spec = spec_new, obs = observed_embed, silent = silent, maxit = 500, precondmethod = precond_method, NNarray = NNarray, tol = tol)

    }
    avg_spec <- spec_new

    # get estimate of the likelihood
    cat("Computing estimate of likelihood \n")
    locsfull <- as.matrix( expand.grid( 1:nvec[1], 1:nvec[2] ) )
    locs <- locsfull[observed_embed,]
    NNarray <- aldodevel::findOrderedNN_kdtree(locs,40)
    covarray <- 1/prod(dim(avg_spec))*Re( fft( avg_spec, inverse = TRUE ) )
    yvec <- y_embed[observed_embed]
    loglik <- vecchiaLik(covarray,yvec,locs,NNarray)

    # get conditional expectation
    condexp <- condexp_spec(y0,avg_spec,observed_embed, silent=silent,maxit=500,precondmethod = precond_method, NNarray = NNarray)

    # get conditional simulations
    ncondsim <- 1
    condsim_array <- array( NA, c(nvec_obs,ncondsim) )
    for(j in 1:ncondsim){
        cursim <- condsim_spec(y0,avg_spec,observed_embed,silent=silent,precondmethod=precond_method,NNarray=NNarray)
        condsim_array[,,j] <- cursim[1:nvec_obs[1],1:nvec_obs[2]]
    }

    return(list( spec = avg_spec, cov = covarray, loglik = loglik, condexp = condexp[1:nvec_obs[1],1:nvec_obs[2]], condsim = condsim_array) )
}



iterate_spec_semipar_old <- function(y, observed, n_iter = 100, par_spec_fun = spec_AR1,
                                 kern_parm = 10, n_avg = 30, n_intraclass = 40, precondmethod = "fft", ...){

    # pass NNarray through "..."

    if(identical(par_spec_fun,FALSE)){
        do_parametric_filter <- FALSE
    } else {
        do_parametric_filter <- TRUE
    }

    # get grid size and define kerenel
    nvec <- dim(observed)
    n <- prod(nvec)
    kern <- sqexp_kern(kern_parm, nvec)

    # impute with mean
    y0 <- y
    y0[!observed] <- mean(y, na.rm = TRUE)

    # define likelihood function (written for AR1 filter)
    likfun <- function(x){
        expitx <- expit(x)/4
        return(-whittle_lik(pgram,expitx,par_spec_fun))
    }

    # set an initial value for optimization
    parm <- 1/8
    logitpar <- logit(parm*4)

    # loop over number of iterations
    specs <- array(NA,c(nvec,n_iter))

    for(k in 1:n_iter){
        print(k)
        # periodogram
        pgram <- 1/n*abs( fft(y0) )^2

        # do the optimization
        if(do_parametric_filter){
            res <- optim(logitpar,likfun,method="Brent",lower=-6,upper=6)
            logitpar <- res$par         # update logit parameter
            parm <- expit(logitpar)/4   # update parameter

            # parametric approximation
            param_spec <- par_spec_fun(parm,nvec)
        } else {
            param_spec <- array(1,nvec)
        }
        #cat(paste("parameter",round(parm,4),"\n"))

        # smooth ratio of periodogram to parametric approximation
        sm_pgram <- smooth_pgram(pgram/param_spec,kern)
        # estimate is product of smoothed ratio and parametric spec
        specs[,,k] <- sm_pgram*param_spec

        # do a conditional simulation with the new spectrum
        avg_indices <- (max(1,k-n_avg+1)):k
        avg_spec <- apply( specs[,,avg_indices], c(1,2), mean)
        #y0 <- condsim_spec(y0,avg_spec,observed,silent = TRUE)
        y0 <- condsim_spec(y0,specs[,,k],observed, silent = FALSE, maxit = 500, precondmethod = precondmethod, ...)
    }

    var_asymptotic <- sum(kern^2)*avg_spec^2

    spec_intraclass <- array(NA,c(nvec,n_intraclass))

    condsims <- array( NA, c(nvec,n_intraclass) )
    for(k in 1:n_intraclass){

        y0 <- condsim_spec(y0,avg_spec,observed,silent=TRUE, maxit = 500, precondmethod = precondmethod, ...)
        condsims[,,k] <- y0
        pgram <- 1/n*abs( fft(y0) )^2

        # do the optimization
        do_optimization <- FALSE
        if(do_optimization){
            if(do_parametric_filter){
                res <- optim(logitpar,likfun,method="Brent",lower=-6,upper=6)
                logitpar <- res$par         # update logit parameter
                parm <- expit(logitpar)/4   # update parameter

                # parametric approximation
                param_spec <- par_spec_fun(parm,nvec)
            } else {
                param_spec <- array(1,nvec)
            }
        }

        # smooth ratio of periodogram to parametric approximation
        sm_pgram <- smooth_pgram(pgram/param_spec,kern)
        # estimate is product of smoothed ratio and parametric spec
        spec_intraclass[,,k] <- sm_pgram*param_spec
    }

    var_intraclass <- apply(spec_intraclass, c(1,2), var)


    return(list(estimate = avg_spec, var_asymptotic = var_asymptotic,
                var_intraclass = var_intraclass,
                condsim = condsims))

}

