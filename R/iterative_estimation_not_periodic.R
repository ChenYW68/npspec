
#' @export
iterate_spec_not_periodic <- function(y, observed, X = NULL, embed_fac = 1.2, burn_iters = 100, par_spec_fun = spec_AR1,
                                 kern_parm = 2*pi/sqrt(sum(observed)), precond_method = "fft", m = 10,
                                 silent = TRUE, max_iter = 200, tol = 1e-6,
                                 converge_tol = 0.05,
                                 ncondsim = 1){

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

    nvec  <- round( nvec_obs*embed_fac )
    nvec_double <- 2*round( nvec_obs*embed_fac )

    onearray <- array(1,nvec)
    odds1 <- 2*(1:nvec[1])-1
    odds2 <- 2*(1:nvec[2])-1
    onearray_double <- array(0,nvec_double)
    onearray_double[odds1,odds2] <- onearray



    y_embed <- array(NA,nvec_double)
    y_embed[1:nvec_obs[1],1:nvec_obs[2]] <- y
    observed_embed <- array(FALSE,nvec_double)
    observed_embed[1:nvec_obs[1],1:nvec_obs[2]] <- observed

    # get setup for vecchia's approximation
    locsfull <- expand.grid( 1:nvec_double[1], 1:nvec_double[2] )
    locs <- locsfull[observed_embed,]
    NNarray <- findOrderedNN_kdtree(locs,m)
    NNarray[m+1,] <- (m+1):1 # need to make sure first part goes in the right order

    # get grid size and define kerenel
    n <- prod(nvec)
    kern <- sqexp_kern(kern_parm, nvec)
    kern2 <- sqexp_kern(0.01, nvec_double)
    normarray_double <- smooth_pgram( onearray_double, kern2, smoothlog = FALSE )


    # do least squares to get initial estimate of mean
    if( !is.null(X) ){
        if( is.na(dim(X)[3]) ){
            n_covar <- 1
            dim(X) <- c(dim(X),1)
        } else {
            n_covar <- dim(X)[3]
        }
        X_embed <- array(NA, c(nvec,n_covar))
        for(j in 1:n_covar){
            X_embed[1:nvec_obs[1], 1:nvec_obs[2],j] <- X[ , , j]
        }
        X_embed_vecs <- X_embed
        dim(X_embed_vecs) <- c( prod(dim(X_embed)[1:2]), n_covar )
        infomat <- crossprod( X_embed_vecs[observed_embed,] )
        xcrossy <- crossprod( X_embed_vecs[observed_embed,], y_embed[observed_embed] )
        betahat <- solve( infomat, xcrossy )
        muvec <- X_embed_vecs %*% betahat
        mumat <- array( muvec, nvec )
        y0 <- y_embed - mumat
    } else {
        y0 <- y_embed - mean( y_embed, na.rm = TRUE )
        mumat <- array( mean( y_embed, na.rm = TRUE ), nvec )
        betahat <- NA
        infomat <- NA
    }

    # impute with mean
    y0[!observed_embed] <- mean(y0, na.rm = TRUE)

    # define likelihood function (written for AR1 filter)
    likfun <- function(x){
        expitx <- expit(x)/4
        return(-whittle_lik(pgram,expitx,par_spec_fun))
    }

    # set an initial value for optimization
    parm <- 1/8
    logitpar <- logit(parm*4)

    # loop over number of iterations

    for(k in 1:burn_iters){

        # periodogram
        # take periodogram over embedded area, but not double embedded area
        pgram <- 1/n*abs( fft(y0[1:nvec[1],1:nvec[2]]) )^2
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
        sm_pgram <- smooth_pgram(pgram/param_spec,kern, smoothlog = FALSE)
        # estimate is product of smoothed ratio and parametric spec
        spec <- sm_pgram*param_spec

        # this is different for not periodic case!
        # we need to interpolate spectrum in order
        # to impute non-periodically on nvec_embed
        spec_double <- array(0,nvec_double)
        spec_double[odds1,odds2] <- spec
        spec_double <- smooth_pgram( spec_double, kern2, smoothlog = FALSE )/normarray_double


        # every 10th iteration, re-estimate the mean
        if( k %% 10 == 0 ){
            if(!is.null(X)){
                Xsolve <- array(NA, c(nvec,n_covar))
                for(j in 1:n_covar){
                    Xsolve[,,j] <- pcg_spec( X_embed[,,j], spec, observed_embed, precond_method, NNarray = NNarray, silent = silent )$x
                }
                Xsolve_vec <- Xsolve
                dim(Xsolve_vec) <- dim(X_embed_vecs)
                infomat <- crossprod( X_embed_vecs[observed_embed,], Xsolve_vec[observed_embed,] )
                xcrossy <- crossprod( Xsolve_vec[observed_embed,], y_embed[observed_embed] )
                betahat <- solve( infomat, xcrossy )
                muvec <- X_embed_vecs %*% betahat
                mumat <- array( muvec, nvec )
                y0 <- y_embed - mumat
            }
        }

        # do a conditional simulation with the new spectrum
        y0 <- condsim_spec(y = y0, spec = spec_double, obs = observed_embed, silent = silent, maxit = 500, precondmethod = precond_method, NNarray = NNarray, tol = tol)
    }

    spec_old <- spec

    for(k in 1:max_iter){
        pgram <- 1/n*abs( fft(y0[1:nvec[1],1:nvec[2]]) )^2
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
        sm_pgram <- smooth_pgram(pgram/param_spec,kern, smoothlog = FALSE)
        # estimate is product of smoothed ratio and parametric spec
        spec <- sm_pgram*param_spec

        # update the spectrum estimate
        spec_new <- (k-1)/k*spec_old + 1/k*spec

        # this is different for not periodic case!
        # we need to interpolate spectrum in order
        # to impute non-periodically on nvec_embed
        spec_double <- array(0,nvec_double)
        spec_double[odds1,odds2] <- spec_new
        spec_double <- smooth_pgram( spec_double, kern2, smoothlog = FALSE )/normarray_double


        # compare difference to tolerance
        spec_sd <- sqrt( smooth_pgram( spec_old^2, kern^2 ) )

        criterion <- max( abs(spec_new - spec_old)/spec_sd )
        #if( k %% 10 == 0 ) cat(paste("Averaging Iteration",k,"Criterion =",round(criterion,4),"\n"))
        if( criterion < converge_tol ){
            break
        }

        # every 10th iteration, re-estimate the mean
        if( k %% 10 == 0 ){
            if(!is.null(X)){
                Xsolve <- array(NA, c(nvec,n_covar))
                for(j in 1:n_covar){
                    Xsolve[,,j] <- pcg_spec( X_embed[,,j], spec_new, observed_embed, precond_method, NNarray = NNarray, silent = silent )$x
                }
                Xsolve_vec <- Xsolve
                dim(Xsolve_vec) <- dim(X_embed_vecs)
                infomat <- crossprod( X_embed_vecs[observed_embed,], Xsolve_vec[observed_embed,] )
                xcrossy <- crossprod( Xsolve_vec[observed_embed,], y_embed[observed_embed] )
                betahat <- solve( infomat, xcrossy )
                muvec <- X_embed_vecs %*% betahat
                mumat <- array( muvec, nvec )
                y0 <- y_embed - mumat
            }
        }

        spec_old <- spec_new
        y0 <- condsim_spec(y = y0, spec = spec_double, obs = observed_embed, silent = silent, maxit = 500, precondmethod = precond_method, NNarray = NNarray, tol = tol)

    }
    avg_spec <- spec_new

    # get estimate of the likelihood
    #cat("Computing estimate of likelihood \n")
    #locsfull <- as.matrix( expand.grid( 1:nvec[1], 1:nvec[2] ) )
    #locs <- locsfull[observed_embed,]
    #NNarray <- findOrderedNN_kdtree(locs,40)
    #covarray <- 1/prod(dim(avg_spec))*Re( fft( avg_spec, inverse = TRUE ) )
    #yvec <- y_embed[observed_embed]
    #loglik <- vecchiaLik(covarray,yvec,locs,NNarray)

    # get conditional expectation
    #cat("Computing Conditional Simulations \n")
    #condexp <- condexp_spec(y0,avg_spec,observed_embed, silent=silent,
    #    maxit=500,precondmethod = precond_method, NNarray = NNarray) +
    #    mumat

    # get conditional simulations
    #condsim_array <- array( NA, c(nvec_obs,ncondsim) )
    #condsim_array <- array( NA, c(nvec,ncondsim) )
    #for(j in 1:ncondsim){
    #    cursim <- condsim_spec(y0,avg_spec,observed_embed,silent=silent,
    #        precondmethod=precond_method,NNarray=NNarray) +
    #        mumat
    #    condsim_array[,,j] <- cursim[1:nvec_obs[1],1:nvec_obs[2]]
        #condsim_array[,,j] <- cursim
    #}

    return(list( spec = avg_spec ))
    #, cov = covarray, loglik = loglik,
    #    condexp = condexp[1:nvec_obs[1],1:nvec_obs[2]], condsim = condsim_array,
    #    betahat = betahat, mumat = mumat[1:nvec_obs[1],1:nvec_obs[2]], betacov = solve(infomat)) )
}
