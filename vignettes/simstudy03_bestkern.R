

library("npspec", lib.loc = "./Rpackages" )

# pick a random seed based on the clock
tm <- as.numeric(Sys.time())
seed <- 1e8 * (tm - floor(tm))
set.seed(seed); print(seed)


# read in arguments from command line
# kernel parameter 'kern1000'
# simulation scenario 'setting'
# twice the smoothness parameter 'twonu'
args=(commandArgs(TRUE))
for(i in 1:length(args)) eval(parse(text=args[[i]]))


# load in best kernels
load("bestkern03.RData")

# set covariance parameters
parms <- c(2,8,twonu/2,0)

if( setting == 1 ){
    nvec_obs <- c(50,50)
    observed <- array(TRUE,nvec_obs)
} else
if( setting == 2 ){
    nvec_obs <- c(80,80)
    pmissing <- 0.3
    observed <- array(TRUE,nvec_obs)
    l1 <- round( nvec_obs[1]*sqrt(pmissing)*nvec_obs[1]/nvec_obs[2] )
    l2 <- round( nvec_obs[2]*sqrt(pmissing)*nvec_obs[2]/nvec_obs[1] )
    s1 <- floor( nvec_obs[1]/2 - l1/2 ) + 1
    e1 <- (s1+l1-1)
    s2 <- floor( nvec_obs[2]/2 - l1/2 ) + 1
    e2 <- (s2+l2-1)
    observed[ s1:e1, s2:e2 ] <- FALSE
} else
if( setting == 3 ){
    nvec_obs <- c(80,80)
    n <- prod(nvec_obs)
    pmissing <- 0.3
    observed <- array( rbinom(n,1,1-pmissing),nvec_obs ) == 1
}

# number of simulations
nsims <- 10

# embedding values and filter options
embed_vec <- c(1.0,1.1,1.2,1.3)
par_filter_list <- list( FALSE, spec_AR1 )

# outputs are different sizes, so need to use lists
est <- vector("list", length = 18 )
for(j in c(1:4,11:12) ){  est[[j]] <- array(NA, c(round(embed_vec[1]*nvec_obs),nsims) ) }
for(j in c(5:6,13:14) ){  est[[j]] <- array(NA, c(round(embed_vec[2]*nvec_obs),nsims) ) }
for(j in c(7:8,15:16) ){  est[[j]] <- array(NA, c(round(embed_vec[3]*nvec_obs),nsims) ) }
for(j in c(9:10,17:18)){ est[[j]] <- array(NA, c(round(embed_vec[4]*nvec_obs),nsims) ) }


for( simnum in 1:nsims ){

cat(paste(simnum))

# simulate data
y_full <- uncond_sim_matern2D(nvec_obs,parms)

# remove missing values
y <- y_full
y[!observed] <- NA
n_obs <- sum(observed)


# estimate spectrum using Fuentes method (plug in zeros)
kernparm <- bestkern[1,setting,twonu]
kern <- sqexp_kern( kernparm, nvec_obs )
y0 <- y
y0[!observed] <- 0
pgram0 <- 1/n_obs*abs( fft(y0) )^2
est[[1]][,,simnum] <- smooth_pgram(pgram0,kern)


# estimate using a taper
taper_prop <- 0.05
if( setting == 1 ){
    taper1 <- spec.taper(x = rep(1,nvec_obs[1]), p = taper_prop)
    taper2 <- spec.taper(x = rep(1,nvec_obs[2]), p = taper_prop)
    taper <- outer(taper1,taper2)
} else
if( setting == 2 ){
    taper_out1 <- spec.taper(x = rep(1,nvec_obs[1]), p = taper_prop)
    taper_out2 <- spec.taper(x = rep(1,nvec_obs[2]), p = taper_prop)
    ntaper1 <- sum( taper_out1 < 1 )
    len1 <- length(s1:e1) + ntaper1
    taper_in1 <- spec.taper( x = rep(1,len1), p = 1/2*ntaper1/len1 )
    ntaper2 <- sum( taper_out2 < 1 )
    len2 <- length(s2:e2) + ntaper2
    taper_in2 <- spec.taper( x = rep(1,len2), p = 1/2*ntaper2/len2 )
    taper_out <- outer(taper_out1,taper_out2)
    taper_in <- 1-outer(taper_in1,taper_in2)
    taper <- taper_out
    taper[ (s1-ntaper1/2):(e1+ntaper1/2), (s2-ntaper2/2):(e2+ntaper2/2) ] <- taper_in
} else
if( setting == 3 ){
    taper1 <- spec.taper(x = rep(1,nvec_obs[1]), p = taper_prop)
    taper2 <- spec.taper(x = rep(1,nvec_obs[2]), p = taper_prop)
    taper <- outer(taper1,taper2)*observed
}

# tapered estimate
kernparm <- bestkern[2,setting,twonu]
kern <- sqexp_kern( kernparm, nvec_obs )
pgram0 <- 1/sum(taper^2)*abs( fft( y0*taper ) )^2
est[[2]][,,simnum] <- smooth_pgram(pgram0,kern)


indexmat1 <- matrix( 3:10,  4, 2, byrow = TRUE )
indexmat2 <- matrix( 11:18, 4, 2, byrow = TRUE )

cat("  ")
# estimate models using these settings
for(j1 in 1:length(embed_vec)){  # embedding size (0% or 10% or 20% or 30%))
    for(j2 in 1:length(par_filter_list)){ # parametric filter (yes or no)

        ind <- indexmat1[j1,j2]

        kernparm <- bestkern[ind,setting,twonu]
        sv <- iterate_spec(
            y, observed,
            embed_fac = embed_vec[j1],
            burn_iters = 100,
            par_spec_fun = par_filter_list[[j2]],
            kern_parm = kernparm,
            converge_tol = 0.01,
    	    precond_method = "Vecchia", silent = TRUE)

        est[[ind]][,,simnum] <- sv$spec


        ind <- indexmat2[j1,j2]
        kernparm <- bestkern[ind,setting,twonu]

        sv <- iterate_spec_not_periodic(
            y, observed,
            embed_fac = embed_vec[j1],
            burn_iters = 100,
            par_spec_fun = par_filter_list[[j2]],
            kern_parm = kernparm,
            converge_tol = 0.01,
    	    precond_method = "Vecchia", silent = TRUE)

        est[[ind]][,,simnum] <- sv$spec

        cat("*")
    }
}
cat("\n")
}

# compile the results
specden1 <- get_matern2D_spec(round(embed_vec[1]*nvec_obs),parms)
specden2 <- get_matern2D_spec(round(embed_vec[2]*nvec_obs),parms)
specden3 <- get_matern2D_spec(round(embed_vec[3]*nvec_obs),parms)
specden4 <- get_matern2D_spec(round(embed_vec[4]*nvec_obs),parms)

est_mean <- vector("list",length=18)
for( j in 1:length(est_mean)){
    est_mean[[j]] <- apply( est[[j]] , c(1,2), mean )
}

msefun <- function( specest , targetspec ){

    sumsq <- array(0,dim(targetspec))
    nsims <- dim(specest)[3]
    for(j in 1:nsims){
        sumsq <- sumsq + 1/nsims*( specest[,,j]/targetspec - 1 )^2
    }
    return(sumsq)
}

est_mse <- vector("list",length=18)
for(j in 1:4){   est_mse[[j]] <- msefun( est[[j]], specden1 )  }
for(j in 5:6){   est_mse[[j]] <- msefun( est[[j]], specden2 )  }
for(j in 7:8){   est_mse[[j]] <- msefun( est[[j]], specden3 )  }
for(j in 9:10){  est_mse[[j]] <- msefun( est[[j]], specden4 )  }
for(j in 11:12){   est_mse[[j]] <- msefun( est[[j]], specden1 )  }
for(j in 13:14){   est_mse[[j]] <- msefun( est[[j]], specden2 )  }
for(j in 15:16){   est_mse[[j]] <- msefun( est[[j]], specden3 )  }
for(j in 17:18){  est_mse[[j]] <- msefun( est[[j]], specden4 )  }


# save the seed and the estimated spectra to object and to file
fname <- paste0("sim03results_bestkern/result_set",setting,"_twonu",twonu,"_batch",batch,".RData")
save(seed, est_mean, est_mse, file = fname)

