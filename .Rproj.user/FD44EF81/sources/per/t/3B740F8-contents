
#install.packages("npspec_0.1.0.tar.gz", repos = NULL, type = "source" )

library("npspec")
library("viridis")

implot <- fields::image.plot

# read in land surface temperature data
load("datasets/AllSatelliteTemps.RData")
tmpr <- matrix( all.sat.temps$MaskTemp, 500, 300 )



# get grid size
n1 <- nrow(tmpr)
n2 <- ncol(tmpr)
nvec_obs <- c(n1,n2)

# get pattern of missing values
y <- tmpr[1:nvec_obs[1],1:nvec_obs[2]]
observed <- !is.na(y)
nobs <- sum(observed)

# define locations
locs <- as.matrix( expand.grid( 1:nvec_obs[1], 1:nvec_obs[2] ) )

# filter out linear mean in lat and lon
lm1 <- lm( y[1:prod(nvec_obs)] ~ locs[,1] + locs[,2] )
y0 <- array( NA, nvec_obs )
y0[observed] <- lm1$residuals

# save the fitted values (for use in making prediction plots below)
fitted_values <- lm1$coefficients[1] + lm1$coefficients[2]*locs[,1] + lm1$coefficients[3]*locs[,2]

# plot the data
colpal <- magma(64)
par(mfrow=c(1,2))
implot(tmpr,col=colpal)
implot(y0[1:nvec_obs[1],1:nvec_obs[2]],col=colpal)



# fit using the kernel parameter that minimized CV prediction error
#kernparm <- kernvec[which.min( crossval_sd)]
kernparm <- .018
t1 <- proc.time()
fit <- iterate_spec(y0, observed, burn_iters = 30, par_spec_fun = spec_AR1, embed_fac = 1.2,
                        kern_parm = kernparm, precond_method = "Vecchia", m = 10,
                        silent = TRUE)
(proc.time() - t1)/60




#
# compute 30 conditional simulations
#
ncondsim <- 30

# pick an embedding factor and define embedded objects
embed_fac <- 1.2
nvec_embed <- round(nvec_obs*embed_fac)
observed_embed <- array(FALSE, nvec_embed)
observed_embed[1:nvec_obs[1],1:nvec_obs[2]] <- observed
y_embed <- array(NA, nvec_embed)
y_embed[1:nvec_obs[1],1:nvec_obs[2]] <- y0

# object to store conditional simulations
condsim_array <- array(NA, c(nvec_obs,ncondsim) )

# get locations and nearest neighbors for use in Vecchia preconditioner
locsfull <- as.matrix( expand.grid( 1:nvec_embed[1], 1:nvec_embed[2] ) )
locs <- locsfull[observed_embed,]
m <- 10
NNarray <- aldodevel::findOrderedNN_kdtree(locs,m)
NNarray[m+1,] <- (m+1):1 # need to make sure first part goes in the right order

# do the conditional simulations and add back in
# the fitted values from the regression
timevec <- rep(NA,ncondsim)
for(j in 1:ncondsim){
    t1 <- proc.time()
    cursim <- condsim_spec(y = y_embed, spec = fit$spec, obs = observed_embed,
                           precondmethod="Vecchia", NNarray=NNarray, silent = FALSE, tol = 1e-6)
    timevec[j] <- proc.time()[3] - t1[3]
    cursim_obs <- cursim[1:nvec_obs[1],1:nvec_obs[2]]
    cursim_obs[1:prod(nvec_obs)] <- cursim_obs[1:prod(nvec_obs)] + fitted_values
    condsim_array[,,j] <- cursim_obs
    print(j)
}
mean(timevec)


# get conditional expectation
condexp <- fit$condexp
condexp[1:prod(nvec_obs)] <- condexp[1:prod(nvec_obs)] + fitted_values

# compute estimate of conditional standard deviation
sumsq <- matrix(0, nvec_obs[1], nvec_obs[2] )
for(j in 1:ncondsim) sumsq <- sumsq + (condsim_array[,,j] - condexp)^2
condsd <- sqrt( sumsq/ncondsim )

# plot the data, conditional expectation, conditional sd,
# and three conditional simulations
cxall <- 0.8; lnall <- 0.5
par(mfrow=c(2,3),mar=c(1,1,2,5))
zlimits <- range(condsim_array)
implot(tmpr[,n2:1],zlim=zlimits,axes=FALSE,col=colpal)
title("")
mtext("Data",side=3,line=lnall,family="serif",cex=cxall)
box()
implot(condexp[,n2:1],zlim=zlimits,axes=FALSE,col=colpal)
title("")
mtext("Conditional Expectation",side=3,line=lnall,family="serif",cex=cxall)
box()
implot(condsd[,n2:1],zlim=c(0,max(condsd)),axes=FALSE,col=colpal)
title("")
mtext("Standard Deviation",side=3,line=lnall,family="serif",cex=cxall)
box()
for(j in 1:3){
    implot(condsim_array[,n2:1,j],zlim=zlimits,axes=FALSE,col=colpal)
    title("")
    mtext(paste("Conditional Simulation",j),side=3,line=lnall,family="serif",cex=cxall)
    box()
}


# plot the estimated spectrum
par(mar=c(2,2,2,4))
implot(log10(fit$spec),axes=FALSE,col=colpal)
box()
title("")
mtext("Log 10 Spectral Density Estimate",side=3,line=0.5,family="serif")
mtext("0",side=1,line=0.5,at=0,family="serif",cex=0.8)
mtext(1,side=1,line=0.5,at=1,family="serif",cex=0.8)
mtext("0",side=2,line=0.5,at=0,family="serif",cex=0.8)
mtext(1,side=2,line=0.5,at=1,family="serif",cex=0.8)




