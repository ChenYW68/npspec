
#install.packages("../../npspec_0.1.0.tar.gz", repos = NULL, type = "source" )
library("npspec")

# read in land surface temperature data and put in a matrix
load("./datasets/SimulatedTemps.RData")
tmpr <- matrix( sim.data$Temp, 500, 300 )

# get grid size
n1 <- nrow(tmpr)
n2 <- ncol(tmpr)
nvec_obs <- c(n1,n2)

# get pattern of missing values
y <- tmpr[1:nvec_obs[1],1:nvec_obs[2]]
observed <- !is.na(y)
nobs <- sum(observed)

# define locations and covariates
locs <- as.matrix( expand.grid( 1:nvec_obs[1], 1:nvec_obs[2] ) )
X <- array(NA, c(nvec_obs,3))
X[,,1] <- 1
X[,,2] <- array( locs[,1], nvec_obs)
X[,,3] <- array( locs[,2], nvec_obs)

# fit the model
t1 <- proc.time()
fit <- iterate_spec(y, observed, X = X, burn_iters = 30, par_spec_fun = spec_AR1, embed_fac = 1.2,
        precond_method = "Vecchia", m = 10,
        silent = TRUE, ncondsim = 50)
(proc.time() - t1)/60

# predictions
pred_mat <- fit$condexp
pred_vec <- c(pred_mat)

# calculate the prediction variances based on the
# conditional simulations
cond_diff <- array(NA, dim(fit$condsim) )
for(j in 1:dim(fit$condsim)[3]) cond_diff[,,j] <- fit$condsim[,,j] - fit$condexp
meansq <- function(x) 1/length(x)*sum(x^2)
predvar_mat <- apply( cond_diff, c(1,2), meansq )
predvar_vec <- c(predvar_mat)

# plots
par(mfrow=c(1,2))
fields::image.plot(pred_mat)
fields::image.plot(predvar_mat)


