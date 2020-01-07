library(fatBVARS)
library(bridgesampling)
###########################################################################
# Default values for simulation
# K = 5, p = 2, t_max = 1000,
# b0 = 0.6, a0 = 0.5, h = 0, nu = 6, gamma = 0.5,
# seednum = 0, burn_in = 0

###########################################################################
datagen <- sim.VAR.SV(dist="Hyper.multiStudent")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Hyper.multiStudent", SV = T)
inits <- get_init(prior, samples = 12000, burnin = 2000, thin = 1)
Chain10 <- BVAR.SV(y, K = 5, p = 2, dist = "Hyper.multiStudent", y0 = NULL, prior = prior, inits = inits)
# plot(Chain10)
lub <- get_lu_param(Chain10$mcmc)
lb <- lub$lb
ub <- lub$ub

log_posterior(Chain10$mcmc[1,], Chain10)
waic_chain10 <- WAIC(Chain10)
recur_chain10 <- recursive_forecast(Chain = Chain10, y0 = y[1:500,], y_future = y[501:1000,], t_pred = 12, reestimated = F)
###########################################################################

prior <- get_prior(y, p = 2, dist="Gaussian", SV = F)
inits <- get_init(prior, samples = 12000, burnin = 2000, thin = 1)
Chain1 <- BVAR.novol(y, K = 5, p = 2, dist = "Gaussian", y0 = NULL, prior = prior, inits = inits)
waic_chain1 <- WAIC(Chain1)

# plot(Chain1)
# plot(Chain1, element = "b")
lub <- get_lu_param(Chain1$mcmc)
lb <- lub$lb
ub <- lub$ub


log_posterior(Chain1$mcmc[1,], Chain1)
bridge_result1 <- bridge_sampler(samples = Chain1$mcmc, log_posterior = log_posterior,
                                 data = Chain1, lb = lb, ub = ub,
                                 method = "normal",
                                 maxiter = 1000, silent = F)

###########################################################################

prior <- get_prior(y, p = 2, dist="Student", SV = F)
inits <- get_init(prior, samples = 12000, burnin = 2000, thin = 1)
Chain2 <- BVAR.novol(y, K = 5, p = 2, dist = "Student", y0 = NULL, prior = prior, inits = inits)
waic_chain2 <- WAIC(Chain2)

# plot(Chain2)
lub <- get_lu_param(Chain2$mcmc)
lb <- lub$lb
ub <- lub$ub

log_posterior(Chain2$mcmc[1,], Chain2)

bridge_result2 <- bridge_sampler(samples = Chain2$mcmc, log_posterior = log_posterior,
                                 data = Chain2, lb = lb, ub = ub,
                                 method = "normal",
                                 maxiter = 1000, silent = F)

###########################################################################

prior <- get_prior(y, p = 2, dist="Hyper.Student", SV = F)
inits <- get_init(prior, samples = 12000, burnin = 2000, thin = 1)
Chain3 <- BVAR.novol(y, K = 5, p = 2, dist = "Hyper.Student", y0 = NULL, prior = prior, inits = inits)
waic_chain3 <- WAIC(Chain3)

# plot(Chain3)
lub <- get_lu_param(Chain3$mcmc)
lb <- lub$lb
ub <- lub$ub

log_posterior(Chain3$mcmc[1,], Chain3)

bridge_result3 <- bridge_sampler(samples = Chain3$mcmc, log_posterior = log_posterior,
                                 data = Chain3, lb = lb, ub = ub,
                                 method = "normal",
                                 maxiter = 1000, silent = F)
###########################################################################

prior <- get_prior(y, p = 2, dist="multiStudent", SV = F)
inits <- get_init(prior, samples = 12000, burnin = 2000, thin = 1)
Chain4 <- BVAR.novol(y, K = 5, p = 2, dist = "multiStudent", y0 = NULL, prior = prior, inits = inits)
waic_chain4 <- WAIC(Chain4)

# plot(Chain4)
lub <- get_lu_param(Chain4$mcmc)
lb <- lub$lb
ub <- lub$ub

log_posterior(Chain4$mcmc[1,], Chain4)

bridge_result4 <- bridge_sampler(samples = Chain4$mcmc, log_posterior = log_posterior,
                                 data = Chain4, lb = lb, ub = ub,
                                 method = "normal",
                                 maxiter = 1000, silent = F)

###########################################################################

prior <- get_prior(y, p = 2, dist="Hyper.multiStudent", SV = F)
inits <- get_init(prior, samples = 12000, burnin = 2000, thin = 1)
Chain5 <- BVAR.novol(y, K = 5, p = 2, dist = "Hyper.multiStudent", y0 = NULL, prior = prior, inits = inits)
waic_chain5 <- WAIC(Chain5)

# plot(Chain5)
lub <- get_lu_param(Chain5$mcmc)
lb <- lub$lb
ub <- lub$ub

log_posterior(Chain5$mcmc[1,], Chain5)

bridge_result5 <- bridge_sampler(samples = Chain5$mcmc, log_posterior = log_posterior,
                                 data = Chain5, lb = lb, ub = ub,
                                 method = "normal",
                                 maxiter = 1000, silent = F)
# Important sampling
Brobdingnag::sum( exp( as.brob(bridge_result5$q21 - bridge_result5$q22))) / 500
# Harmonic mean
1/ (Brobdingnag::sum( exp( as.brob(bridge_result5$q12 - bridge_result5$q11))) / 500)
error_measures(bridge_result5)

###########################################################################

prior <- get_prior(y, p = 2, dist="Gaussian", SV = T)
inits <- get_init(prior, samples = 12000, burnin = 2000, thin = 1)
Chain6 <- BVAR.SV(y, K = 5, p = 2, dist = "Gaussian", y0 = NULL, prior = prior, inits = inits)
waic_chain6 <- WAIC(Chain6)

# plot(Chain6)
lub <- get_lu_param(Chain6$mcmc)
lb <- lub$lb
ub <- lub$ub

log_posterior(Chain6$mcmc[1,], Chain6)
bridge_result6 <- bridge_sampler(samples = Chain6$mcmc, log_posterior = log_posterior,
                                 data = Chain6, lb = lb, ub = ub,
                                 method = "normal",
                                 maxiter = 1000, silent = F)

###########################################################################

prior <- get_prior(y, p = 2, dist="Student", SV = T)
inits <- get_init(prior, samples = 12000, burnin = 2000, thin = 1)
Chain7 <- BVAR.SV(y, K = 5, p = 2, dist = "Student", y0 = NULL, prior = prior, inits = inits)
waic_chain7 <- WAIC(Chain7)

# plot(Chain7)
lub <- get_lu_param(Chain7$mcmc)
lb <- lub$lb
ub <- lub$ub

log_posterior(Chain7$mcmc[1,], Chain7)

###########################################################################

prior <- get_prior(y, p = 2, dist="Hyper.Student", SV = T)
inits <- get_init(prior, samples = 12000, burnin = 2000, thin = 1)
Chain8 <- BVAR.SV(y, K = 5, p = 2, dist = "Hyper.Student", y0 = NULL, prior = prior, inits = inits)
waic_chain8 <- WAIC(Chain8)

# plot(Chain8)
lub <- get_lu_param(Chain8$mcmc)
lb <- lub$lb
ub <- lub$ub

log_posterior(Chain8$mcmc[1,], Chain8)

###########################################################################

prior <- get_prior(y, p = 2, dist="multiStudent", SV = T)
inits <- get_init(prior, samples = 12000, burnin = 2000, thin = 1)
Chain9 <- BVAR.SV(y, K = 5, p = 2, dist = "multiStudent", y0 = NULL, prior = prior, inits = inits)
waic_chain9 <- WAIC(Chain9)

# plot(Chain9)
lub <- get_lu_param(Chain9$mcmc)
lb <- lub$lb
ub <- lub$ub

log_posterior(Chain9$mcmc[1,], Chain9)

matrix(c(waic_chain1$waic, waic_chain2$waic, waic_chain3$waic, waic_chain4$waic, waic_chain5$waic,
         waic_chain6$waic, waic_chain7$waic, waic_chain8$waic, waic_chain9$waic, waic_chain10$waic), ncol = 2)

