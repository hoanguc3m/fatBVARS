library(fatBVARS)
library(bridgesampling)
# library(readr)
# FRED <- read_csv("~/Downloads/FRED.csv",
#                  col_types = cols(sasdate = col_date(format = "%m/%d/%Y")))
data(fred)
max_length <- 729
lagnum <- 12
# IndustrialInd <- diff(log(FRED$INDPRO), lag = lagnum)*100
# Inflation <- diff(log(FRED$PCEPI), lag = lagnum)*100
# Interest <- FRED$TB3MS
# SP500lgr <- diff(log(FRED$`S&P 500`), lag = lagnum)*100

IndustrialInd <- diff(FRED$INDPRO, lag = lagnum)/ head(FRED$INDPRO, max_length-lagnum) *100
Inflation <- diff(FRED$PCEPI, lag = lagnum)/ head(FRED$PCEPI, max_length-lagnum)*100
Interest <- FRED$TB3MS
SP500lgr <- diff(FRED$`S&P 500`, lag = 1)/ head(FRED$`S&P 500`, max_length-1)*100


time_id <- 634
date <- zoo::as.Date(FRED$sasdate)

date[time_id]
p <- 13
y0 = cbind(IndustrialInd[c(1:p)],
           Inflation[c(1:p)],
           Interest[c(1:p)],
           SP500lgr[c(1:p)])

y <- cbind(IndustrialInd[c((p+1):time_id)],
           Inflation[c((p+1):time_id)],
           Interest[c((p+1):time_id)],
           SP500lgr[c((p+1):time_id)])


###########################################################################
prior <- get_prior(y, p = 13, dist="Gaussian", SV = F)
inits <- get_init(prior,  samples = 11000, burnin = 1000, thin = 10)
Chain1 <- BVAR.novol(y, K = 4, p = 13, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)
# plot(Chain1)
K = 4; p = 13; M = K + p * K^2;numa = 0.5* K * (K-1)

lb <- rep(-Inf, ncol(Chain1$mcmc))
lb[(M+numa+1):(M+numa+K)] <- 0
ub <- rep(Inf, ncol(Chain1$mcmc))
names(lb) <- names(ub) <- colnames(Chain1$mcmc)

log_posterior(Chain1$mcmc[1,], Chain1)
bridge_result1 <- bridge_sampler(samples = Chain1$mcmc, log_posterior = log_posterior,
                                data = Chain1, lb = lb, ub = ub,
                                method = "normal", repetitions = 5,
                                maxiter = 900, silent = F)

# Important sampling
Brobdingnag::sum( exp( as.brob(bridge_result$q21 - bridge_result$q22))) / 500
# Harmonic mean
1/ (Brobdingnag::sum( exp( as.brob(bridge_result$q12 - bridge_result$q11))) / 500)
error_measures(bridge_result1)
###########################################################################
prior <- get_prior(y, p = 13, dist="Student", SV = F)
inits <- get_init(prior,  samples = 11000, burnin = 1000, thin = 10)
Chain2 <- BVAR.novol(y, K = 4, p = 13, dist = "Student", y0 = y0, prior = prior, inits = inits)
# plot(Chain2)

K = 4; p = 13; M = K + p * K^2;numa = 0.5* K * (K-1); t_max = nrow(Chain2$y)

lb <- rep(-Inf, ncol(Chain2$mcmc))
lb[(M+numa+1):(M+numa+K)] <- 0
lb[(M+numa+K+1):(M+numa+K+1)] <- 2
lb[(M+numa+K+2):(M+numa+K+t_max+1)] <- 0

ub <- rep(Inf, ncol(Chain2$mcmc))
ub[(M+numa+K+1):(M+numa+K+1)] <- 100

names(lb) <- names(ub) <- colnames(Chain2$mcmc)

log_posterior(Chain2$mcmc[1,], Chain2)

bridge_result2 <- bridge_sampler(samples = Chain2$mcmc, log_posterior = log_posterior,
                                 data = Chain2, lb = lb, ub = ub,
                                 method = "normal",
                                 maxiter = 1000, silent = F)
# Important sampling
Brobdingnag::sum( exp( as.brob(bridge_result2$q21 - bridge_result2$q22))) / 500
# Harmonic mean
1/ (Brobdingnag::sum( exp( as.brob(bridge_result2$q12 - bridge_result2$q11))) / 500)
error_measures(bridge_result2)
###########################################################################
prior <- get_prior(y, p = 13, dist="Hyper.Student", SV = F)
inits <- get_init(prior,  samples = 11000, burnin = 1000, thin = 10)
Chain3 <- BVAR.novol(y, K = 4, p = 13, dist = "Hyper.Student", y0 = y0, prior = prior, inits = inits)
# plot(Chain3)

###########################################################################
prior <- get_prior(y, p = 13, dist="multiStudent", SV = F)
inits <- get_init(prior,  samples = 11000, burnin = 1000, thin = 10)
Chain4 <- BVAR.novol(y, K = 4, p = 13, dist = "multiStudent", y0 = y0, prior = prior, inits = inits)
# plot(Chain4)

###########################################################################
prior <- get_prior(y, p = 13, dist="Hyper.multiStudent", SV = F)
inits <- get_init(prior,  samples = 11000, burnin = 1000, thin = 10)
Chain5 <- BVAR.novol(y, K = 4, p = 13, dist = "Hyper.multiStudent", y0 = y0, prior = prior, inits = inits)
# plot(Chain5)

###########################################################################
prior <- get_prior(y, p = 13, dist="Gaussian", SV = T)
inits <- get_init(prior,  samples = 11000, burnin = 1000, thin = 10)
Chain6 <- BVAR.SV(y, K = 4, p = 13, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)
# plot(Chain6)

###########################################################################
prior <- get_prior(y, p = 13, dist="Student", SV = T)
inits <- get_init(prior,  samples = 11000, burnin = 1000, thin = 10)
Chain7 <- BVAR.SV(y, K = 4, p = 13, dist = "Student", y0 = y0, prior = prior, inits = inits)
# plot(Chain7)

###########################################################################
prior <- get_prior(y, p = 13, dist="Hyper.Student", SV = T)
inits <- get_init(prior,  samples = 11000, burnin = 1000, thin = 10)
Chain8 <- BVAR.SV(y, K = 4, p = 13, dist = "Hyper.Student", y0 = y0, prior = prior, inits = inits)
# plot(Chain8)

###########################################################################
prior <- get_prior(y, p = 13, dist="multiStudent", SV = T)
inits <- get_init(prior,  samples = 11000, burnin = 1000, thin = 10)
Chain9 <- BVAR.SV(y, K = 4, p = 13, dist = "multiStudent", y0 = y0, prior = prior, inits = inits)
# plot(Chain9)

###########################################################################
prior <- get_prior(y, p = 13, dist="Hyper.multiStudent", SV = T)
inits <- get_init(prior,  samples = 11000, burnin = 1000, thin = 10)
Chain10 <- BVAR.SV(y, K = 4, p = 13, dist = "Hyper.multiStudent", y0 = y0, prior = prior, inits = inits)
# plot(Chain10)
