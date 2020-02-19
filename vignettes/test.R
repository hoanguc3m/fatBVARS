library(fatBVARS)
library(bridgesampling)
K = 5; p = 2
###########################################################################
# Default values for simulation
# K = 5, p = 2, t_max = 1000,
# b0 = 0.6, a0 = 0.5, h = 0, nu = 6, gamma = 0.5,
# seednum = 0, burn_in = 0

datagen <- sim.VAR.novol(dist="Gaussian")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Gaussian", SV = F)
inits <- get_init(prior)
Chain1 <- BVAR.novol(y, K = 5, p = 2, dist = "Gaussian", y0 = NULL, prior = prior, inits = inits)
ML_chain1 <- marginalLL(Chain1)
waic_chain1 <- WAIC(Chain1)
recursive_model1 <- recursive_forecast(y = y, t_start = 750, t_pred = 12, K = K, p = p, dist="Gaussian", SV = F)
recursive_model1v1 <- recursive_forecast(y = head(y,112), t_start = 100, t_pred = 12, K = K, p = p, dist="Gaussian", SV = F)
recursive_model1v2 <- recursive_forecast(y = head(y,113), t_start = 100, t_pred = 12, K = K, p = p, dist="Gaussian", SV = F)
recursive_model1v3 <- recursive_forecast(y = head(y,114), t_start = 100, t_pred = 12, K = K, p = p, dist="Gaussian", SV = F)
###########################################################################
datagen <- sim.VAR.novol(dist="Student")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Student", SV = F)
inits <- get_init(prior)
Chain2 <- BVAR.novol(y, K = 5, p = 2, dist = "Student", y0 = NULL, prior = prior, inits = inits)
ML_chain2 <- marginalLL(Chain2)
waic_chain2 <- WAIC(Chain2)

###########################################################################
datagen <- sim.VAR.novol(dist="Hyper.Student")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Hyper.Student", SV = F)
inits <- get_init(prior)
Chain3 <- BVAR.novol(y, K = 5, p = 2, dist = "Hyper.Student", y0 = NULL, prior = prior, inits = inits)
ML_chain3 <- marginalLL(Chain3)
waic_chain3 <- WAIC(Chain3)

###########################################################################
datagen <- sim.VAR.novol(dist="multiStudent")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="multiStudent", SV = F)
inits <- get_init(prior)
Chain4 <- BVAR.novol(y, K = 5, p = 2, dist = "multiStudent", y0 = NULL, prior = prior, inits = inits)
ML_chain4 <- marginalLL(Chain4)
waic_chain4 <- WAIC(Chain4)

###########################################################################
datagen <- sim.VAR.novol(dist="Hyper.multiStudent")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Hyper.multiStudent", SV = F)
inits <- get_init(prior)
Chain5 <- BVAR.novol(y, K = 5, p = 2, dist = "Hyper.multiStudent", y0 = NULL, prior = prior, inits = inits)
ML_chain5 <- marginalLL(Chain5)
waic_chain5 <- WAIC(Chain5)

###########################################################################
datagen <- sim.VAR.SV(dist="Gaussian")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Gaussian", SV = T)
inits <- get_init(prior)
Chain6 <- BVAR.SV(y, K = 5, p = 2, dist = "Gaussian", y0 = NULL, prior = prior, inits = inits)
ML_chain6 <- marginalLL(Chain6)
waic_chain6 <- WAIC(Chain6)

###########################################################################
datagen <- sim.VAR.SV(dist="Student")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Student", SV = T)
inits <- get_init(prior)
Chain7 <- BVAR.SV(y, K = 5, p = 2, dist = "Student", y0 = NULL, prior = prior, inits = inits)
ML_chain7 <- marginalLL(Chain7)
waic_chain7 <- WAIC(Chain7)

###########################################################################
datagen <- sim.VAR.SV(dist="Hyper.Student")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Hyper.Student", SV = T)
inits <- get_init(prior)
Chain8 <- BVAR.SV(y, K = 5, p = 2, dist = "Hyper.Student", y0 = NULL, prior = prior, inits = inits)
ML_chain8 <- marginalLL(Chain8)
waic_chain8 <- WAIC(Chain8)

###########################################################################
datagen <- sim.VAR.SV(dist="multiStudent")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="multiStudent", SV = T)
inits <- get_init(prior, samples = 1100, burnin = 100, thin = 1)
Chain9 <- BVAR.SV(y, K = 5, p = 2, dist = "multiStudent", y0 = NULL, prior = prior, inits = inits)
ML_chain9 <- marginalLL(Chain9)
waic_chain9 <- WAIC(Chain9)

###########################################################################
datagen <- sim.VAR.SV(dist="Hyper.multiStudent")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Hyper.multiStudent", SV = T)
inits <- get_init(prior, samples = 1100, burnin = 100, thin = 1)
Chain10 <- BVAR.SV(y, K = 5, p = 2, dist = "Hyper.multiStudent", y0 = NULL, prior = prior, inits = inits)
ML_chain10 <- marginalLL(Chain10)
waic_chain10 <- WAIC(Chain10)
recursive_model10 <- recursive_forecast(y = y, t_start = 985, t_pred = 12, K = K, p = p, dist="Hyper.multiStudent", SV = T)
###########################################################################

datagen <- sim.VAR.novol(dist="multiOrthStudent")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="multiStudent", SV = F)
inits <- get_init(prior, samples = 1100, burnin = 100, thin = 1)
Chain11 <- BVAR.novol(y, K = 5, p = 2, dist = "multiOrthStudent", y0 = NULL, prior = prior, inits = inits)
ML_chain11 <- marginalLL(Chain11)
###########################################################################
datagen <- sim.VAR.novol(dist="Hyper.multiOrthStudent")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Hyper.multiStudent", SV = F)
inits <- get_init(prior, samples = 1100, burnin = 100, thin = 1)
Chain12 <- BVAR.novol(y, K = 5, p = 2, dist = "Hyper.multiOrthStudent", y0 = NULL, prior = prior, inits = inits)
ML_chain12 <- marginalLL(Chain12)
###########################################################################
datagen <- sim.VAR.SV(dist="multiOrthStudent")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="multiStudent", SV = T)
inits <- get_init(prior, samples = 1100, burnin = 100, thin = 1)
Chain13 <- BVAR.SV(y, K = 5, p = 2, dist = "multiOrthStudent", y0 = NULL, prior = prior, inits = inits)
ML_chain13 <- marginalLL(Chain13)
###########################################################################
datagen <- sim.VAR.SV(dist="Hyper.multiOrthStudent")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Hyper.multiStudent", SV = T)
inits <- get_init(prior, samples = 1100, burnin = 100, thin = 1)
Chain14 <- BVAR.SV(y, K = 5, p = 2, dist = "Hyper.multiOrthStudent", y0 = NULL, prior = prior, inits = inits)
ML_chain14 <- marginalLL(Chain14)
