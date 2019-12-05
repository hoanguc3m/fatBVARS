library(fatBVARS)
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
# plot(Chain1)

###########################################################################
datagen <- sim.VAR.novol(dist="Student")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Student", SV = F)
inits <- get_init(prior)
Chain2 <- BVAR.novol(y, K = 5, p = 2, dist = "Student", y0 = NULL, prior = prior, inits = inits)
# plot(Chain2)

###########################################################################
datagen <- sim.VAR.novol(dist="Hyper.Student")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Hyper.Student", SV = F)
inits <- get_init(prior)
Chain3 <- BVAR.novol(y, K = 5, p = 2, dist = "Hyper.Student", y0 = NULL, prior = prior, inits = inits)
# plot(Chain3)

###########################################################################
datagen <- sim.VAR.novol(dist="multiStudent")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="multiStudent", SV = F)
inits <- get_init(prior)
Chain4 <- BVAR.novol(y, K = 5, p = 2, dist = "multiStudent", y0 = NULL, prior = prior, inits = inits)
# plot(Chain4)

###########################################################################
datagen <- sim.VAR.novol(dist="Hyper.multiStudent")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Hyper.multiStudent", SV = F)
inits <- get_init(prior)
Chain5 <- BVAR.novol(y, K = 5, p = 2, dist = "Hyper.multiStudent", y0 = NULL, prior = prior, inits = inits)
# plot(Chain5)

###########################################################################
datagen <- sim.VAR.SV(dist="Gaussian")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Gaussian", SV = T)
inits <- get_init(prior)
Chain6 <- BVAR.SV(y, K = 5, p = 2, dist = "Gaussian", y0 = NULL, prior = prior, inits = inits)
# plot(Chain6)

###########################################################################
datagen <- sim.VAR.SV(dist="Student")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Student", SV = T)
inits <- get_init(prior)
Chain7 <- BVAR.SV(y, K = 5, p = 2, dist = "Student", y0 = NULL, prior = prior, inits = inits)
# plot(Chain7)

###########################################################################
datagen <- sim.VAR.SV(dist="Hyper.Student")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Hyper.Student", SV = T)
inits <- get_init(prior)
Chain8 <- BVAR.SV(y, K = 5, p = 2, dist = "Hyper.Student", y0 = NULL, prior = prior, inits = inits)
# plot(Chain8)

###########################################################################
datagen <- sim.VAR.SV(dist="multiStudent")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="multiStudent", SV = T)
inits <- get_init(prior, samples = 11000, burnin = 1000, thin = 10)
Chain9 <- BVAR.SV(y, K = 5, p = 2, dist = "multiStudent", y0 = NULL, prior = prior, inits = inits)
# plot(Chain9)

###########################################################################
datagen <- sim.VAR.SV(dist="Hyper.multiStudent")
y <- datagen$y
prior <- get_prior(y, p = 2, dist="Hyper.multiStudent", SV = T)
inits <- get_init(prior, samples = 11000, burnin = 1000, thin = 10)
Chain10 <- BVAR.SV(y, K = 5, p = 2, dist = "Hyper.multiStudent", y0 = NULL, prior = prior, inits = inits)
# plot(Chain10)
