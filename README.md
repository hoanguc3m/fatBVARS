## Bayesian VAR with Stochastic volatility and fat tails 
### Package installation
```
devtools::install_github("hoanguc3m/fatBVARS")
```
### Simulation
```
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
```

