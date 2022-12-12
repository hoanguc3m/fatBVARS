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

### Empirical illustration
See **vignettes** folder for more examples

```
library(fatBVARS)
data(MonthlyData)
view(datafin)

K <- ncol(datafin)
p <- 4
Time <- Monthly$sasdate[133:732]
numCores = 16
ndraws = 100000

###########################################################################

path = ""
paste(path, "MonthlyData.RData",sep = "")
#save.image(paste(path, "MonthlyData.RData",sep = ""))

time_id = nrow(datafin)
y0 = datafin[c(1:p),]
time_id <- nrow(datafin)
y <- datafin[c((p+1):time_id),]

cat("time_id", time_id, "\n")
###########################################################################
prior <- get_prior(y, p = p, dist="Gaussian", SV = F)
inits <- get_init(prior, samples = 110000, burnin = 10000, thin = 10)
Chain1 <- BVAR.novol(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)
paste(path, "MonthlyData_Chain1.RData",sep = "")
save.image(paste(path, "MonthlyData_Chain1.RData",sep = ""))

###########################################################################

```
### References

Karlsson, S., Mazur, S., & Nguyen, H. (2022). Vector autoregression models with skewness and heavy tails. Journal of Economic Dynamics and Control, 104580.