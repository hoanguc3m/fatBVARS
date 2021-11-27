library(fatBVARS)
load("MonthlyData.RData")

###########################################################################

time_id = nrow(datafin)
y0 = datafin[c(1:p),]
y <- datafin[c((p+1):time_id),]

burn <- 0
reps <- 10000

times <- matrix(0,2,7)
colnames(times) <- c("Gaussian", "Stud-t", "Skew-t", "OT", "MT", "OST", "MST")
rownames(times) <- c("Const. var", "SV")
###########################################################################

  prior <- get_prior(y, p = p, dist="Gaussian", SV = F)
  inits <- get_init(prior, samples = reps, burnin = burn, thin = 5)
  ptm <- proc.time()
  Chain <- BVAR.novol(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)
  ptm <- proc.time() - ptm
  times[1,1] <- ptm[3]

###########################################################################

  prior <- get_prior(y, p = p, dist="Student", SV = F)
  inits <- get_init(prior, samples = reps, burnin = burn, thin = 5)
  ptm <- proc.time()
  Chain <- BVAR.novol(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
  ptm <- proc.time() - ptm
  times[1,2] <- ptm[3]

###########################################################################

  prior <- get_prior(y, p = p, dist="Skew.Student", SV = F)
  inits <- get_init(prior, samples = reps, burnin = burn, thin = 5)
  ptm <- proc.time()
  Chain <- BVAR.novol(y, K = K, p = p, dist = "Skew.Student", y0 = y0, prior = prior, inits = inits)
  ptm <- proc.time() - ptm
  times[1,3] <- ptm[3]

###########################################################################

  prior <- get_prior(y, p = p, dist="OT", SV = F)
  inits <- get_init(prior, samples = reps, burnin = burn, thin = 5)
  ptm <- proc.time()
  Chain <- BVAR.novol(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
  ptm <- proc.time() - ptm
  times[1,4] <- ptm[3]

###########################################################################

  prior <- get_prior(y, p = p, dist="MT", SV = F)
  inits <- get_init(prior, samples = reps, burnin = burn, thin = 5)
  ptm <- proc.time()
  Chain <- BVAR.novol(y, K = K, p = p, dist = "MT", y0 = y0, prior = prior, inits = inits)
  ptm <- proc.time() - ptm
  times[1,5] <- ptm[3]

###########################################################################

  prior <- get_prior(y, p = p, dist="OST", SV = F)
  inits <- get_init(prior, samples = reps, burnin = burn, thin = 5)
  ptm <- proc.time()
  Chain <- BVAR.novol(y, K = K, p = p, dist = "OST", y0 = y0, prior = prior, inits = inits)
  ptm <- proc.time() - ptm
  times[1,6] <- ptm[3]

###########################################################################

  prior <- get_prior(y, p = p, dist="MST", SV = F)
  inits <- get_init(prior, samples = reps, burnin = burn, thin = 5)
  ptm <- proc.time()
  Chain <- BVAR.novol(y, K = K, p = p, dist = "MST", y0 = y0, prior = prior, inits = inits)
  ptm <- proc.time() - ptm
  times[1,7] <- ptm[3]

###########################################################################

  prior <- get_prior(y, p = p, dist="Gaussian", SV = T)
  inits <- get_init(prior, samples = reps, burnin = burn, thin = 5)
  ptm <- proc.time()
  Chain <- BVAR.SV(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)
  ptm <- proc.time() - ptm
  times[2,1] <- ptm[3]

###########################################################################

  prior <- get_prior(y, p = p, dist="Student", SV = T)
  inits <- get_init(prior, samples = reps, burnin = burn, thin = 5)
  ptm <- proc.time()
  Chain <- BVAR.SV(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
  ptm <- proc.time() - ptm
  times[2,2] <- ptm[3]

###########################################################################

  prior <- get_prior(y, p = p, dist="Skew.Student", SV = T)
  inits <- get_init(prior, samples = reps, burnin = burn, thin = 5)
  ptm <- proc.time()
  Chain <- BVAR.SV(y, K = K, p = p, dist = "Skew.Student", y0 = y0, prior = prior, inits = inits)
  ptm <- proc.time() - ptm
  times[2,3] <- ptm[3]

###########################################################################

  prior <- get_prior(y, p = p, dist="OT", SV = T)
  inits <- get_init(prior, samples = reps, burnin = burn, thin = 5)
  ptm <- proc.time()
  Chain <- BVAR.SV(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
  ptm <- proc.time() - ptm
  times[2,4] <- ptm[3]

###########################################################################

  prior <- get_prior(y, p = p, dist="MT", SV = T)
  inits <- get_init(prior, samples = reps, burnin = burn, thin = 5)
  ptm <- proc.time()
  Chain <- BVAR.SV(y, K = K, p = p, dist = "MT", y0 = y0, prior = prior, inits = inits)
  ptm <- proc.time() - ptm
  times[2,5] <- ptm[3]

###########################################################################

  prior <- get_prior(y, p = p, dist="OST", SV = T)
  inits <- get_init(prior, samples = reps, burnin = burn, thin = 5)
  ptm <- proc.time()
  Chain <- BVAR.SV(y, K = K, p = p, dist = "OST", y0 = y0, prior = prior, inits = inits)
  ptm <- proc.time() - ptm
  times[2,6] <- ptm[3]

###########################################################################

  prior <- get_prior(y, p = p, dist="MST", SV = T)
  inits <- get_init(prior,  samples = reps, burnin = burn, thin = 5)
  ptm <- proc.time()
  Chain <- BVAR.SV(y, K = K, p = p, dist = "MST", y0 = y0, prior = prior, inits = inits)
  ptm <- proc.time() - ptm
  times[2,7] <- ptm[3]

###########################################################################

cat( "\n\n", reps, "draws with", burn, "burnin\n\nTimes in seconds\n")
print(times)
cat("\nTimes relative to Gaussian SV\n")
print(times/times[2,1])
