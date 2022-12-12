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
  prior <- get_prior(y, p = p, dist="Student", SV = F)
  inits <- get_init(prior, samples = 110000, burnin = 10000, thin = 10)
  Chain2 <- BVAR.novol(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
  save.image(paste(path, "MonthlyData_Chain2.RData",sep = ""))

###########################################################################
  prior <- get_prior(y, p = p, dist="Skew.Student", SV = F)
  inits <- get_init(prior, samples = 110000, burnin = 10000, thin = 10)
  Chain3 <- BVAR.novol(y, K = K, p = p, dist = "Skew.Student", y0 = y0, prior = prior, inits = inits)
  save.image(paste(path, "MonthlyData_Chain3.RData",sep = ""))

###########################################################################
  prior <- get_prior(y, p = p, dist="MT", SV = F)
  inits <- get_init(prior, samples = 110000, burnin = 10000, thin = 10)
  Chain4 <- BVAR.novol(y, K = K, p = p, dist = "MT", y0 = y0, prior = prior, inits = inits)
  save.image(paste(path, "MonthlyData_Chain4.RData",sep = ""))

###########################################################################
  prior <- get_prior(y, p = p, dist="MST", SV = F)
  inits <- get_init(prior, samples = 110000, burnin = 10000, thin = 10)
  Chain5 <- BVAR.novol(y, K = K, p = p, dist = "MST", y0 = y0, prior = prior, inits = inits)
  save.image(paste(path, "MonthlyData_Chain5.RData",sep = ""))

###########################################################################
  prior <- get_prior(y, p = p, dist="Gaussian", SV = T)
  inits <- get_init(prior, samples = 110000, burnin = 10000, thin = 10)
  Chain6 <- BVAR.SV(y, K = K, p = p, dist = "Gaussian", y0 = y0, prior = prior, inits = inits)
  save.image(paste(path, "MonthlyData_Chain6.RData",sep = ""))

###########################################################################
  prior <- get_prior(y, p = p, dist="Student", SV = T)
  inits <- get_init(prior, samples = 110000, burnin = 10000, thin = 10)
  Chain7 <- BVAR.SV(y, K = K, p = p, dist = "Student", y0 = y0, prior = prior, inits = inits)
  save.image(paste(path, "MonthlyData_Chain7.RData",sep = ""))

###########################################################################
  prior <- get_prior(y, p = p, dist="Skew.Student", SV = T)
  inits <- get_init(prior, samples = 110000, burnin = 10000, thin = 10)
  Chain8 <- BVAR.SV(y, K = K, p = p, dist = "Skew.Student", y0 = y0, prior = prior, inits = inits)
  save.image(paste(path, "MonthlyData_Chain8.RData",sep = ""))

###########################################################################
  prior <- get_prior(y, p = p, dist="MT", SV = T)
  inits <- get_init(prior, samples = 110000, burnin = 10000, thin = 10)
  Chain9 <- BVAR.SV(y, K = K, p = p, dist = "MT", y0 = y0, prior = prior, inits = inits)
  save.image(paste(path, "MonthlyData_Chain9.RData",sep = ""))

###########################################################################
  prior <- get_prior(y, p = p, dist="MST", SV = T)
  inits <- get_init(prior,  samples = 110000, burnin = 10000, thin = 10)
  Chain10 <- BVAR.SV(y, K = K, p = p, dist = "MST", y0 = y0, prior = prior, inits = inits)
  save.image(paste(path, "MonthlyData_Chain10.RData",sep = ""))

###########################################################################
  prior <- get_prior(y, p = p, dist="OT", SV = F)
  inits <- get_init(prior, samples = 110000, burnin = 10000, thin = 10)
  Chain11 <- BVAR.novol(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
  save.image(paste(path, "MonthlyData_Chain11.RData",sep = ""))

###########################################################################
  prior <- get_prior(y, p = p, dist="OST", SV = F)
  inits <- get_init(prior, samples = 110000, burnin = 10000, thin = 10)
  Chain12 <- BVAR.novol(y, K = K, p = p, dist = "OST", y0 = y0, prior = prior, inits = inits)
  save.image(paste(path, "MonthlyData_Chain12.RData",sep = ""))

###########################################################################
  prior <- get_prior(y, p = p, dist="OT", SV = T)
  inits <- get_init(prior, samples = 110000, burnin = 10000, thin = 10)
  Chain13 <- BVAR.SV(y, K = K, p = p, dist = "OT", y0 = y0, prior = prior, inits = inits)
  save.image(paste(path, "MonthlyData_Chain13.RData",sep = ""))

###########################################################################
  prior <- get_prior(y, p = p, dist="OST", SV = T)
  inits <- get_init(prior, samples = 110000, burnin = 10000, thin = 10)
  Chain14 <- BVAR.SV(y, K = K, p = p, dist = "OST", y0 = y0, prior = prior, inits = inits)
  save.image(paste(path, "MonthlyData_Chain14.RData",sep = ""))
