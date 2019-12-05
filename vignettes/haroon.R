library(fatBVARS)
# library(readr)
# FRED <- read_csv("~/Downloads/FRED.csv",
#                  col_types = cols(sasdate = col_date(format = "%m/%d/%Y")))
data(fred)
max_length <- 729
lagnum <- 12
IndustrialInd <- diff(log(FRED$INDPRO), lag = lagnum)*100
Inflation <- diff(log(FRED$PCEPI), lag = lagnum)*100
Interest <- FRED$TB3MS
SP500lgr <- diff(log(FRED$`S&P 500`), lag = lagnum)*100
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

###########################################################################
prior <- get_prior(y, p = 13, dist="Student", SV = F)
inits <- get_init(prior,  samples = 11000, burnin = 1000, thin = 10)
Chain2 <- BVAR.novol(y, K = 4, p = 13, dist = "Student", y0 = y0, prior = prior, inits = inits)
# plot(Chain2)

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
