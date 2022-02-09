##########################################################################
# Adatptive metropolis hasting #
##########################################################################
TARGACCEPT = 0.3
batchlength = 10

#' @export
adaptamount <- function(iteration){
  return( min( 0.01, 1.0 / sqrt(iteration) ) );
}

# This function is borrowed from https://github.com/FK83/bvarsv
#' @export
getmix <- function(){
  # 7 components
  q <- c(0.00730, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.25750)      # probabilities
  m <- c(-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -1.08819) # means
  u2 <- c(5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261)    #variances

  # # 10 components
  # q <- c(0.00609, 0.04775, 0.13057, 0.20674, 0.22715, 0.18842, 0.12047, 0.05591, 0.01575, 0.00115)      # probabilities
  # m <- c(3.13417, 2.55484, 1.94244, 1.23006, 0.35567, -0.76538, -2.26048, -4.34506, -7.47644, -13.44260)
  # u2 <- c(0.11265, 0.17788, 0.26768, 0.40611, 0.62699, 0.98583, 1.57469, 2.54498, 4.16591, 7.33342)    #variances
  # # m <- c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000) # means in the Omori paper
  return(list(q=q,m=m,u2=u2))
}
##########################################################################

##########################################################################
# Utility functions  #
##########################################################################
#' @export
makeRegressor <- function(y, y0, t_max, K, p){
  xt <- NULL
  yt = t(y)

  if (is.null(y0)){
    y0 <- matrix(0, ncol = K, nrow = p)

  }
  y0 <- rbind(y0, y)
  for (i in c(1:p)){
    xt <- cbind(xt, rbind(1, vec( as.matrix( t(y0)[,(p+i-1):i]))) )
  }
  for (i in c( p:(t_max-1))){
    xt <- cbind(xt, rbind(1, vec( as.matrix( yt[,i:(i-p+1)]))) )
  }
  return(xt)
}

#' @export
Vec_to_Mat <- function(b, K, p){
  return(matrix(b, nrow = K))
}

#' @export
reprow<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

#' @export
repcol<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

#' @export
sample_A_ele <- function(ysub, xsub, a_sub, V_a_sub){
  t_max = length(ysub)
  n = nrow(xsub)

  a_post = rep(0, n)
  V_a_sub_inv = solve(V_a_sub)

  V_a_post_inv <- V_a_sub_inv + xsub %*% t(xsub)
  V_a_post <- solve(V_a_post_inv)
  a_post <- V_a_post %*% ( V_a_sub_inv %*% a_sub + xsub %*% ysub)
  return(a_post  + t(chol(V_a_post)) %*% rnorm(n))
}

#' @export
sample_h_ele <- function(ytilde, sigma_h = 0.0001*diag(K), h0_mean = rep(0,K),
                         h = matrix(0, nrow = t_max, ncol = K), K, t_max, prior){
  tmp <- getmix()
  q <- tmp$q
  m_mean <- tmp$m
  u2 <- tmp$u2
  # {
  #   Zs <- matrix(1,t_max,1) %x% diag(K)
  #   sigma_prmean <- h0_mean # mean h_0
  #   sigma_prvar <- 4*diag(K)   # variance h_0
  #
  #
  #   aux <- sigmahelper4(t(ytilde^2), q, m_mean, u2, h, Zs, sigma_h, sigma_prmean, sigma_prvar)
  #   h <- aux$Sigtdraw
  #   h0 <- as.numeric(aux$h0)
  # }

  {
    h0 <- rep(0,K)
    cond_var_sigma <- rep(0,K)
    cond_mean_sigma <- rep(0,K)
    Zs <- matrix(1,t_max,1) %x% diag(1)
    for (i in c(1:K)){
      sigma_prmean <- h0_mean[i] # mean h_0
      sigma_prvar <- matrix(4)   # variance h_0
      aux <- sigmahelper4(t(ytilde[ i,, drop =FALSE]^2), q, m_mean, u2, h[ i,, drop =FALSE], Zs, matrix(sigma_h[i,i]), sigma_prmean, sigma_prvar)
      h[i,] <- aux$Sigtdraw
      h0[i] <- as.numeric(aux$h0)
      # vart <- aux$vart # Variance of components
      # yss1 <- aux$yss1 # demean of components
      # h_tilde <- (h[i,1:t_max] - h0[i])/sqrt(sigma_h[i,i])
      # cond_var_sigma[i] <- 1/( 1 + sum(h_tilde^2/vart))
      # cond_mean_sigma[i] <- cond_var_sigma[i] * ( 0 + sum(h_tilde * yss1 / vart))
    }
  }

  # sqrtvol <- aux$sigt
  # [TODO] fix this
  # sse_2 <- apply( (h[,2:t_max] - h[,1:(t_max-1)])^2, MARGIN = 1, FUN = sum)
  if (K>1) {
    sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
  } else {
    sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
  }

  # # Normal prior
  # # Equation 9 in https://doi.org/10.1016/j.csda.2013.01.002
  # sigma_post_a <- rep(t_max,K) # prior of sigma_h Gamma(1,0.0001)
  # sigma_post_b <- sse_2 # prior of sigma_h
  #
  # for (i in c(1:K)){
  #   sigma_new <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
  #   alpha = (sigma_h[i,i] - sigma_new) / 2 / prior$sigma_S0 + 0.5 * (log(sigma_new) - log(sigma_h[i,i])) # B_sigma = 1
  #   temp = log(runif(1))
  #   if (alpha > temp){
  #     sigma_h[i,i] <- sigma_new
  #   }
  #   #log_sigma_den[]
  # }

  sigma_h <- diag(mapply( GIGrvg::rgig, n = 1, lambda = - (t_max - 1)*0.5, chi = sse_2,
                       psi = 1/prior$sigma_S0 ) , nrow = K)



  # # Invgamma conjugate prior
  # sigma_post_a <- 1 + rep(t_max,K) # prior of sigma_h Gamma(1,0.0001)
  # sigma_post_b <- 0.01 + sse_2 # prior of sigma_h
  #
  #   for (i in c(1:K)){
  #     sigma_h[i,i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
  #   }

  aux <- list(sigma_h = sigma_h,
              h0 = h0,
              Sigtdraw = h,
              sigt = exp(0.5*h))
  return(aux)
}

#' @export
sample_h_mod <- function(ytilde, sigma_h = 0.0001*diag(K), h0_mean = rep(0,K),
                         h = matrix(0, nrow = t_max, ncol = K), K, t_max, prior){
  tmp <- getmix()
  q <- tmp$q
  m_mean <- tmp$m
  u2 <- tmp$u2

    h0 <- rep(0,K)
    h0mean <- rep(0,K)
    h0var <- rep(0,K)

    cond_var_sigma <- rep(0,K)
    cond_mean_sigma <- rep(0,K)
    Zs <- matrix(1,t_max,1) %x% diag(1)
    for (i in c(1:K)){
      sigma_prmean <- h0_mean[i] # mean h_0
      sigma_prvar <- matrix(4)   # variance h_0
      aux <- sigmahelper4(t(ytilde[ i,, drop =FALSE]^2), q, m_mean, u2, h[ i,, drop =FALSE], Zs, matrix(sigma_h[i,i]), sigma_prmean, sigma_prvar)
      h[i,] <- aux$Sigtdraw
      h0[i] <- as.numeric(aux$h0)
      h0mean[i] <- as.numeric(aux$h0mean)
      h0var[i] <- as.numeric(aux$h0var)

    }

  aux <- list(h0 = h0,
              Sigtdraw = h,
              sigt = exp(0.5*h),
              h0mean = h0mean,
              h0var = h0var)
  return(aux)
}
##########################################################################
# Plot functions  #
##########################################################################

#' @export
plot.fatBVARSV <- function(fatBVARSVobj, element = NULL){
  if (is.null(element)) {
    plot(fatBVARSVobj$mcmc)
  } else {
    plot(get_post.fatBVARSV(fatBVARSVobj, element))
  }
}

##########################################################################
# Get posterior functions  #
##########################################################################
#' Get Posterior samples of BVAR model
#'
#' This function returns Posterior samples of BVAR-SV-fatTail model.
#' @param obj The Chain/mcmc obtained by fatBVARSVobj$mcmc
#' @param element The name of parameters.
#' @return The Posterior samples of BVAR-SV-fatTail model in mcmc object
#' @export
#' @examples
#' \dontrun{
#' B_samples <- get_post(Chain, element = "B")
#' }
get_post <- function(obj, element = NULL, ...) {
  UseMethod("get_post", obj)
}

#' @export
get_post.fatBVARSV <- function(obj, element = NULL){
  if (is.null(element)) {
    return(obj$mcmc)
  } else {
    mcmc_name <- substr(colnames(obj$mcmc), start = 1, stop = nchar(element))
    mcmc_id = (mcmc_name == element)
    return(obj$mcmc[,mcmc_id, drop = FALSE])
  }
}

#' @export
get_post.numeric <- function(obj, element = NULL){
  if (is.null(element)) {
    return(obj)
  } else {
    mcmc_name <- substr(names(obj), start = 1, stop = nchar(element))
    mcmc_id = (mcmc_name == element)
    return(obj[mcmc_id])
  }
}

#' @export
get_post.matrix <- function(obj, element = NULL){
  if (is.null(element)) {
    return(obj)
  } else {
    mcmc_name <- substr(colnames(obj), start = 1, stop = nchar(element))
    mcmc_id = (mcmc_name == element)
    return(as.matrix(obj[,mcmc_id], nrow = nrow(obj)))
  }
}

#' @export
get_post.mcmc <- function(obj, element = NULL){
  if (is.null(element)) {
    return(obj)
  } else {
    mcmc_name <- substr(colnames(obj), start = 1, stop = nchar(element))
    mcmc_id = (mcmc_name == element)
    return(obj[,mcmc_id, drop = FALSE])
  }
}


PF_GaussianSV <- function(ytilde, h0, sigma_h, noParticles = 1000) {

  t_len <- length(ytilde) - 1

  particles <- matrix(0, nrow = noParticles, ncol = t_len + 1)
  ancestorIndices <- matrix(0, nrow = noParticles, ncol = t_len + 1)
  weights <- matrix(1, nrow = noParticles, ncol = t_len + 1)
  normalisedWeights <- matrix(0, nrow = noParticles, ncol = t_len + 1)
  xHatFiltered <- matrix(0, nrow = t_len, ncol = 1)
  logLikelihood <- 0

  ancestorIndices[, 1] <- 1:noParticles
  normalisedWeights[, 1] = 1 / noParticles

  # Generate initial state
  particles[, 1] <- h0 + rnorm(noParticles, mean = 0, sd = sigma_h^0.5)

  for (t in 2:(t_len + 1)) {
    # Resample ( multinomial )
    newAncestors <- sample(noParticles, replace = TRUE, prob = normalisedWeights[, t - 1])
    ancestorIndices[, 1:(t - 1)] <- ancestorIndices[newAncestors, 1:(t - 1)]
    ancestorIndices[, t] <- newAncestors

    # Propagate
    part1 <- particles[newAncestors, t - 1]
    particles[, t] <- part1 + rnorm(noParticles, mean = 0, sd = sigma_h^0.5)

    # Compute weights
    yhatMean <- 0
    yhatVariance <- exp(particles[, t] / 2)
    weights[, t] <- dnorm(ytilde[t - 1], yhatMean, yhatVariance, log = TRUE)

    maxWeight <- max(weights[, t])
    weights[, t] <- exp(weights[, t] - maxWeight)

    sumWeights <- sum(weights[, t])
    normalisedWeights[, t] <- weights[, t] / sumWeights

    # Estimate the log-likelihood
    logLikelihood <- logLikelihood + maxWeight + log(sumWeights) - log(noParticles)

  }

  # Sample the state estimate using the weights at t=t_len
  ancestorIndex  <- sample(noParticles, 1, prob = normalisedWeights[, t_len])
  xHatFiltered <- particles[cbind(ancestorIndices[ancestorIndex, ], 1:(t_len + 1))]

  list(xHatFiltered = xHatFiltered, logLikelihood = logLikelihood)
}

PF_uniStudentSV <- function(ytilde, h0, sigma_h, nu, noParticles = 1000) {

  t_len <- length(ytilde) - 1

  particles <- matrix(0, nrow = noParticles, ncol = t_len + 1)
  ancestorIndices <- matrix(0, nrow = noParticles, ncol = t_len + 1)
  weights <- matrix(1, nrow = noParticles, ncol = t_len + 1)
  normalisedWeights <- matrix(0, nrow = noParticles, ncol = t_len + 1)
  xHatFiltered <- matrix(0, nrow = t_len, ncol = 1)
  logLikelihood <- 0

  ancestorIndices[, 1] <- 1:noParticles
  normalisedWeights[, 1] = 1 / noParticles

  # Generate initial state
  particles[, 1] <- h0 + rnorm(noParticles, mean = 0, sd = sigma_h^0.5)

  for (t in 2:(t_len + 1)) {
    # Resample ( multinomial )
    newAncestors <- sample(noParticles, replace = TRUE, prob = normalisedWeights[, t - 1])
    ancestorIndices[, 1:(t - 1)] <- ancestorIndices[newAncestors, 1:(t - 1)]
    ancestorIndices[, t] <- newAncestors

    # Propagate
    part1 <- particles[newAncestors, t - 1]
    particles[, t] <- part1 + rnorm(noParticles, mean = 0, sd = sigma_h^0.5)

    # Compute weights
    yhatMean <- 0
    yhatVariance <- exp(particles[, t] / 2)
    #weights[, t] <- LaplacesDemon::dst(ytilde[t - 1], mu = yhatMean, sigma = yhatVariance, nu = nu, log = TRUE)
    weights[, t] <- dt(ytilde[t - 1]/yhatVariance, df = nu, log = T) - log(yhatVariance)

    maxWeight <- max(weights[, t])
    weights[, t] <- exp(weights[, t] - maxWeight)

    sumWeights <- sum(weights[, t])
    normalisedWeights[, t] <- weights[, t] / sumWeights

    # Estimate the log-likelihood
    logLikelihood <- logLikelihood + maxWeight + log(sumWeights) - log(noParticles)

  }

  # Sample the state estimate using the weights at t=t_len
  ancestorIndex  <- sample(noParticles, 1, prob = normalisedWeights[, t_len])
  xHatFiltered <- particles[cbind(ancestorIndices[ancestorIndex, ], 1:(t_len + 1))]

  list(xHatFiltered = xHatFiltered, logLikelihood = logLikelihood)
}

PF_uniSkewSV <- function(ytilde, h0, sigma_h, nu, gamma, noParticles = 1000) {

  t_len <- length(ytilde) - 1

  particles <- matrix(0, nrow = noParticles, ncol = t_len + 1)
  ancestorIndices <- matrix(0, nrow = noParticles, ncol = t_len + 1)
  weights <- matrix(1, nrow = noParticles, ncol = t_len + 1)
  normalisedWeights <- matrix(0, nrow = noParticles, ncol = t_len + 1)
  xHatFiltered <- matrix(0, nrow = t_len, ncol = 1)
  logLikelihood <- 0

  ancestorIndices[, 1] <- 1:noParticles
  normalisedWeights[, 1] = 1 / noParticles

  # Generate initial state
  particles[, 1] <- h0 +  rnorm(noParticles, mean = 0, sd = sigma_h^0.5)
  for (t in 2:(t_len + 1)) {
    # Resample ( multinomial )
    newAncestors <- sample(noParticles, replace = TRUE, prob = normalisedWeights[, t - 1])
    ancestorIndices[, 1:(t - 1)] <- ancestorIndices[newAncestors, 1:(t - 1)]
    ancestorIndices[, t] <- newAncestors

    # Propagate
    part1 <- particles[newAncestors, t - 1]
    particles[, t] <- part1 + rnorm(noParticles, mean = 0, sd = sigma_h^0.5)

    # Compute weights
    yhatMean <- 0
    yhatVariance <- exp(particles[, t] / 2)
    weights[, t] <- mapply(FUN = SkewHyperbolic::dskewhyp, x = ytilde[t - 1],
                           mu = yhatMean, delta = sqrt(nu) * yhatVariance, beta = gamma, nu = nu, log = TRUE)

    #SkewHyperbolic::dskewhyp(x = ytilde[t - 1], mu = yhatMean, delta = sqrt(nu) * yhatVariance, beta = gamma, nu = nu, log = TRUE)

    maxWeight <- max(weights[, t])
    weights[, t] <- exp(weights[, t] - maxWeight)

    sumWeights <- sum(weights[, t])
    normalisedWeights[, t] <- weights[, t] / sumWeights

    # Estimate the log-likelihood
    logLikelihood <- logLikelihood + maxWeight + log(sumWeights) - log(noParticles)

  }

  # Sample the state estimate using the weights at t=t_len
  ancestorIndex  <- sample(noParticles, 1, prob = normalisedWeights[, t_len])
  xHatFiltered <- particles[cbind(ancestorIndices[ancestorIndex, ], 1:(t_len + 1))]

  list(xHatFiltered = xHatFiltered, logLikelihood = logLikelihood)
}

PF_StudentSV <- function(ytilde, h0, sigma_h, nu, noParticles = 1000) {

  t_len <- length(ytilde) - 1

  particles <- matrix(0, nrow = noParticles, ncol = t_len + 1)
  ancestorIndices <- matrix(0, nrow = noParticles, ncol = t_len + 1)
  weights <- matrix(1, nrow = noParticles, ncol = t_len + 1)
  normalisedWeights <- matrix(0, nrow = noParticles, ncol = t_len + 1)
  xHatFiltered <- matrix(0, nrow = t_len, ncol = 1)
  logLikelihood <- 0

  ancestorIndices[, 1] <- 1:noParticles
  normalisedWeights[, 1] = 1 / noParticles

  # Generate initial state
  particles[, 1] <- h0 + rnorm(noParticles, mean = 0, sd = sigma_h^0.5)

  for (t in 2:(t_len + 1)) {
    # Resample ( multinomial )
    newAncestors <- sample(noParticles, replace = TRUE, prob = normalisedWeights[, t - 1])
    ancestorIndices[, 1:(t - 1)] <- ancestorIndices[newAncestors, 1:(t - 1)]
    ancestorIndices[, t] <- newAncestors

    # Propagate
    part1 <- particles[newAncestors, t - 1]
    particles[, t] <- part1 + rnorm(noParticles, mean = 0, sd = sigma_h^0.5)

    # Compute weights
    yhatMean <- 0
    yhatVariance <- exp(particles[, t] / 2)
    #weights[, t] <- LaplacesDemon::dst(ytilde[t - 1], mu = yhatMean, sigma = yhatVariance, nu = nu, log = TRUE)
    weights[, t] <- dt(ytilde[t - 1]/yhatVariance, df = nu, log = T) - log(yhatVariance)

    maxWeight <- max(weights[, t])
    weights[, t] <- exp(weights[, t] - maxWeight)

    sumWeights <- sum(weights[, t])
    normalisedWeights[, t] <- weights[, t] / sumWeights

    # Estimate the log-likelihood
    logLikelihood <- logLikelihood + maxWeight + log(sumWeights) - log(noParticles)

  }

  # Sample the state estimate using the weights at t=t_len
  ancestorIndex  <- sample(noParticles, 1, prob = normalisedWeights[, t_len])
  xHatFiltered <- particles[cbind(ancestorIndices[ancestorIndex, ], 1:(t_len + 1))]

  list(xHatFiltered = xHatFiltered, logLikelihood = logLikelihood)
}

#' @export
logmeanexp <- function(sum_log){
  if (is(sum_log, "numeric")){
    max_sumlog <- max(sum_log)
    out <- log( mean(exp(sum_log - max_sumlog))) + max_sumlog
  } else {
    max_sumlog = apply(sum_log, MARGIN = 2 , FUN = max)

    out = log( apply(exp(sum_log-reprow(max_sumlog,nrow(sum_log))), MARGIN = 2, FUN = mean )) + max_sumlog

  }
  return(out)
}
