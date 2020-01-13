#' Marginal log likelihood of BVAR model
#'
#' This function returns a marginal log likelihood density of BVAR-SV-fatTail model.
#' @param Chain The fatBVARSV object from command BVAR.
#' @return The marginal log likelihood density of of BVAR model
#' @export
#' @examples
#' \dontrun{
#' marginal_dens1 <- marginalLL(Chain1)
#' }
#'
#' @export
marginalLL <- function(Chain){
  K <- Chain$K
  p <- Chain$p

  m <- K + K*K*p

  dist <- Chain$dist
  SV <- Chain$prior$SV

  y <- Chain$y
  t_max <- nrow(y)
  if (is.null(Chain$y0)){
    y0 <- matrix(0, ncol = K, nrow = p)
  } else {
    y0 <- Chain$y0
  }

  # y_combine <- rbind( tail(y0, p) , y)

  mcmc <- Chain$mcmc
  ndraws <- nrow(mcmc)

  B_mat <- get_post(mcmc, element = "B")
  A_mat <- get_post(mcmc, element = "a")
  Gamma_mat <- get_post(mcmc, element = "gamma")
  BAG_samples <- Normal_approx(cbind(B_mat,A_mat, Gamma_mat), ndraws = ndraws)
  B_gen <- BAG_samples$B_gen
  A_gen <- BAG_samples$A_gen
  Gamma_gen <- BAG_samples$Gamma_gen
  sum_log_prop <- BAG_samples$sum_log_prop

  Sigma_mat <- get_post(mcmc, element = "sigma") # No SV
  Sigma_H <- sqrt(get_post(mcmc, element = "sigma_h")) # SV
  H_mat <- get_post(mcmc, element = "h")
  Nu_mat <- as.matrix(get_post(mcmc, element = "nu"))
  colnames(Nu_mat) <- rep("nu", ncol(Nu_mat))
  W_mat <- get_post(mcmc, element = "w")

  sum_log <- rep(0, ndraws)


  if (SV){

  } else {

    Sigma_gen_list <- InvGamma_approx(Sigma_mat, ndraws = ndraws)
    Sigma_gen <- Sigma_gen_list$new_samples
    sum_log_prop <- sum_log_prop + Sigma_gen_list$sum_log_prop

    if (dist == "Gaussian") {
      param <- cbind(B_gen, A_gen, Sigma_gen)
      for (j in c(1:ndraws)){
        sum_log[j] <- log_posterior(param[j,], Chain)
      }
    } else {
      Nu_gen_list <- Gamma_approx(Nu_mat, ndraws = ndraws)
      Nu_gen <- Nu_gen_list$new_samples
      sum_log_prop <- sum_log_prop + Nu_gen_list$sum_log_prop

      W_gen_list <- InvGamma_approx(W_mat, ndraws = ndraws)
      W_gen <- W_gen_list$new_samples
      sum_log_prop <- sum_log_prop + W_gen_list$sum_log_prop

      param <- cbind(B_gen, A_gen, Gamma_gen, Sigma_gen, Nu_gen, W_gen)

      if (dist == "Student") {
        for (j in c(1:ndraws)){
          sum_log[j] <- log_posterior(param[j,], Chain)
        }
      }

      if (dist == "Hyper.Student") {
        for (j in c(1:ndraws)){
          sum_log[j] <- log_posterior(param[j,], Chain)
        }
      }

      if (dist == "multiStudent") {
        for (j in c(1:ndraws)){
          sum_log[j] <- log_posterior(param[j,], Chain)
        }
      }

      if (dist == "Hyper.multiStudent") {
        for (j in c(1:ndraws)){
          sum_log[j] <- log_posterior(param[j,], Chain)
        }
      }


    }


  }


  sum_log <- sum_log - sum_log_prop # log posterior - log proposal

  short_sumlog = matrix(sum_log, nrow = ndraws/20, ncol = 20)
  max_sumlog = apply(short_sumlog, MARGIN = 2 , FUN = max)
  bigml = log( apply(exp(short_sumlog-reprow(max_sumlog,ndraws/20)), MARGIN = 2, FUN = mean )) + max_sumlog
  ml = mean(bigml)
  mlstd = std(bigml)/sqrt(20)


  # return( list( LL = mean(sum_log),
  #               std = std(sum_log)/sqrt(ndraws)))
  return( list( LL = ml,
                std = mlstd))
}

#' @export
Normal_approx <- function(mcmc_sample, ndraws){
  mcmc_mean <- apply(mcmc_sample, 2, mean)
  mcmc_Sigma <- cov(mcmc_sample)
  nElements <- length(mcmc_mean)
  # mcmc_chol <- t(chol(cov(mcmc_sample)))
  new_samples <- mvnfast::rmvn(ndraws, mu = mcmc_mean, sigma = mcmc_Sigma)
  colnames(new_samples) <- colnames(mcmc_sample)
  nB = sum(substr(colnames(mcmc_sample),1,1) == "B")
  nA = sum(substr(colnames(mcmc_sample),1,1) == "a")
  nG = sum(substr(colnames(mcmc_sample),1,5) == "gamma")
  B_gen <- new_samples[,1:nB]
  A_gen <- new_samples[,(nB+1):(nB+nA)]
  Gamma_gen <- NULL
  if (nG > 0) {
    Gamma_gen <- new_samples[,(nB+nA+1):(nB+nA+nG)]
  }
  sum_log_prop <- mvnfast::dmvn(X = new_samples,
                                mu = as.vector(mcmc_mean), sigma = mcmc_Sigma, log = T)
  return(list(B_gen = B_gen,
              A_gen = A_gen,
              Gamma_gen = Gamma_gen,
              sum_log_prop = sum_log_prop))
}

#' @export
InvGamma_approx <- function(mcmc_sample, ndraws){
  nElements <- ncol(mcmc_sample)
  new_samples <- matrix(NA, ncol = nElements, nrow = ndraws, dimnames = list(c(), colnames(mcmc_sample)))
  Density_prop <-  matrix(NA, ncol = nElements, nrow = ndraws)
  shape_param <- rep(0, nElements)
  rate_param <- rep(0, nElements)
  for (i in c(1:nElements)){
    fit.gamma <- fitdistrplus::fitdist(1/as.numeric(mcmc_sample[,i]), distr = "gamma", method = "mle")
    shape_param[i] <- fit.gamma$estimate[1]
    rate_param[i] <- fit.gamma$estimate[2]
    new_samples[,i] <- rinvgamma(ndraws, shape = shape_param[i], rate_param[i])
    Density_prop[,i] <- invgamma::dinvgamma(new_samples[,i],
                                               shape = shape_param[i],
                                               rate = rate_param[i], log = T)
  }
  return(list(new_samples = new_samples,
              sum_log_prop = apply(Density_prop, 1, sum)))
}

#' @export
Gamma_approx <- function(mcmc_sample, ndraws){
  nElements <- ncol(mcmc_sample)
  new_samples <- matrix(NA, ncol = nElements, nrow = ndraws, dimnames = list(c(), colnames(mcmc_sample)))
  Density_prop <-  matrix(NA, ncol = nElements, nrow = ndraws)
  shape_param <- rep(0, nElements)
  rate_param <- rep(0, nElements)
  for (i in c(1:nElements)){
    fit.gamma <- fitdistrplus::fitdist(as.numeric(mcmc_sample[,i]), distr = "gamma", method = "mle")
    shape_param[i] <- fit.gamma$estimate[1]
    rate_param[i] <- fit.gamma$estimate[2]
    new_samples[,i] <- rgamma(ndraws, shape = shape_param[i], rate_param[i])
    Density_prop[,i] <- dgamma(new_samples[,i], shape = shape_param[i], rate = rate_param[i], log = T)
  }
  return(list(new_samples = new_samples,
              sum_log_prop = apply(Density_prop, 1, sum)))
}
