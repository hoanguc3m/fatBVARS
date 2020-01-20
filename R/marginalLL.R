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
  prior <- Chain$prior
  SV <- Chain$prior$SV

  y <- Chain$y
  yt <- t(y)
  t_max <- nrow(y)
  if (is.null(Chain$y0)){
    y0 <- matrix(0, ncol = K, nrow = p)
  } else {
    y0 <- Chain$y0
  }
  xt <- makeRegressor(y, y0, t_max, K, p)

  # y_combine <- rbind( tail(y0, p) , y)

  mcmc <- Chain$mcmc
  ndraws <- nrow(mcmc)

  B_mat <- get_post(mcmc, element = "B")
  A_mat <- get_post(mcmc, element = "a")
  Gamma_mat <- get_post(mcmc, element = "gamma")
  BAG_samples <- Normal_approx(cbind(B_mat,A_mat, Gamma_mat), ndraws = ndraws)
  B_gen <- get_post(BAG_samples$new_samples, element = "B")
  A_gen <- get_post(BAG_samples$new_samples, element = "a")
  Gamma_gen <- get_post(BAG_samples$new_samples, element = "gamma")
  sum_log_prop <- BAG_samples$sum_log_prop

  Sigma_mat <- get_post(mcmc, element = "sigma") # No SV is sigma / SV is sigma^2
  # Sigma_H <- sqrt(get_post(mcmc, element = "sigma_h")) # SV
  if (SV) {
    Sigma_gen_list <- InvGamma_approx(Sigma_mat, ndraws = ndraws)
    Sigma_gen <- Sigma_gen_list$new_samples
  } else {
    Sigma_gen_list <- InvGamma_approx(Sigma_mat^2, ndraws = ndraws)
    Sigma_gen <- sqrt(Sigma_gen_list$new_samples)
  }
  sum_log_prop <- sum_log_prop + Sigma_gen_list$sum_log_prop

  H_mat <- get_post(mcmc, element = "h")
  Nu_mat <- as.matrix(get_post(mcmc, element = "nu"))
  if (ncol(Nu_mat) == 1) colnames(Nu_mat) <- "nu"
  W_mat <- get_post(mcmc, element = "w")

  sum_log <- rep(0, ndraws)


  if (SV){
    H_mat <- H_mat[,1:K]
    H0_gen_list <- Normal_approx(H_mat, ndraws = ndraws)
    H0_gen <- H0_gen_list$new_samples
    sum_log_prop <- sum_log_prop + H0_gen_list$sum_log_prop

    #################### prior ###########################

    if (dist == "Gaussian") {
      # param <- cbind(B_gen, A_gen, Sigma_gen)

      for (j in c(1:ndraws)){
        B <- matrix(B_gen[j,], nrow = K)
        A <- matrix(0, nrow = K, ncol = K)
        A[upper.tri(A)] <- A_gen[j,]
        A <- t(A)
        diag(A) <- 1

        sigma_h <- Sigma_gen[j,]
        h0 <- H0_gen[j,]

        ytilde <- as.vector(A %*% (yt - B %*% xt))

        sum_log[j] <- int_h(ytilde = ytilde, h0 = h0, sigma_h = sigma_h, t_max = t_max, K = K) +
          mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
          mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
          sum(invgamma::dinvgamma(sigma_h, shape = 1* 0.5, rate = 0.0001 * 0.5, log = T)) +
          sum(dnorm(h0, mean = 0, sd = 1, log = T))
      }
    }



  } else {


    if (dist == "Gaussian") {
      param <- cbind(B_gen, A_gen, Sigma_gen)
      for (j in c(1:ndraws)){
        sum_log[j] <- log_posterior(param[j,], Chain)
      }
    } else {
      Nu_gen_list <- Gamma_approx(Nu_mat, ndraws = ndraws)
      Nu_gen <- Nu_gen_list$new_samples
      sum_log_prop <- sum_log_prop + Nu_gen_list$sum_log_prop

      # W_gen_list <- InvGamma_approx(W_mat, ndraws = ndraws)
      # W_gen <- W_gen_list$new_samples
      # sum_log_prop <- sum_log_prop + W_gen_list$sum_log_prop
      #
      # param <- cbind(B_gen, A_gen, Gamma_gen, Sigma_gen, Nu_gen, W_gen)

      if (dist == "Student") {
        for (j in c(1:ndraws)){
          B <- matrix(B_gen[j,], nrow = K)
          A <- matrix(0, nrow = K, ncol = K)
          A[upper.tri(A)] <- A_gen[j,]
          A <- t(A)
          diag(A) <- 1

          sigma <- Sigma_gen[j,]
          nu <- Nu_gen[j]

          llw <- int_w_Student(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu,
                               t_max = t_max, K = K, R = 100)
          # int_w_Student2(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu, t_max = t_max, K = K, R = 100)


          sum_log[j] <- llw +
            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
            sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5,
                                    rate = prior$sigma_S0 * 0.5, log = T)) +
            dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)
        }
      }

      if (dist == "Hyper.Student") {
        for (j in c(1:ndraws)){
          B <- matrix(B_gen[j,], nrow = K)
          A <- matrix(0, nrow = K, ncol = K)
          A[upper.tri(A)] <- A_gen[j,]
          A <- t(A)
          diag(A) <- 1
          gamma <- Gamma_gen[j,]

          sigma <- Sigma_gen[j,]
          nu <- Nu_gen[j]

          llw <- int_w_HyperStudent(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu, gamma = gamma,
                                    t_max = t_max, K = K, R = 100)

          sum_log[j] <- llw +
            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
            sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5,
                                    rate = prior$sigma_S0 * 0.5, log = T)) +
            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
            mvnfast::dmvn(X = gamma, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)

          # sum_log[j] <- log_posterior(param[j,], Chain)
        }
      }

      if (dist == "multiStudent") {
        for (j in c(1:ndraws)){
          B <- matrix(B_gen[j,], nrow = K)
          A <- matrix(0, nrow = K, ncol = K)
          A[upper.tri(A)] <- A_gen[j,]
          A <- t(A)
          diag(A) <- 1

          sigma <- Sigma_gen[j,]
          nu <- Nu_gen[j,]

          llw <- int_w_MultiStudent(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu,
                                    t_max = t_max, K = K, R = 100)

          sum_log[j] <- llw +
            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
            sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5,
                                    rate = prior$sigma_S0 * 0.5, log = T)) +
            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T))


          # sum_log[j] <- log_posterior(param[j,], Chain)
        }
      }

      if (dist == "Hyper.multiStudent") {
        for (j in c(1:ndraws)){
          B <- matrix(B_gen[j,], nrow = K)
          A <- matrix(0, nrow = K, ncol = K)
          A[upper.tri(A)] <- A_gen[j,]
          A <- t(A)
          diag(A) <- 1
          gamma <- Gamma_gen[j,]

          sigma <- Sigma_gen[j,]
          nu <- Nu_gen[j,]

          llw <- int_w_HyperMultiStudent(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu, gamma = gamma,
                                         t_max = t_max, K = K, R = 100)

          sum_log[j] <- llw +
            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
            sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5,
                                    rate = prior$sigma_S0 * 0.5, log = T)) +
            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
            mvnfast::dmvn(X = gamma, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)


          # sum_log[j] <- log_posterior(param[j,], Chain)
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

int_w_Student <- function(y, xt, A, B, sigma, nu, t_max, K, R = 100){
  u <- (y - t(xt) %*% t(B))
  #ytilde <- (A %*% (yt - B %*% xt))
  ytilde <- (A %*% t(u))
  shape <- nu*0.5 + K*0.5
  rate <- nu*0.5 + 0.5 * apply(ytilde^2, 2, sum)/sigma^2

  store_llw <- matrix(NA, nrow = t_max, ncol = R)
  for (i in c(1:R)){
    w <- mapply(FUN = rinvgamma, n = 1, shape = shape, rate = rate)
    dens_w <- mapply(FUN = dinvgamma,x = w, shape = shape, rate = rate, log = T)
    Prior_dens_w <- mapply(FUN = dinvgamma,x = w, shape = nu * 0.5, rate = nu * 0.5, log = T)
    store_llw[,i] <- mvnfast::dmvn(X = u/repcol(sqrt(w),K),
                                      mu = rep(0,K),
                                      sigma = solve(A) %*% diag(sigma), log = T, isChol = T) -
                    0.5 * K * log(w) + Prior_dens_w - dens_w
  }
  maxllike = apply(store_llw, MARGIN = 1, max)
  llw = log( apply(exp(store_llw-maxllike), MARGIN = 1, FUN = mean) ) + maxllike
  return(sum(llw))
}

int_w_Student2 <- function(y, xt, A, B, sigma, nu, t_max, K, R = 100){
  u <- (y - t(xt) %*% t(B))
  allw <- mvnfast::dmvt(X = u,
                                 mu = rep(0,K),
                                 sigma = solve(A) %*% diag(sigma), df = nu, log = T, isChol = T)

  return(sum(allw))
}

int_w_HyperStudent <- function(y, xt, A, B, sigma, nu, gamma, t_max, K, R = 100){
  u <- (y - t(xt) %*% t(B))
  #ytilde <- (A %*% (yt - B %*% xt))
  ytilde <- (A %*% t(u))
  shape <- nu*0.5 + K*0.5
  rate <- nu*0.5 + 0.5 * apply(ytilde^2, 2, sum)/sigma^2

  store_llw <- matrix(NA, nrow = t_max, ncol = R)
  for (i in c(1:R)){
    w <- matrix(mapply(FUN = rinvgamma, n = 1, shape = shape, rate = rate), ncol = 1)
    dens_w <- mapply(FUN = dinvgamma,x = w, shape = shape, rate = rate, log = T)
    Prior_dens_w <- mapply(FUN = dinvgamma,x = w, shape = nu * 0.5, rate = nu * 0.5, log = T)
    store_llw[,i] <- mvnfast::dmvn(X = (u - repcol(w,K) * reprow(gamma, t_max))/repcol(sqrt(w),K),
                                   mu = rep(0,K),
                                   sigma = solve(A) %*% diag(sigma), log = T, isChol = T) -
      0.5 * K * log(w) + Prior_dens_w - dens_w # proposal
  }
  maxllike = apply(store_llw, MARGIN = 1, max)
  llw = log( apply(exp(store_llw-maxllike), MARGIN = 1, FUN = mean) ) + maxllike
  return(sum(llw))
}


int_w_MultiStudent <- function(y, xt, A, B, sigma, nu, t_max, K, R = 100){
  # u <- (y - t(xt) %*% t(B))
  # u_proposal <- (A %*% t(u))

  u <- (t(y) - B %*%xt)
  u_proposal <- A %*% (u)

  a_target <- (nu*0.5 + 1*0.5)
  b_target <- (nu*0.5 + 0.5 * u_proposal^2 / sigma^2)


  store_llw <- matrix(NA, nrow = t_max, ncol = R)
  for (i in c(1:R)){
    w <- matrix(mapply(FUN = rinvgamma, n = 1, shape = a_target, rate = b_target), nrow = K)
    dens_w <- matrix(mapply(FUN = dinvgamma, x = w, shape = a_target, rate = b_target, log = T), nrow = K)
    Prior_dens_w <- matrix(mapply(FUN = dinvgamma, x = w, shape = nu*0.5, rate = nu*0.5, log = T), nrow = K)
    store_llw[,i] <- mvnfast::dmvn(X = t(u/sqrt(w)),
                                   mu = rep(0,K),
                                   sigma = solve(A) %*% diag(sigma), log = T, isChol = T) -
      0.5 * apply(log(w), MARGIN = 2, FUN = sum) + apply(Prior_dens_w, MARGIN = 2, FUN = sum) - apply(dens_w, MARGIN = 2, FUN = sum)
  }
  maxllike = apply(store_llw, MARGIN = 1, max)
  llw = log( apply(exp(store_llw-maxllike), MARGIN = 1, FUN = mean) ) + maxllike
  return(sum(llw))
}

int_w_HyperMultiStudent <- function(y, xt, A, B, sigma, nu, gamma, t_max, K, R = 100){
  # u <- (y - t(xt) %*% t(B))
  # u_proposal <- (A %*% t(u))

  u <- (t(y) - B %*%xt)
  u_proposal <- A %*% (u)

  a_target <- (nu*0.5 + 1*0.5)
  b_target <- (nu*0.5 + 0.5 * u_proposal^2 / sigma^2)


  store_llw <- matrix(NA, nrow = t_max, ncol = R)
  for (i in c(1:R)){
    w <- matrix(mapply(FUN = rinvgamma, n = 1, shape = a_target, rate = b_target), nrow = K)
    dens_w <- matrix(mapply(FUN = dinvgamma, x = w, shape = a_target, rate = b_target, log = T), nrow = K)
    Prior_dens_w <- matrix(mapply(FUN = dinvgamma, x = w, shape = nu*0.5, rate = nu*0.5, log = T), nrow = K)
    store_llw[,i] <- mvnfast::dmvn(X = t( (u - w * repcol(gamma, t_max) )  /sqrt(w)),
                                   mu = rep(0,K),
                                   sigma = solve(A) %*% diag(sigma), log = T, isChol = T) -
      0.5 * apply(log(w), MARGIN = 2, FUN = sum) + apply(Prior_dens_w, MARGIN = 2, FUN = sum) - apply(dens_w, MARGIN = 2, FUN = sum)
  }
  maxllike = apply(store_llw, MARGIN = 1, max)
  llw = log( apply(exp(store_llw-maxllike), MARGIN = 1, FUN = mean) ) + maxllike
  return(sum(llw))
}

# This function evaluates the integrated likelihood of the VAR-SV model in
# Chan and Eisenstat (2018)
#' @export
int_h <- function(ytilde, h0, sigma_h, t_max, K){
    s2 = ytilde^2
    max_loop = 100
    Hh = sparseMatrix(i = 1:(t_max*K),
                      j = 1:(t_max*K),
                      x = rep(1,t_max*K)) -
      sparseMatrix( i = (K+1):(t_max*K),
                    j = 1:((t_max-1)*K),
                    x = rep(1,(t_max-1)*K),
                    dims =  c(t_max*K, t_max*K))
    SH = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = rep(1./sigma_h, t_max))
    HinvSH_h = Matrix::t(Hh) %*% SH %*% Hh
    alph = Matrix::solve(Hh, sparseMatrix(i = 1:K, j = rep(1,K), x = h0, dims = c(t_max*K,1)))

    e_h = 1
    ht = rep(h0,t_max)
    count = 0
    while ( e_h> .01 & count < max_loop){
      einvhts2 = exp(-ht)*s2
      gh = - HinvSH_h %*% (ht-alph) - 0.5 * (1-einvhts2)
      Gh = - HinvSH_h -.5*sparseMatrix(i = 1:(t_max*K),j = 1:(t_max*K), x = einvhts2)
      newht = ht - Matrix::solve(Gh,gh)
      e_h = max(abs(newht-ht));
      ht = newht;
      count = count + 1;
    }
    if (count == max_loop){
      ht = rep(h0,t_max)
      einvhts2 = exp(-ht)*s2
      Gh = - HinvSH_h -.5*sparseMatrix(i = 1:(t_max*K),j = 1:(t_max*K), x = einvhts2)
    }

    Kh = -Gh
    CKh = Matrix::t(Matrix::chol(Kh))

    c_pri = -t_max*K*0.5*log(2*pi) -.5*t_max*sum(log(sigma_h))
    c_IS = -t_max*K*0.5*log(2*pi) + sum(log(Matrix::diag(CKh)))

    R = 10
    store_llike = rep(0,R)
    for (i in c(1:R)){
      hc = ht + Matrix::solve(Matrix::t(CKh), rnorm(t_max*K))
      llike = -t_max*K*0.5*log(2*pi) - 0.5*sum(hc) -
        0.5*t(ytilde) %*% sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = exp(-hc)) %*% ytilde
      store_llike[i] = as.numeric(llike + (c_pri -.5*Matrix::t(hc-alph)%*%HinvSH_h%*%(hc-alph)) -
                                    (c_IS -.5*Matrix::t(hc-ht)%*%Kh%*%(hc-ht)))

    }

  maxllike = max(store_llike)
  llk = log(mean(exp(store_llike-maxllike))) + maxllike
  return(llk)
}

#' @export
Normal_approx <- function(mcmc_sample, ndraws){
  mcmc_mean <- apply(mcmc_sample, 2, mean)
  mcmc_Sigma <- cov(mcmc_sample)
  nElements <- length(mcmc_mean)
  # mcmc_chol <- t(chol(cov(mcmc_sample)))
  new_samples <- mvnfast::rmvn(ndraws, mu = mcmc_mean, sigma = mcmc_Sigma)
  colnames(new_samples) <- colnames(mcmc_sample)
  # nB = sum(substr(colnames(mcmc_sample),1,1) == "B")
  # nA = sum(substr(colnames(mcmc_sample),1,1) == "a")
  # nG = sum(substr(colnames(mcmc_sample),1,5) == "gamma")
  # B_gen <- new_samples[,1:nB]
  # A_gen <- new_samples[,(nB+1):(nB+nA)]
  # Gamma_gen <- NULL
  # if (nG > 0) {
  #   Gamma_gen <- new_samples[,(nB+nA+1):(nB+nA+nG)]
  # }
  sum_log_prop <- mvnfast::dmvn(X = new_samples,
                                mu = as.vector(mcmc_mean), sigma = mcmc_Sigma, log = T)
  return(list(new_samples = new_samples,
              sum_log_prop = sum_log_prop))
}

#' @export
H0_approx <- function(mcmc_sample, ndraws){
  mcmc_mean <- apply(mcmc_sample, 2, mean)
  mcmc_Sigma <- cov(mcmc_sample)
  nElements <- length(mcmc_mean)
  # mcmc_chol <- t(chol(cov(mcmc_sample)))
  new_samples <- mvnfast::rmvn(ndraws, mu = mcmc_mean, sigma = mcmc_Sigma)
  colnames(new_samples) <- colnames(mcmc_sample)
  sum_log_prop <- mvnfast::dmvn(X = new_samples,
                                mu = as.vector(mcmc_mean), sigma = mcmc_Sigma, log = T)
  return(list(H0_gen = new_samples,
              sum_log_prop = sum_log_prop))
}

#' @export
InvGamma_approx <- function(mcmc_sample, ndraws){
  nElements <- ncol(mcmc_sample)
  new_samples <- matrix(NA, ncol = nElements, nrow = ndraws, dimnames = list(c(), colnames(mcmc_sample)))
  Density_prop <-  matrix(NA, ncol = nElements, nrow = ndraws)
  shape_param <- rep(0, nElements)
  rate_param <- rep(0, nElements)
  mcmc_mean <- apply(mcmc_sample, 2, mean)
  mcmc_sd <- apply(mcmc_sample, 2, sd)
  if (mean(mcmc_mean) < 0.01)  mcmc_sample <- mcmc_sample * 100

  for (i in c(1:nElements)){
    if (mcmc_sd[i] > 0){
      fit.gamma <- fitdistrplus::fitdist(1/as.numeric(mcmc_sample[,i]), distr = "gamma", method = "mle")
      shape_param[i] <- fit.gamma$estimate[1]
      rate_param[i] <- fit.gamma$estimate[2]
      new_samples[,i] <- rinvgamma(ndraws, shape = shape_param[i], rate_param[i])
      Density_prop[,i] <- invgamma::dinvgamma(new_samples[,i],
                                              shape = shape_param[i],
                                              rate = rate_param[i], log = T)
    } else {
      new_samples[,i] <- as.numeric(mcmc_sample[,i])
      Density_prop[,i] <- 0
    }

  }

  if (mean(mcmc_mean) < 0.01)  new_samples <- new_samples / 100
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

