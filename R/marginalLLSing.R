#' Marginal log likelihood of BVAR model
#'
#' This function returns a marginal log likelihood density of BVAR-SV-fatTail model.
#' @param Chain The fatBVARSV object from command BVAR.
#' @return The marginal log likelihood density of of BVAR model
#' @export
#' @examples
#' \dontrun{
#' marginal_dens1 <- marginalLLSing(Chain1)
#' }
#'
#' @export
marginalLLSing <- function(Chain, ndraws = NULL){
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
  if (is.null(ndraws)) ndraws <- nrow(mcmc)

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
  H0_mat <- get_post(mcmc, element = "lh0")
  Nu_mat <- as.matrix(get_post(mcmc, element = "nu"))
  if (ncol(Nu_mat) == 1) colnames(Nu_mat) <- "nu"
  W_mat <- get_post(mcmc, element = "w")

  sum_log <- rep(0, ndraws)


  if (SV){
    h_mean <- matrix(apply(H_mat, MARGIN = 2, mean), nrow = K)
    H0_gen_list <- Normal_approx(H0_mat, ndraws = ndraws)
    H0_gen <- H0_gen_list$new_samples
    sum_log_prop <- sum_log_prop + H0_gen_list$sum_log_prop

    #################### prior ###########################

    if (dist == "Gaussian") {
      # param <- cbind(B_gen, A_gen, Sigma_gen)

      for (j in c(1:ndraws)){
        B <- matrix(B_gen[j,], nrow = K)
        A <- a0toA(A_gen[j,], K)
        sigma_h <- Sigma_gen[j,]
        h0 <- H0_gen[j,]
        ytilde <- as.vector(A %*% (yt - B %*% xt)) # from dim K * t_max to (Kt_max) * 1

        sum_log[j] <- int_h(ytilde = ytilde, h0 = h0, sigma_h = sigma_h, t_max = t_max, K = K) +
          mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
          mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
          sum(invgamma::dinvgamma(sigma_h, shape = 1* 0.5, rate = 0.0001 * 0.5, log = T)) +
          sum(dnorm(h0, mean = 0, sd = sqrt(10), log = T))
      }
    } else {
      Nu_gen_list <- Gamma_approx(Nu_mat, ndraws = ndraws)
      Nu_gen <- Nu_gen_list$new_samples
      sum_log_prop <- sum_log_prop + Nu_gen_list$sum_log_prop

      if (dist == "Student") {
        for (j in c(1:ndraws)){
          B <- matrix(B_gen[j,], nrow = K)
          A <- a0toA(A_gen[j,], K)
          sigma_h <- Sigma_gen[j,]
          h0 <- H0_gen[j,]
          nu <- Nu_gen[j]

          llw <- int_h_Student(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                               h_mean = h_mean, t_max = t_max, K = K, R = 10)

          sum_log[j] <- llw +
            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
            sum(invgamma::dinvgamma(sigma_h, shape = 1* 0.5, rate = 0.0001 * 0.5, log = T)) +
            sum(dnorm(h0, mean = 0, sd = sqrt(10), log = T)) +
            dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)

          }
      }

      if (dist == "Skew.Student") {
        for (j in c(1:ndraws)){
          B <- matrix(B_gen[j,], nrow = K)
          A <- a0toA(A_gen[j,], K)
          gamma <- Gamma_gen[j,]

          sigma_h <- Sigma_gen[j,]
          h0 <- H0_gen[j,]
          nu <- Nu_gen[j]


          llw <- int_h_Student(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                               gamma = gamma, h_mean = h_mean, t_max = t_max, K = K, R = 10)

          sum_log[j] <- llw +
            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
            sum(invgamma::dinvgamma(sigma_h, shape = 1* 0.5, rate = 0.0001 * 0.5, log = T)) +
            sum(dnorm(h0, mean = 0, sd = sqrt(10), log = T)) +
            dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T) +
            mvnfast::dmvn(X = gamma, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
        }
      }

      if (dist == "MT") {
        for (j in c(1:ndraws)){
          B <- matrix(B_gen[j,], nrow = K)
          A <- a0toA(A_gen[j,], K)
          sigma_h <- Sigma_gen[j,]
          h0 <- H0_gen[j,]
          nu <- Nu_gen[j,]


          llw <- int_h_MultiStudent(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                    h_mean = h_mean, t_max = t_max, K = K, R = 100)

          sum_log[j] <- llw +
            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
            sum(invgamma::dinvgamma(sigma_h, shape = 1* 0.5, rate = 0.0001 * 0.5, log = T)) +
            sum(dnorm(h0, mean = 0, sd = sqrt(10), log = T)) +
            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T))
        }
      }

      if (dist == "MST") {
        for (j in c(1:ndraws)){
          B <- matrix(B_gen[j,], nrow = K)
          A <- a0toA(A_gen[j,], K)
          gamma <- Gamma_gen[j,]

          sigma_h <- Sigma_gen[j,]
          h0 <- H0_gen[j,]
          nu <- Nu_gen[j,]


          llw <- int_h_MultiStudent(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                    gamma = gamma, h_mean = h_mean, t_max = t_max, K = K, R = 100)

          sum_log[j] <- llw +
            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
            sum(invgamma::dinvgamma(sigma_h, shape = 1* 0.5, rate = 0.0001 * 0.5, log = T)) +
            sum(dnorm(h0, mean = 0, sd = sqrt(10), log = T)) +
            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
            mvnfast::dmvn(X = gamma, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
        }
      }

    }



  } else {


    if (dist == "Gaussian") {
      param <- cbind(B_gen, A_gen, Sigma_gen)
      for (j in c(1:ndraws)){
        sum_log[j] <- log_posterior(param[j,], Chain)
      }
    } else {
      # Fat tail
      Nu_gen_list <- Gamma_approx(Nu_mat, ndraws = ndraws)
      Nu_gen <- Nu_gen_list$new_samples
      sum_log_prop <- sum_log_prop + Nu_gen_list$sum_log_prop

      if (dist == "Student") {
        for (j in c(1:ndraws)){
          B <- matrix(B_gen[j,], nrow = K)
          A <- a0toA(A_gen[j,], K)
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

      if (dist == "Skew.Student") {
        for (j in c(1:ndraws)){
          B <- matrix(B_gen[j,], nrow = K)
          A <- a0toA(A_gen[j,], K)
          gamma <- Gamma_gen[j,]

          sigma <- Sigma_gen[j,]
          nu <- Nu_gen[j]

          llw <- int_w_Student(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu, gamma = gamma,
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

      if (dist == "MT") {
        for (j in c(1:ndraws)){
          B <- matrix(B_gen[j,], nrow = K)
          A <- a0toA(A_gen[j,], K)
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

      if (dist == "MST") {
        for (j in c(1:ndraws)){
          B <- matrix(B_gen[j,], nrow = K)
          A <- a0toA(A_gen[j,], K)
          gamma <- Gamma_gen[j,]
          sigma <- Sigma_gen[j,]
          nu <- Nu_gen[j,]

          llw <- int_w_MultiStudent(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu, gamma = gamma,
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
  mlstd = sd(bigml)/sqrt(20)

  # return( list( LL = mean(sum_log),
  #               std = sd(sum_log)/sqrt(ndraws)))
  return( list( LL = ml,
                std = mlstd,
                sum_log = sum_log))
}
