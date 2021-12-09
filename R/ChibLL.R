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
ChibLLP <- function(Chain, ndraws = NULL, numCores = NULL){
    RhpcBLASctl::blas_set_num_threads(1)

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

    B_gen <- apply(B_mat, MARGIN = 2, FUN = mean)
    A_gen <- apply(A_mat, MARGIN = 2, FUN = mean)
    Gamma_gen <- apply(Gamma_mat, MARGIN = 2, FUN = mean)


    Sigma_mat <- get_post(mcmc, element = "sigma") # No SV is sigma / SV is sigma^2
    Sigma_gen <- apply(Sigma_mat, MARGIN = 2, FUN = mean)

    H_mat <- get_post(mcmc, element = "h")
    H0_mat <- get_post(mcmc, element = "lh0")
    Nu_mat <- as.matrix(get_post(mcmc, element = "nu"))

    W_mat <- get_post(mcmc, element = "w")

    H_gen <- apply(H_mat, MARGIN = 2, FUN = mean)
    H0_gen <- apply(H0_mat, MARGIN = 2, FUN = mean)
    Nu_gen <- apply(Nu_mat, MARGIN = 2, FUN = mean)


    sum_log <- rep(0, ndraws)
    lprior <- 0

    if (is.null(numCores)) numCores <- parallel::detectCores()*0.5

    if (SV){
      h_mean <- H_gen
      #################### End prior ###########################

      if (dist == "Gaussian") {
        lprior <- mvnfast::dmvn(X = B_gen, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                  mvnfast::dmvn(X = A_gen, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                  # sum(dnorm(H0_gen, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                  sum(dgamma(Sigma_gen, shape = 0.5, rate = 0.5 * prior$sigma_S0, log = T))
        sum_log <- parallel::mclapply(1:ndraws,
                                      FUN = function(j) {
                                        B <- matrix(B_gen, nrow = K)
                                        A <- a0toA(A_gen, K)
                                        sigma_h <- Sigma_gen
                                        h0 <- H0_gen

                                        ytilde <- as.numeric(A %*% (yt - B %*% xt)) # from dim K * t_max to (Kt_max) * 1
                                        int_h_Gaussian(ytilde = ytilde, h0 = h0, sigma_h = sigma_h, t_max = t_max, K = K)
                                      },
                                      mc.cores = numCores)


      } else {

        if (dist == "Student") {
          w_mean <- apply(W_mat, MARGIN = 2, mean)
          lprior <- mvnfast::dmvn(X = B_gen, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                    mvnfast::dmvn(X = A_gen, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                    sum(dgamma(Sigma_gen, shape = 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                    # sum(dnorm(H0_gen, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                    dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          sigma_h <- Sigma_gen
                                          h0 <- H0_gen
                                          nu <- Nu_gen

                                          llw <- int_h_StudentSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                                h_mean = h_mean, w_mean = w_mean, t_max = t_max, K = K, R = 100)

                                          llw
                                        },
                                        mc.cores = numCores)


          }

        if (dist == "Skew.Student") {
          w_mean <- apply(W_mat, MARGIN = 2, mean)
          lprior <- mvnfast::dmvn(X = B_gen, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                    mvnfast::dmvn(X = A_gen, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                    sum(dgamma(Sigma_gen, shape = 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                    # sum(dnorm(H0_gen, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                    dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T) +
                    mvnfast::dmvn(X = Gamma_gen, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)

          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          gamma <- Gamma_gen

                                          sigma_h <- Sigma_gen
                                          h0 <- H0_gen
                                          nu <- Nu_gen


                                          llw <- int_h_SkewSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                               gamma = gamma, h_mean = h_mean, w_mean = w_mean, t_max = t_max, K = K, R = 1000)


                                          llw
                                        },
                                        mc.cores = numCores)

          }

        if (dist == "MT") {
          w_mean <- matrix(apply(W_mat, MARGIN = 2, mean), nrow = K)
          lprior <- mvnfast::dmvn(X = B_gen, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                    mvnfast::dmvn(X = A_gen, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                    sum(dgamma(Sigma_gen, shape = 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                    # sum(dnorm(H0_gen, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                    sum(dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T))

          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          sigma_h <- Sigma_gen
                                          h0 <- H0_gen
                                          nu <- Nu_gen


                                          llw <- int_h_MSTSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                                    h_mean = h_mean, w_mean = w_mean, t_max = t_max, K = K, R = 100)

                                          llw
                                          },
                                        mc.cores = numCores)

          }

        if (dist == "MST") {
          w_mean <- matrix(apply(W_mat, MARGIN = 2, mean), nrow = K)
          lprior <- mvnfast::dmvn(X = B_gen, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                    mvnfast::dmvn(X = A_gen, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                    sum(dgamma(Sigma_gen, shape = 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                    # sum(dnorm(H0_gen, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                    sum(dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
                    mvnfast::dmvn(X = Gamma_gen, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          gamma <- Gamma_gen

                                          sigma_h <- Sigma_gen
                                          h0 <- H0_gen
                                          nu <- Nu_gen

                                          llw <- int_h_MSTSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                                    gamma = gamma, h_mean = h_mean, w_mean = w_mean, t_max = t_max, K = K, R = 100)

                                          llw
                                        },
                                        mc.cores = numCores)

        }

        if (dist == "OT") {
          w_mean <- matrix(apply(W_mat, MARGIN = 2, mean), nrow = K)
          lprior <- mvnfast::dmvn(X = B_gen, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                    mvnfast::dmvn(X = A_gen, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                    sum(dgamma(Sigma_gen, shape = 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                    # sum(dnorm(H0_gen, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                    sum(dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T))

          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          sigma_h <- Sigma_gen
                                          h0 <- H0_gen
                                          nu <- Nu_gen

                                          llw <- int_h_OTSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                                    h_mean = h_mean, w_mean = w_mean, t_max = t_max, K = K, R = 100)

                                          llw
                                        },
                                        mc.cores = numCores)


        }

        if (dist == "OST") {
          w_mean <- matrix(apply(W_mat, MARGIN = 2, mean), nrow = K)
          lprior <- mvnfast::dmvn(X = B_gen, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                    mvnfast::dmvn(X = A_gen, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                    sum(dgamma(Sigma_gen, shape = 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                    # sum(dnorm(H0_gen, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                    sum(dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
                    mvnfast::dmvn(X = Gamma_gen, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          gamma <- Gamma_gen

                                          sigma_h <- Sigma_gen
                                          h0 <- H0_gen
                                          nu <- Nu_gen


                                          llw <- int_h_OSTSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                                    gamma = gamma, h_mean = h_mean, w_mean = w_mean, t_max = t_max, K = K, R = 100)

                                          llw
                                          },
                                        mc.cores = numCores)

        }


      } # end if student



    } else {


      if (dist == "Gaussian") {
        lprior <- mvnfast::dmvn(X = B_gen, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                  mvnfast::dmvn(X = A_gen, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                  sum(invgamma::dinvgamma(Sigma_gen^2, shape = prior$sigma_T0 * 0.5, rate = 0.5 * prior$sigma_S0, log = T))

        sum_log <- parallel::mclapply(1:ndraws,
                                       FUN = function(j) {
                                         B <- matrix(B_gen, nrow = K)
                                         A <- a0toA(A_gen, K)
                                         sigma <- Sigma_gen
                                         sum(mvnfast::dmvn(X = (y - t(xt) %*% t(B)),
                                                                             mu = rep(0,K),
                                                                             sigma = t(solve(A) %*% diag(sigma, nrow = K)), log = T, isChol = T))
                                       },
                                       mc.cores = numCores)

      } else {
        # Fat tail
        if (dist == "Student") {
          lprior <- mvnfast::dmvn(X = B_gen, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                    mvnfast::dmvn(X = A_gen, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                    sum(invgamma::dinvgamma(Sigma_gen^2, shape = prior$sigma_T0 * 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                    dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          sigma <- Sigma_gen
                                          nu <- Nu_gen

                                          llw <- int_w_TnonSV(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu, t_max = t_max, K = K)

                                          llw
                                        },
                                        mc.cores = numCores)
          }

        if (dist == "Skew.Student") {
          lprior <- mvnfast::dmvn(X = B_gen, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                    mvnfast::dmvn(X = A_gen, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                    sum(invgamma::dinvgamma(Sigma_gen^2, shape = prior$sigma_T0 * 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                    sum(dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
                    mvnfast::dmvn(X = Gamma_gen, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          gamma <- Gamma_gen
                                          sigma <- Sigma_gen
                                          nu <- Nu_gen

                                          llw <- int_w_STnonSV(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu, gamma = gamma,
                                                     t_max = t_max, K = K)
                                          llw
                                        },
                                        mc.cores = numCores)
          }

        if (dist == "MT") {
          lprior <- mvnfast::dmvn(X = B_gen, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                    mvnfast::dmvn(X = A_gen, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                    sum(invgamma::dinvgamma(Sigma_gen^2, shape = prior$sigma_T0 * 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                    sum(dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T))
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          sigma <- Sigma_gen
                                          nu <- Nu_gen

                                          llw <- int_w_MSTnonSV(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu,
                                                                    t_max = t_max, K = K, R = 100)
                                          llw
                                        },
                                        mc.cores = numCores)

          }

        if (dist == "MST") {
          lprior <- mvnfast::dmvn(X = B_gen, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                    mvnfast::dmvn(X = A_gen, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                    sum(invgamma::dinvgamma(Sigma_gen^2, shape = prior$sigma_T0 * 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                    sum(dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
                    mvnfast::dmvn(X = Gamma_gen, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          gamma <- Gamma_gen
                                          sigma <- Sigma_gen
                                          nu <- Nu_gen

                                          llw <- int_w_MSTnonSV(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu, gamma = gamma,
                                                                    t_max = t_max, K = K, R = 100)

                                          llw
                                        },
                                        mc.cores = numCores)
          }

        if (dist == "OT") {
          lprior <- mvnfast::dmvn(X = B_gen, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                    mvnfast::dmvn(X = A_gen, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                    sum(invgamma::dinvgamma(Sigma_gen^2, shape = prior$sigma_T0 * 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                    sum(dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T))
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          sigma <- Sigma_gen
                                          nu <- Nu_gen

                                          llw <- int_w_OSTnonSV(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu,
                                                                 t_max = t_max, K = K)

                                          llw
                                        },
                                        mc.cores = numCores)
        }

        if (dist == "OST") {
          lprior <- mvnfast::dmvn(X = B_gen, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                    mvnfast::dmvn(X = A_gen, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                    sum(invgamma::dinvgamma(Sigma_gen^2, shape = prior$sigma_T0 * 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                    sum(dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
                    mvnfast::dmvn(X = Gamma_gen, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          gamma <- Gamma_gen
                                          sigma <- Sigma_gen
                                          nu <- Nu_gen

                                          llw <- int_w_OSTnonSV(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu, gamma = gamma,
                                                                         t_max = t_max, K = K)

                                          llw
                                        },
                                        mc.cores = numCores)
        }

      } # end Student


    } # end SV

    sum_log <- unlist(sum_log)

    short_sumlog = matrix(sum_log, nrow = ndraws/20, ncol = 20)
    max_sumlog = apply(short_sumlog, MARGIN = 2 , FUN = max)
    bigml = log( apply(exp(short_sumlog-reprow(max_sumlog,ndraws/20)), MARGIN = 2, FUN = mean )) + max_sumlog
    lc = mean(bigml) # log conditional
    mlstd = sd(bigml)/sqrt(20)
  return( list( lprior = as.numeric(lprior),
                LL = as.numeric(lc + lprior),
                lc = lc,
                std = mlstd,
                sum_log = sum_log))
}
