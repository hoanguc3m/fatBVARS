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
                                        h0 <- 2 * log(prior$sigma)
                                        # h0 <- H0_gen
                                        ytilde <- as.numeric(A %*% (yt - B %*% xt)) # from dim K * t_max to (Kt_max) * 1
                                        Chibint_h_Gaussian(ytilde = ytilde, h0 = h0, sigma_h = sigma_h, t_max = t_max, K = K)
                                      },
                                      mc.cores = numCores)


      } else {

        if (dist == "Student") {
          w_mean <- apply(W_mat, MARGIN = 2, mean)
          lprior <- mvnfast::dmvn(X = B_gen, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                    mvnfast::dmvn(X = A_gen, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                    sum(dgamma(Sigma_gen, shape = 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                    # sum(dnorm(H0_gen, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                    dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T) +
                    log( 1/(pgamma(100, shape = 2, rate = 0.1) - pgamma(4, shape = 2, rate = 0.1) ) ) # Truncated
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          sigma_h <- Sigma_gen
                                          h0 <- 2 * log(prior$sigma)
                                          # h0 <- H0_gen
                                          nu <- Nu_gen

                                          llw <- Chibint_h_StudentSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
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
                    log( 1/(pgamma(100, shape = 2, rate = 0.1) - pgamma(4, shape = 2, rate = 0.1) ) ) + # Truncated
                    mvnfast::dmvn(X = Gamma_gen, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)

          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          gamma <- Gamma_gen

                                          sigma_h <- Sigma_gen
                                          h0 <- 2 * log(prior$sigma)
                                          # h0 <- H0_gen
                                          nu <- Nu_gen


                                          llw <- Chibint_h_SkewSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
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
                    sum(dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
                    K * log( 1/(pgamma(100, shape = 2, rate = 0.1) - pgamma(4, shape = 2, rate = 0.1) ) ) # Truncated

          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          sigma_h <- Sigma_gen
                                          h0 <- 2 * log(prior$sigma)
                                          # h0 <- H0_gen
                                          nu <- Nu_gen


                                          llw <- Chibint_h_MSTSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
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
                    K * log( 1/(pgamma(100, shape = 2, rate = 0.1) - pgamma(4, shape = 2, rate = 0.1) ) ) + # Truncated
                    mvnfast::dmvn(X = Gamma_gen, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          gamma <- Gamma_gen

                                          sigma_h <- Sigma_gen
                                          h0 <- 2 * log(prior$sigma)
                                          # h0 <- H0_gen
                                          nu <- Nu_gen

                                          llw <- Chibint_h_MSTSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
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
                    sum(dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
                    K * log( 1/(pgamma(100, shape = 2, rate = 0.1) - pgamma(4, shape = 2, rate = 0.1) ) ) # Truncated

          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          sigma_h <- Sigma_gen
                                          h0 <- 2 * log(prior$sigma)
                                          # h0 <- H0_gen
                                          nu <- Nu_gen

                                          llw <- Chibint_h_OTSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
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
                    K * log( 1/(pgamma(100, shape = 2, rate = 0.1) - pgamma(4, shape = 2, rate = 0.1) ) ) + # Truncated
                    mvnfast::dmvn(X = Gamma_gen, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen, nrow = K)
                                          A <- a0toA(A_gen, K)
                                          gamma <- Gamma_gen

                                          sigma_h <- Sigma_gen
                                          h0 <- 2 * log(prior$sigma)
                                          # h0 <- H0_gen
                                          nu <- Nu_gen


                                          llw <- Chibint_h_OSTSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
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
                    dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T) +
                    log( 1/(pgamma(100, shape = 2, rate = 0.1) - pgamma(4, shape = 2, rate = 0.1) ) ) # Truncated
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
                    log( 1/(pgamma(100, shape = 2, rate = 0.1) - pgamma(4, shape = 2, rate = 0.1) ) ) + # Truncated
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
                    sum(dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
                    K * log( 1/(pgamma(100, shape = 2, rate = 0.1) - pgamma(4, shape = 2, rate = 0.1) ) ) # Truncated
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
                    K * log( 1/(pgamma(100, shape = 2, rate = 0.1) - pgamma(4, shape = 2, rate = 0.1) ) ) + # Truncated
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
                    sum(dgamma(Nu_gen, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
                    K * log( 1/(pgamma(100, shape = 2, rate = 0.1) - pgamma(4, shape = 2, rate = 0.1) ) ) # Truncated
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
                    K * log( 1/(pgamma(100, shape = 2, rate = 0.1) - pgamma(4, shape = 2, rate = 0.1) ) ) + # Truncated
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


#' @export
Chibint_h_Gaussian <- function(ytilde, h0, sigma_h, t_max, K, R = 100){
  s2 = ytilde^2
  max_loop = 100
  prior_var_h0 <- 4
  Hh = sparseMatrix(i = 1:(t_max*K),
                    j = 1:(t_max*K),
                    x = rep(1,t_max*K)) -
    sparseMatrix( i = (K+1):(t_max*K),
                  j = 1:((t_max-1)*K),
                  x = rep(1,(t_max-1)*K),
                  dims =  c(t_max*K, t_max*K))
  SH = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = c(1./ (sigma_h + prior_var_h0), rep(1./sigma_h, t_max-1))) # Prior for h1 \sim N(2 log sigma, sigmah^2 + 4 )
  HinvSH_h = Matrix::t(Hh) %*% SH %*% Hh
  alph = Matrix::solve(Hh, sparseMatrix(i = 1:K, j = rep(1,K), x = h0, dims = c(t_max*K,1)))

  e_h = 1
  ht = log(s2+0.001)
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

  c_pri = -t_max*K*0.5*log(2*pi) -.5*(t_max-1)*sum(log(sigma_h)) -.5*sum(log(sigma_h + prior_var_h0))
  c_IS = -t_max*K*0.5*log(2*pi) + sum(log(Matrix::diag(CKh)))

  store_llike = rep(0,R)
  for (i in c(1:R)){
    hc = ht + Matrix::solve(Matrix::t(CKh), rnorm(t_max*K))


    llike = -t_max*K*0.5*log(2*pi) - 0.5*sum(hc) - 0.5 * sum(ytilde^2 * exp(-hc))
    store_llike[i] = as.numeric(llike + (c_pri -.5*Matrix::t(hc-alph)%*%HinvSH_h%*%(hc-alph)) -
                                  (c_IS -.5*Matrix::t(hc-ht)%*%Kh%*%(hc-ht)))

  }

  maxllike = max(store_llike)
  llk = log(mean(exp(store_llike-maxllike))) + maxllike
  return(llk)
}

#' @export
Chibint_h_StudentSV <- function(yt, xt, B, A, h0, sigma_h, nu, gamma = rep(0,K), h_mean, w_mean, t_max, K, R = 100){
  y_tilde = A %*% ( (yt - B %*% xt) )  # Student dist with Identity correlation matrix
  s2 = as.numeric(y_tilde)^2
  max_loop = 10000
  prior_var_h0 <- 4
  Hh = sparseMatrix(i = 1:(t_max*K),
                    j = 1:(t_max*K),
                    x = rep(1,t_max*K)) -
    sparseMatrix( i = (K+1):(t_max*K),
                  j = 1:((t_max-1)*K),
                  x = rep(1,(t_max-1)*K),
                  dims =  c(t_max*K, t_max*K))
  SH = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = c(1./ (sigma_h + prior_var_h0), rep(1./sigma_h, t_max-1))) # Prior for h1 \sim N(2 log sigma, sigmah^2 + 4 )
  HinvSH_h = Matrix::t(Hh) %*% SH %*% Hh
  alph = Matrix::solve(Hh, sparseMatrix(i = 1:K, j = rep(1,K), x = h0, dims = c(t_max*K,1)))

  # e_h = 1
  ht = log(s2+0.001)
  # count = 0

  errh_out = 1;
  while (errh_out> 0.001){
    # E-step
    s2_mod <- as.numeric(reprow(colSums(matrix(s2 * exp(-ht), nrow = K)), K))  # Multivariate student rowSum(y^2 exp(-h))
    Eilam = (nu+K)/(nu + s2_mod)
    s2Eilam = s2*Eilam
    # M-step
    htt = ht
    errh_in = 1
    while (errh_in> 0.001){
      eht = exp(htt)
      sieht = s2Eilam/eht
      fh = -.5 + .5*sieht
      Gh = .5*sieht
      Kh = HinvSH_h + sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = Gh)
      newht = Matrix::solve(Kh, fh+Gh*htt+ HinvSH_h %*% alph)
      errh_in = max(abs(newht-htt))
      htt = newht
    }

    errh_out = max(abs(ht-htt))
    ht = htt
  }

  # compute negative Hessian
  Gh = (nu+K)/(2*nu) * (s2*exp(ht)) / ((exp(ht) + s2/nu )^2)
  Kh = HinvSH_h + sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = Gh)
  CKh = Matrix::t(Matrix::chol(Kh))


  c_pri = -t_max*K*0.5*log(2*pi) -.5*(t_max-1)*sum(log(sigma_h)) -.5*sum(log(sigma_h + prior_var_h0))
  c_IS = -t_max*K*0.5*log(2*pi) + sum(log(Matrix::diag(CKh)))

  store_llike = rep(0,R)

  c_LL <- lgamma(0.5*(nu+K)) - lgamma(0.5*nu) - 0.5*K*log(nu) - 0.5 * K * log(pi)

  for (i in c(1:R)){
    hc = ht + Matrix::solve(Matrix::t(CKh), rnorm(t_max*K))
    hnew <- matrix(hc, nrow = K)
    # allw <- sum( unlist(lapply(1:t_max, FUN = function(t) {
    #   mvnfast::dmvt(X = y_tilde[,t],
    #                 mu = rep(0,K),
    #                 sigma = diag(exp(hnew[,t]/2), nrow = K), df = nu, log = T, isChol = T)
    # })) )
    allw <- c_LL * t_max + sum(- 0.5 * colSums(hnew) - 0.5*(nu+K) * log(1 + colSums(y_tilde^2 * exp(-hnew)) / nu ))
    store_llike[i] = as.numeric(allw + (c_pri -.5*Matrix::t(hc-alph)%*%HinvSH_h%*%(hc-alph)) -
                                  (c_IS -.5*Matrix::t(hc-ht)%*%Kh%*%(hc-ht)))

  }

  maxllike = max(store_llike)
  llk = log(mean(exp(store_llike-maxllike))) + maxllike
  return(llk)
}

#' @export
Chibint_h_SkewSV <- function(yt, xt, B, A, h0, sigma_h, nu, gamma = rep(0,K), h_mean, w_mean, t_max, K, R = 100){
  u <- yt - B %*% xt + nu/(nu-2) * gamma
  y_tilde = A %*% (( u - reprow(w_mean,K)*gamma ) ) # Approximation multivariate Student dist
  s2 = as.numeric(y_tilde)^2
  max_loop = 10000
  prior_var_h0 <- 4
  Hh = sparseMatrix(i = 1:(t_max*K),
                    j = 1:(t_max*K),
                    x = rep(1,t_max*K)) -
    sparseMatrix( i = (K+1):(t_max*K),
                  j = 1:((t_max-1)*K),
                  x = rep(1,(t_max-1)*K),
                  dims =  c(t_max*K, t_max*K))
  SH = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = c(1./ (sigma_h + prior_var_h0), rep(1./sigma_h, t_max-1))) # Prior for h1 \sim N(2 log sigma, sigmah^2 + 4 )
  HinvSH_h = Matrix::t(Hh) %*% SH %*% Hh
  alph = Matrix::solve(Hh, sparseMatrix(i = 1:K, j = rep(1,K), x = h0, dims = c(t_max*K,1)))

  # e_h = 1
  ht = log(s2+0.001)
  # count = 0

  errh_out = 1;
  while (errh_out> 0.001){
    # E-step
    s2_mod <- as.numeric(reprow(colSums(matrix(s2 * exp(-ht), nrow = K)), K))  # Multivariate student rowSum(y^2 exp(-h))
    Eilam = (nu+K)/(nu + s2_mod)
    s2Eilam = s2*Eilam
    # M-step
    htt = ht
    errh_in = 1
    while (errh_in> 0.001){
      eht = exp(htt)
      sieht = s2Eilam/eht
      fh = -.5 + .5*sieht
      Gh = .5*sieht
      Kh = HinvSH_h + sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = Gh)
      newht = Matrix::solve(Kh, fh+Gh*htt+ HinvSH_h %*% alph)
      errh_in = max(abs(newht-htt))
      htt = newht
    }

    errh_out = max(abs(ht-htt))
    ht = htt
  }

  # compute negative Hessian
  Gh = (nu+K)/(2*nu) * (s2*exp(ht)) / ((exp(ht) + s2/nu )^2)
  Kh = HinvSH_h + sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = Gh)
  CKh = Matrix::t(Matrix::chol(Kh))

  c_pri = -t_max*K*0.5*log(2*pi) -.5*(t_max-1)*sum(log(sigma_h)) -.5*sum(log(sigma_h + prior_var_h0))
  c_IS = -t_max*K*0.5*log(2*pi) + sum(log(Matrix::diag(CKh)))

  store_llike = rep(0,R)
  #y_til = A %*% ( (yt - B %*% xt) ) # is Skew Student distribution

  for (i in c(1:R)){
    hc = ht + Matrix::solve(Matrix::t(CKh), rnorm(t_max*K))
    hnew <- matrix(hc, nrow = K)
    # if (K == 1){
    #   allw <-sum( unlist(lapply(1:t_max, FUN = function(t) {
    #     ghyp::dghyp(u[,t], object = ghyp::ghyp(lambda = -0.5*nu, chi = nu, psi = 0, mu = rep(0,K),
    #                                                sigma = diag(exp(hnew[,t]/2), nrow = K),
    #                                                gamma = gamma), logvalue = T) # use sigma for univariable
    #
    #
    #   })) )
    #
    # } else {
    #
    #     allw <-sum( unlist(lapply(1:t_max, FUN = function(t) {
    #       ghyp::dghyp(u[,t], object = ghyp::ghyp(lambda = -0.5*nu, chi = nu, psi = 0, mu = rep(0,K),
    #                                              sigma = solve(A) %*% diag(exp(hnew[,t]), nrow = K) %*% t(solve(A)),
    #                                              gamma = gamma), logvalue = T) # use Sigma = AA' for multivariable
    #     })) )
    # }
    lambda = -0.5*nu; chi = nu; psi = 0;
    inv.vola <- exp(-0.5*hnew)
    xtilde <- A %*% u * inv.vola
    gtilde <- as.vector(A %*% gamma) * inv.vola
    Q.vec <- colSums(xtilde^2)
    skewnorm.vec <- colSums(gtilde^2)
    skewness.scaled.vec <- colSums(xtilde*gtilde)
    interm.vec <- sqrt((chi + Q.vec) * skewnorm.vec)
    lambda.min.d.2 <- lambda - K/2
    log.const.top <- -lambda * log(chi) - lambda.min.d.2 * log(skewnorm.vec)
    log.const.bottom <- K/2 * log(2 * pi) + 0.5 * colSums(hnew) + lgamma(-lambda) - (lambda + 1) * log(2)
    log.top <- log(besselK(interm.vec, lambda.min.d.2, expon.scaled = TRUE)) - interm.vec + skewness.scaled.vec
    log.bottom <- - lambda.min.d.2 * log(interm.vec)
    allw <- sum(log.const.top + log.top - log.const.bottom - log.bottom)



    #
    # mvnfast::dmvt(X = u[,t],
    #               mu = rep(0,K),
    #               sigma = diag(exp(hnew[,t]/2), nrow = K), df = nu, log = T, isChol = T)
    # dt(x = u[,t]/diag(exp(hnew[,t]/2), nrow = K), df = nu, log = T) - hnew[,t]/2


    store_llike[i] = as.numeric(allw + (c_pri -.5*Matrix::t(hc-alph)%*%HinvSH_h%*%(hc-alph)) -
                                  (c_IS -.5*Matrix::t(hc-ht)%*%Kh%*%(hc-ht)))

  }

  maxllike = max(store_llike)
  llk = log(mean(exp(store_llike-maxllike))) + maxllike
  return(llk)
}

#' @export
Chibint_h_MTSV <- function(yt, xt, B, A, h0, sigma_h, nu, h_mean, w_mean, t_max, K, R = 100){
  u <- (yt - B %*% xt)
  y_tilde <- (A %*% u) # Mimic univariate Student
  s2 = as.numeric(y_tilde)^2
  max_loop = 10000
  prior_var_h0 <- 4
  Hh = sparseMatrix(i = 1:(t_max*K),
                    j = 1:(t_max*K),
                    x = rep(1,t_max*K)) -
    sparseMatrix( i = (K+1):(t_max*K),
                  j = 1:((t_max-1)*K),
                  x = rep(1,(t_max-1)*K),
                  dims =  c(t_max*K, t_max*K))
  SH = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = c(1./ (sigma_h + prior_var_h0), rep(1./sigma_h, t_max-1))) # Prior for h1 \sim N(2 log sigma, sigmah^2 + 4 )
  HinvSH_h = Matrix::t(Hh) %*% SH %*% Hh
  alph = Matrix::solve(Hh, sparseMatrix(i = 1:K, j = rep(1,K), x = h0, dims = c(t_max*K,1)))

  ht = log(s2+0.001)
  nu_vec <- rep(nu, t_max)
  errh_out = 1;
  while (errh_out> 0.001){
    # E-step
    s2_mod <- s2 * exp(-ht)  # Univariate student (y^2 exp(-h))
    Eilam = (nu_vec+1)/(nu_vec + s2_mod)
    s2Eilam = s2*Eilam
    # M-step
    htt = ht
    errh_in = 1
    while (errh_in> 0.001){
      eht = exp(htt)
      sieht = s2Eilam/eht
      fh = -.5 + .5*sieht
      Gh = .5*sieht
      Kh = HinvSH_h + sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = Gh)
      newht = Matrix::solve(Kh, fh+Gh*htt+ HinvSH_h %*% alph)
      errh_in = max(abs(newht-htt))
      htt = newht
    }

    errh_out = max(abs(ht-htt))
    ht = htt
  }

  # compute negative Hessian
  Gh = (nu_vec+1)/(2*nu_vec) * (s2*exp(ht)) / ((exp(ht) + s2/nu_vec )^2)
  Kh = HinvSH_h + sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = Gh)
  CKh = Matrix::t(Matrix::chol(Kh))

  c_pri = -t_max*K*0.5*log(2*pi) -.5*(t_max-1)*sum(log(sigma_h)) -.5*sum(log(sigma_h + prior_var_h0))
  c_IS = -t_max*K*0.5*log(2*pi) + sum(log(Matrix::diag(CKh)))

  store_llw = rep(0,R)
  for (i in c(1:R)){
    hc = ht + Matrix::solve(Matrix::t(CKh), rnorm(t_max*K))
    hnew <- matrix(hc, nrow = K)
    allw <- int_w_MTnonSV(y = t(yt), xt = xt, A = A, B = B, sigma = exp(hnew/2), nu = nu,
                          t_max = t_max, K = K, R = 100)
    store_llw[i] = as.numeric(allw + (c_pri -.5*Matrix::t(hc-alph)%*%HinvSH_h%*%(hc-alph)) -
                                (c_IS -.5*Matrix::t(hc-ht)%*%Kh%*%(hc-ht)))
  }

  maxllike = max(store_llw)
  llk = log(mean(exp(store_llw-maxllike))) + maxllike
  return(llk)
}

#' @export
Chibint_h_MSTSV <- function(yt, xt, B, A, h0, sigma_h, nu, gamma = rep(0,K), h_mean, w_mean, t_max, K, R = 100){
  u <- (yt - B %*% xt - w_mean * gamma + nu/(nu-2)*gamma)
  y_tilde <- (A %*% u) # Mimic Univariate Student
  s2 = as.numeric(y_tilde)^2
  max_loop = 10000
  prior_var_h0 <- 4
  Hh = sparseMatrix(i = 1:(t_max*K),
                    j = 1:(t_max*K),
                    x = rep(1,t_max*K)) -
    sparseMatrix( i = (K+1):(t_max*K),
                  j = 1:((t_max-1)*K),
                  x = rep(1,(t_max-1)*K),
                  dims =  c(t_max*K, t_max*K))
  SH = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = c(1./ (sigma_h + prior_var_h0), rep(1./sigma_h, t_max-1))) # Prior for h1 \sim N(2 log sigma, sigmah^2 + 4 )
  HinvSH_h = Matrix::t(Hh) %*% SH %*% Hh
  alph = Matrix::solve(Hh, sparseMatrix(i = 1:K, j = rep(1,K), x = h0, dims = c(t_max*K,1)))

  ht = log(s2+0.001)
  nu_vec <- rep(nu, t_max)
  errh_out = 1;
  while (errh_out> 0.001){
    # E-step
    s2_mod <- s2 * exp(-ht)  # Univariate student (y^2 exp(-h))
    Eilam = (nu_vec+1)/(nu_vec + s2_mod)
    s2Eilam = s2*Eilam
    # M-step
    htt = ht
    errh_in = 1
    while (errh_in> 0.001){
      eht = exp(htt)
      sieht = s2Eilam/eht
      fh = -.5 + .5*sieht
      Gh = .5*sieht
      Kh = HinvSH_h + sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = Gh)
      newht = Matrix::solve(Kh, fh+Gh*htt+ HinvSH_h %*% alph)
      errh_in = max(abs(newht-htt))
      htt = newht
    }

    errh_out = max(abs(ht-htt))
    ht = htt
  }

  # compute negative Hessian
  Gh = (nu_vec+1)/(2*nu_vec) * (s2*exp(ht)) / ((exp(ht) + s2/nu_vec )^2)
  Kh = HinvSH_h + sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = Gh)
  CKh = Matrix::t(Matrix::chol(Kh))

  c_pri = -t_max*K*0.5*log(2*pi) -.5*(t_max-1)*sum(log(sigma_h)) -.5*sum(log(sigma_h + prior_var_h0))
  c_IS = -t_max*K*0.5*log(2*pi) + sum(log(Matrix::diag(CKh)))

  store_llw = rep(0,R)
  for (i in c(1:R)){
    hc = ht + Matrix::solve(Matrix::t(CKh), rnorm(t_max*K))
    hnew <- matrix(hc, nrow = K)
    allw <- int_w_MSTnonSV(y = t(yt), xt = xt, A = A, B = B, sigma = exp(hnew/2), nu = nu, gamma = gamma,
                           t_max = t_max, K = K, R = 100)
    store_llw[i] = as.numeric(allw + (c_pri -.5*Matrix::t(hc-alph)%*%HinvSH_h%*%(hc-alph)) -
                                (c_IS -.5*Matrix::t(hc-ht)%*%Kh%*%(hc-ht)))
  }

  maxllike = max(store_llw)
  llk = log(mean(exp(store_llw-maxllike))) + maxllike
  return(llk)
}

#' @export
Chibint_h_OTSV <- function(yt, xt, B, A, h0, sigma_h, nu, h_mean, w_mean, t_max, K, R = 100){
  u <- (yt - B %*% xt)
  y_tilde <- (A %*% u) # Univarite Student
  s2 = as.numeric(y_tilde)^2
  max_loop = 10000
  prior_var_h0 <- 4
  Hh = sparseMatrix(i = 1:(t_max*K),
                    j = 1:(t_max*K),
                    x = rep(1,t_max*K)) -
    sparseMatrix( i = (K+1):(t_max*K),
                  j = 1:((t_max-1)*K),
                  x = rep(1,(t_max-1)*K),
                  dims =  c(t_max*K, t_max*K))
  SH = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = c(1./ (sigma_h + prior_var_h0), rep(1./sigma_h, t_max-1))) # Prior for h1 \sim N(2 log sigma, sigmah^2 + 4 )
  HinvSH_h = Matrix::t(Hh) %*% SH %*% Hh
  alph = Matrix::solve(Hh, sparseMatrix(i = 1:K, j = rep(1,K), x = h0, dims = c(t_max*K,1)))

  # e_h = 1
  ht = log(s2+0.001)
  # count = 0
  nu_vec <- rep(nu, t_max)
  errh_out = 1;
  while (errh_out> 0.001){
    # E-step
    s2_mod <- s2 * exp(-ht)  # Univariate student (y^2 exp(-h))
    Eilam = (nu_vec+1)/(nu_vec + s2_mod)
    s2Eilam = s2*Eilam
    # M-step
    htt = ht
    errh_in = 1
    while (errh_in> 0.001){
      eht = exp(htt)
      sieht = s2Eilam/eht
      fh = -.5 + .5*sieht
      Gh = .5*sieht
      Kh = HinvSH_h + sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = Gh)
      newht = Matrix::solve(Kh, fh+Gh*htt+ HinvSH_h %*% alph)
      errh_in = max(abs(newht-htt))
      htt = newht
    }

    errh_out = max(abs(ht-htt))
    ht = htt
  }

  # compute negative Hessian
  Gh = (nu_vec+1)/(2*nu_vec) * (s2*exp(ht)) / ((exp(ht) + s2/nu_vec )^2)
  Kh = HinvSH_h + sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = Gh)
  CKh = Matrix::t(Matrix::chol(Kh))


  c_pri = -t_max*K*0.5*log(2*pi) -.5*(t_max-1)*sum(log(sigma_h)) -.5*sum(log(sigma_h + prior_var_h0))
  c_IS = -t_max*K*0.5*log(2*pi) + sum(log(Matrix::diag(CKh)))

  store_llike = rep(0,R)
  y_tilde = A %*% ( (yt - B %*% xt) ) # is ortho-Student distribution
  c_LL <- lgamma(0.5*(nu_vec+1)) - lgamma(0.5*nu_vec) - 0.5 * log(nu_vec) - 0.5 * log(pi)
  for (i in c(1:R)){
    hc = ht + Matrix::solve(Matrix::t(CKh), rnorm(t_max*K))
    hnew <- matrix(hc, nrow = K)
    # allw <- rep(0,K)
    # for (k in c(1:K)){
    #   allw[k] <- sum(dt(y_tilde[k,]/exp(hnew[k,]/2), df = nu[k], log = T)) - 0.5 * sum(hnew[k,])
    # }
    store_llike[i] = as.numeric( sum(c_LL - 0.5 * (nu_vec+1) * log(1 + y_tilde^2*exp(-hnew) / nu_vec ) - 0.5 * hnew ) +
                                   (c_pri -.5*Matrix::t(hc-alph)%*%HinvSH_h%*%(hc-alph)) -
                                   (c_IS -.5*Matrix::t(hc-ht)%*%Kh%*%(hc-ht)))

  }

  maxllike = max(store_llike)
  llk = log(mean(exp(store_llike-maxllike))) + maxllike

  return(llk)
}

#' @export
Chibint_h_OSTSV <- function(yt, xt, B, A, h0, sigma_h, nu, gamma = rep(0,K), h_mean, w_mean, t_max, K, R = 100){
  u <- (yt - B %*% xt)
  y_tilde <- (A %*% u - w_mean * gamma + nu/(nu-2)*gamma) # Mimic univariate Student distribution
  s2 = as.numeric(y_tilde)^2
  max_loop = 10000
  prior_var_h0 <- 4
  Hh = sparseMatrix(i = 1:(t_max*K),
                    j = 1:(t_max*K),
                    x = rep(1,t_max*K)) -
    sparseMatrix( i = (K+1):(t_max*K),
                  j = 1:((t_max-1)*K),
                  x = rep(1,(t_max-1)*K),
                  dims =  c(t_max*K, t_max*K))
  SH = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = c(1./ (sigma_h + prior_var_h0), rep(1./sigma_h, t_max-1))) # Prior for h1 \sim N(2 log sigma, sigmah^2 + 4 )
  HinvSH_h = Matrix::t(Hh) %*% SH %*% Hh
  alph = Matrix::solve(Hh, sparseMatrix(i = 1:K, j = rep(1,K), x = h0, dims = c(t_max*K,1)))

  ht = log(s2+0.001)
  nu_vec <- rep(nu, t_max)
  errh_out = 1;
  while (errh_out> 0.001){
    # E-step
    s2_mod <- s2 * exp(-ht)  # Univariate student (y^2 exp(-h))
    Eilam = (nu_vec+1)/(nu_vec + s2_mod)
    s2Eilam = s2*Eilam
    # M-step
    htt = ht
    errh_in = 1
    while (errh_in> 0.001){
      eht = exp(htt)
      sieht = s2Eilam/eht
      fh = -.5 + .5*sieht
      Gh = .5*sieht
      Kh = HinvSH_h + sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = Gh)
      newht = Matrix::solve(Kh, fh+Gh*htt+ HinvSH_h %*% alph)
      errh_in = max(abs(newht-htt))
      htt = newht
    }

    errh_out = max(abs(ht-htt))
    ht = htt
  }

  # compute negative Hessian
  Gh = (nu_vec+1)/(2*nu_vec) * (s2*exp(ht)) / ((exp(ht) + s2/nu_vec )^2)
  Kh = HinvSH_h + sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = Gh)
  CKh = Matrix::t(Matrix::chol(Kh))

  c_pri = -t_max*K*0.5*log(2*pi) -.5*(t_max-1)*sum(log(sigma_h)) -.5*sum(log(sigma_h + prior_var_h0))
  c_IS = -t_max*K*0.5*log(2*pi) + sum(log(Matrix::diag(CKh)))

  store_llike = rep(0,R)
  y_tilde = A %*% ( (yt - B %*% xt) ) + nu/(nu-2)*gamma # is ortho-Skew Student distribution
  for (i in c(1:R)){
    hc = ht + Matrix::solve(Matrix::t(CKh), rnorm(t_max*K))
    hnew <- matrix(hc, nrow = K)
    # allw <- matrix(NA, nrow = K, ncol = t_max)
    # for (k in c(1:K)){
    #   for (t in c(1:t_max)){
    #     allw[k,t] <- ghyp::dghyp(y_tilde[k,t], object = ghyp::ghyp(lambda = -0.5*nu[k], chi = nu[k], psi = 0, mu = 0,
    #                                                                gamma = gamma[k], sigma = exp(hnew[k,t]/2)),
    #                              logvalue = T)
    #   }
    # }

    # Skew Student density
    lambda = -0.5*nu; chi = nu; psi = 0;
    sigma = exp(hnew/2);
    x = as.vector(y_tilde)

    sigma <- as.vector(sigma)

    Q <- (x/sigma)^2
    d <- 1
    n <- length(x)
    inv.sigma <- 1/sigma^2

    skewness.scaled <- x * (inv.sigma * gamma)
    skewness.norm <- gamma^2 * inv.sigma
    lambda.min.0.5 <- lambda - 0.5
    interm <- sqrt((chi + Q) * as.vector(skewness.norm))
    log.const.top <- -lambda * log(chi) - lambda.min.0.5 * log(as.vector(skewness.norm))
    log.const.bottom <- 0.5 * log(2 * pi) + log(sigma) + lgamma(-lambda) - (lambda + 1) * log(2)
    log.top <- log(besselK(interm, lambda.min.0.5, expon.scaled = TRUE)) - interm + skewness.scaled
    log.bottom <- -lambda.min.0.5 * log(interm)
    allw <- log.const.top + log.top - log.const.bottom - log.bottom

    store_llike[i] = as.numeric( sum(allw) + (c_pri -.5*Matrix::t(hc-alph)%*%HinvSH_h%*%(hc-alph)) -
                                   (c_IS -.5*Matrix::t(hc-ht)%*%Kh%*%(hc-ht)))

  }

  maxllike = max(store_llike)
  llk = log(mean(exp(store_llike-maxllike))) + maxllike

  return(llk)
}
