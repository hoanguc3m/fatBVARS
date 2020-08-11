#'
#' @export
harmonicLL <- function(Chain, ndraws = NULL, numCores = NULL){
  if(.Platform$OS.type == "unix") {
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

    B_gen <- get_post(mcmc, element = "B")
    A_gen <- get_post(mcmc, element = "a")
    Gamma_gen <- get_post(mcmc, element = "gamma")
    Sigma_gen <- get_post(mcmc, element = "sigma") # No SV is sigma / SV is sigma^2

    if (SV){
      Sigma_transform <- log(Sigma_gen)
    } else {
      Sigma_transform <- log(Sigma_gen^2)
    }


    H_gen <- get_post(mcmc, element = "h")
    H0_gen <- get_post(mcmc, element = "lh0")
    Nu_gen <- as.matrix(get_post(mcmc, element = "nu"))
    if (ncol(Nu_gen) == 1) colnames(Nu_gen) <- "nu"
    Nu_transform <- log(Nu_gen)

    W_mat <- get_post(mcmc, element = "w")
    W_gen <- get_post(mcmc, element = "w")

    sum_log <- rep(0, ndraws)

    if (is.null(numCores)) numCores <- parallel::detectCores()*0.75

    if (SV){
      h_mean <- matrix(apply(H_gen, MARGIN = 2, mean), nrow = K)

      #################### prior ###########################

      if (dist == "Gaussian") {
        param_trans <- cbind(B_gen, A_gen, H0_gen, Sigma_transform)
        sum_log_prop <- Laplace_approx(param_trans) -
          apply(Sigma_transform, MARGIN = 1, FUN = sum) # Jacob trans

        sum_log <- parallel::mclapply(1:ndraws,
                                      FUN = function(j) {
                                        B <- matrix(B_gen[j,], nrow = K)
                                        A <- a0toA(A_gen[j,], K)
                                        sigma_h <- Sigma_gen[j,]
                                        h0 <- H0_gen[j,]
                                        ytilde <- as.numeric(A %*% (yt - B %*% xt)) # from dim K * t_max to (Kt_max) * 1

                                        int_h(ytilde = ytilde, h0 = h0, sigma_h = sigma_h, t_max = t_max, K = K) +
                                          mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                          mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                          sum(dnorm(sqrt(sigma_h), log = T)) +
                                          sum(dnorm(h0, mean = log(prior$sigma), sd = sqrt(4), log = T))

                                      },
                                      mc.cores = numCores)
      } else {
        w_mean <- matrix(apply(W_mat, MARGIN = 2, mean), nrow = K)

        if (dist == "Student") {
          param_trans <- cbind(B_gen, A_gen, H0_gen, Sigma_transform, Nu_transform)
          sum_log_prop <- Laplace_approx(param_trans) -
            apply(Sigma_transform, MARGIN = 1, FUN = sum) -
            apply(Nu_transform, MARGIN = 1, FUN = sum) # Jacob trans


          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          sigma_h <- Sigma_gen[j,]
                                          h0 <- H0_gen[j,]
                                          nu <- Nu_gen[j]

                                          # llw <- int_h_Student(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                          #                      h_mean = h_mean, t_max = t_max, K = K, R = 100)
                                          llw <- int_h_Student2(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                                h_mean = h_mean, w_mean = w_mean, t_max = t_max, K = K, R = 100)
                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(dnorm(sqrt(sigma_h), log = T)) +
                                            sum(dnorm(h0, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                                            dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)
                                        },
                                        mc.cores = numCores)
        }


        if (dist == "Hyper.Student") {
          param_trans <- cbind(B_gen, A_gen, Gamma_gen, H0_gen, Sigma_transform, Nu_transform)
          sum_log_prop <- Laplace_approx(param_trans) -
            apply(Sigma_transform, MARGIN = 1, FUN = sum) -
            apply(Nu_transform, MARGIN = 1, FUN = sum) # Jacob trans

          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          gamma <- Gamma_gen[j,]

                                          sigma_h <- Sigma_gen[j,]
                                          h0 <- H0_gen[j,]

                                          nu <- Nu_gen[j]

                                          llw <- int_h_Student3(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                               gamma = gamma, h_mean = h_mean, t_max = t_max, K = K, R = 100)

                                          # int_h_Student(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                          #               gamma = gamma, h_mean = h_mean, t_max = t_max, K = K, R = 100)

                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(dnorm(sqrt(sigma_h), log = T)) +
                                            sum(dnorm(h0, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                                            dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T) +
                                            mvnfast::dmvn(X = gamma, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
                                        },
                                        mc.cores = numCores)
        }

        if (dist == "multiStudent") {
          param_trans <- cbind(B_gen, A_gen, H0_gen, Sigma_transform, Nu_transform)
          sum_log_prop <- Laplace_approx(param_trans) -
            apply(Sigma_transform, MARGIN = 1, FUN = sum) -
            apply(Nu_transform, MARGIN = 1, FUN = sum) # Jacob trans

          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          sigma_h <- Sigma_gen[j,]
                                          h0 <- H0_gen[j,]
                                          nu <- Nu_gen[j,]


                                          llw <- int_h_MultiStudent2(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                             h_mean = h_mean, t_max = t_max, K = K, R = 100)

                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(dnorm(sqrt(sigma_h), log = T)) +
                                            sum(dnorm(h0, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                                            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T))
                                        },
                                        mc.cores = numCores)


        }

        if (dist == "Hyper.multiStudent") {
          param_trans <- cbind(B_gen, A_gen, Gamma_gen, H0_gen, Sigma_transform, Nu_transform)
          sum_log_prop <- Laplace_approx(param_trans) -
            apply(Sigma_transform, MARGIN = 1, FUN = sum) -
            apply(Nu_transform, MARGIN = 1, FUN = sum) # Jacob trans
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          gamma <- Gamma_gen[j,]

                                          sigma_h <- Sigma_gen[j,]
                                          h0 <- H0_gen[j,]
                                          nu <- Nu_gen[j,]

                                          llw <- int_h_MultiStudent2(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                                    gamma = gamma, h_mean = h_mean, t_max = t_max, K = K, R = 100)

                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(dnorm(sqrt(sigma_h), log = T)) +
                                            sum(dnorm(h0, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                                            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
                                            mvnfast::dmvn(X = gamma, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
                                        },
                                        mc.cores = numCores)

        }

        if (dist == "multiOrthStudent") {
          param_trans <- cbind(B_gen, A_gen, H0_gen, Sigma_transform, Nu_transform)
          sum_log_prop <- Laplace_approx(param_trans) -
            apply(Sigma_transform, MARGIN = 1, FUN = sum) -
            apply(Nu_transform, MARGIN = 1, FUN = sum) # Jacob trans
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          sigma_h <- Sigma_gen[j,]
                                          h0 <- H0_gen[j,]
                                          nu <- Nu_gen[j,]


                                          llw <- int_h_OrthStudent2(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                                   h_mean = h_mean, t_max = t_max, K = K, R = 100)

                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(dnorm(sqrt(sigma_h), log = T)) +
                                            sum(dnorm(h0, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                                            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T))
                                        },
                                        mc.cores = numCores)


        }

        if (dist == "Hyper.multiOrthStudent") {
          param_trans <- cbind(B_gen, A_gen, Gamma_gen, H0_gen, Sigma_transform, Nu_transform)
          sum_log_prop <- Laplace_approx(param_trans) -
            apply(Sigma_transform, MARGIN = 1, FUN = sum) -
            apply(Nu_transform, MARGIN = 1, FUN = sum) # Jacob trans
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          gamma <- Gamma_gen[j,]

                                          sigma_h <- Sigma_gen[j,]
                                          h0 <- H0_gen[j,]
                                          nu <- Nu_gen[j,]


                                          llw <- int_h_OrthStudent3(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                                   gamma = gamma, h_mean = h_mean, t_max = t_max, K = K, R = 100)

                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(dnorm(sqrt(sigma_h), log = T)) +
                                            sum(dnorm(h0, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                                            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
                                            mvnfast::dmvn(X = gamma, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
                                        },
                                        mc.cores = numCores)

        }


      } # end if student



    } else {


      if (dist == "Gaussian") {
        param_trans <- cbind(B_gen, A_gen, Sigma_transform)
        sum_log_prop <- Laplace_approx(param_trans) -
                        apply(Sigma_transform, MARGIN = 1, FUN = sum) # Jacob trans
        sum_log <- parallel::mclapply(1:ndraws,
                                      FUN = function(j) {
                                        B <- matrix(B_gen[j,], nrow = K)
                                        A <- a0toA(A_gen[j,], K)
                                        sigma <- Sigma_gen[j,]
                                        sum(mvnfast::dmvn(X = (y - t(xt) %*% t(B)),
                                                          mu = rep(0,K),
                                                          sigma = t(solve(A) %*% diag(sigma, nrow = K)), log = T, isChol = T)) +
                                          mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                          mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                          sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5, rate = prior$sigma_S0 * 0.5, log = T))
                                      },
                                      mc.cores = numCores)
      } else {
        # Fat tail

        if (dist == "Student") {
          param_trans <- cbind(B_gen, A_gen, Sigma_transform, Nu_transform)
          sum_log_prop <- Laplace_approx(param_trans) -
            apply(Sigma_transform, MARGIN = 1, FUN = sum) -
            apply(Nu_transform, MARGIN = 1, FUN = sum) # Jacob trans

          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          sigma <- Sigma_gen[j,]
                                          nu <- Nu_gen[j]

                                          llw <- int_w_Student2(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu, t_max = t_max, K = K, R = 100)

                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5,
                                                                    rate = prior$sigma_S0 * 0.5, log = T)) +
                                            dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)
                                        },
                                        mc.cores = numCores)
        }

        if (dist == "Hyper.Student") {
          param_trans <- cbind(B_gen, A_gen, Gamma_gen, Sigma_transform, Nu_transform)
          sum_log_prop <- Laplace_approx(param_trans) -
            apply(Sigma_transform, MARGIN = 1, FUN = sum) -
            apply(Nu_transform, MARGIN = 1, FUN = sum) # Jacob trans

          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          gamma <- Gamma_gen[j,]

                                          sigma <- Sigma_gen[j,]
                                          nu <- Nu_gen[j]

                                          llw <- int_w_Skew(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu, gamma = gamma,
                                                            t_max = t_max, K = K)

                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5,
                                                                    rate = prior$sigma_S0 * 0.5, log = T)) +
                                            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
                                            mvnfast::dmvn(X = gamma, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
                                        },
                                        mc.cores = numCores)
        }

        if (dist == "multiStudent") {
          param_trans <- cbind(B_gen, A_gen, Sigma_transform, Nu_transform)
          sum_log_prop <- Laplace_approx(param_trans) -
            apply(Sigma_transform, MARGIN = 1, FUN = sum) -
            apply(Nu_transform, MARGIN = 1, FUN = sum) # Jacob trans

          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          sigma <- Sigma_gen[j,]
                                          nu <- Nu_gen[j,]

                                          llw <- int_w_MultiStudent2(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu,
                                                                    t_max = t_max, K = K, R = 100)

                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5,
                                                                    rate = prior$sigma_S0 * 0.5, log = T)) +
                                            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T))
                                        },
                                        mc.cores = numCores)
        }

        if (dist == "Hyper.multiStudent") {
          param_trans <- cbind(B_gen, A_gen, Gamma_gen, Sigma_transform, Nu_transform)
          sum_log_prop <- Laplace_approx(param_trans) -
            apply(Sigma_transform, MARGIN = 1, FUN = sum) -
            apply(Nu_transform, MARGIN = 1, FUN = sum) # Jacob trans

          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          gamma <- Gamma_gen[j,]
                                          sigma <- Sigma_gen[j,]
                                          nu <- Nu_gen[j,]

                                          llw <- int_w_MultiStudent2(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu, gamma = gamma,
                                                                    t_max = t_max, K = K, R = 100)

                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5,
                                                                    rate = prior$sigma_S0 * 0.5, log = T)) +
                                            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
                                            mvnfast::dmvn(X = gamma, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
                                        },
                                        mc.cores = numCores)
        }

        if (dist == "multiOrthStudent") {
          param_trans <- cbind(B_gen, A_gen, Sigma_transform, Nu_transform)
          sum_log_prop <- Laplace_approx(param_trans) -
            apply(Sigma_transform, MARGIN = 1, FUN = sum) -
            apply(Nu_transform, MARGIN = 1, FUN = sum) # Jacob trans

          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          sigma <- Sigma_gen[j,]
                                          nu <- Nu_gen[j,]

                                          llw <- int_w_MultiOrthStudent2(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu,
                                                                         t_max = t_max, K = K)

                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5,
                                                                    rate = prior$sigma_S0 * 0.5, log = T)) +
                                            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T))
                                        },
                                        mc.cores = numCores)
        }

        if (dist == "Hyper.multiOrthStudent") {
          param_trans <- cbind(B_gen, A_gen, Gamma_gen, Sigma_transform, Nu_transform)
          sum_log_prop <- Laplace_approx(param_trans) -
            apply(Sigma_transform, MARGIN = 1, FUN = sum) -
            apply(Nu_transform, MARGIN = 1, FUN = sum) # Jacob trans

          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          gamma <- Gamma_gen[j,]
                                          sigma <- Sigma_gen[j,]
                                          nu <- Nu_gen[j,]

                                          llw <- int_w_MultiOrthStudent2(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu, gamma = gamma,
                                                                 t_max = t_max, K = K)

                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5,
                                                                    rate = prior$sigma_S0 * 0.5, log = T)) +
                                            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
                                            mvnfast::dmvn(X = gamma, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
                                        },
                                        mc.cores = numCores)
        }

      } # end Student


    } # end SV

    sum_log <- unlist(sum_log)
    sum_log <- sum_log_prop - sum_log # log posterior - log proposal

    short_sumlog = matrix(sum_log, nrow = ndraws/20, ncol = 20)
    max_sumlog = apply(short_sumlog, MARGIN = 2 , FUN = max)
    bigml = log( apply(exp(short_sumlog-reprow(max_sumlog,ndraws/20)), MARGIN = 2, FUN = mean )) + max_sumlog
    ml = - mean(bigml) # log hamonic mean
    mlstd = sd(bigml)/sqrt(20)

    return( list( LL = ml,
                  std = mlstd))
  } else {
    return( marginalLLSing(Chain, ndraws) )
  }
}

#' @export
Laplace_approx <- function(mcmc_sample, ndraws){
  mcmc_mean <- apply(mcmc_sample, 2, mean)
  mcmc_Sigma <- cov(mcmc_sample)

  sum_log_prop <- mvnfast::dmvn(X = mcmc_sample,
                                mu = as.numeric(mcmc_mean), sigma = mcmc_Sigma, log = T)
  return(sum_log_prop)
}
