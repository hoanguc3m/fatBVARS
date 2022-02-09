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
marginalLL <- function(Chain, ndraws = NULL, numCores = NULL){
  if(.Platform$OS.type == "unix") {
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
    BAG_samples <- Normal_approx(cbind(B_mat,A_mat, Gamma_mat), ndraws = ndraws)
    B_gen <- get_post(BAG_samples$new_samples, element = "B")
    A_gen <- get_post(BAG_samples$new_samples, element = "a")
    Gamma_gen <- get_post(BAG_samples$new_samples, element = "gamma")
    sum_log_prop <- BAG_samples$sum_log_prop

    Sigma_mat <- get_post(mcmc, element = "sigma") # No SV is sigma / SV is sigma^2
    # Sigma_H <- sqrt(get_post(mcmc, element = "sigma_h")) # SV
    if (SV) {

      # Gamma proposal
      Sigma_gen_list <- Gamma_approx(mcmc_sample = Sigma_mat, ndraws = ndraws)
      Sigma_gen <- Sigma_gen_list$new_samples

      # InvGamma proposal
      # Sigma_gen_list_inv <- InvGamma_approx(Sigma_mat, ndraws = ndraws)
      # Sigma_gen_inv <- Sigma_gen_list_inv$new_samples

      # Normal proposal
      # Sigma_gen_list <- Normal_approx(log(Sigma_mat), ndraws = ndraws) # Change to normal
      # Sigma_gen <- exp(Sigma_gen_list$new_samples) # Change to square
      # Sigma_gen_list$sum_log_prop <- Sigma_gen_list$sum_log_prop - apply(Sigma_gen_list$new_samples, MARGIN = 1, FUN = sum) # Jacobian


    } else {
      Sigma_gen_list <- Gamma_approx(mcmc_sample = Sigma_mat^2, ndraws = ndraws)
      Sigma_gen <- sqrt(Sigma_gen_list$new_samples)

      # InvGamma proposal
      # Sigma_gen_list <- InvGamma_approx(mcmc_sample = Sigma_mat^2, ndraws = ndraws)
      # Sigma_gen <- sqrt(Sigma_gen_list$new_samples)

      # Normal proposal
      # Sigma_gen_list <- Normal_approx(2*log(Sigma_mat), ndraws = ndraws) # Change to normal
      # Sigma_gen <- exp(0.5*Sigma_gen_list$new_samples)
      # Sigma_gen_list$sum_log_prop <- Sigma_gen_list$sum_log_prop - apply(Sigma_gen_list$new_samples, MARGIN = 1, FUN = sum) # Jacobian
    }
    sum_log_prop <- sum_log_prop + Sigma_gen_list$sum_log_prop

    H_mat <- get_post(mcmc, element = "h")
    H0_mat <- get_post(mcmc, element = "lh0")
    Nu_mat <- as.matrix(get_post(mcmc, element = "nu"))
    if (ncol(Nu_mat) == 1) colnames(Nu_mat) <- "nu"
    W_mat <- get_post(mcmc, element = "w")

    sum_log <- rep(0, ndraws)

    if (is.null(numCores)) numCores <- parallel::detectCores()*0.5

    if (SV){
      h_mean <- matrix(apply(H_mat, MARGIN = 2, mean), nrow = K)
      H0_gen_list <- Normal_approx(H0_mat, ndraws = ndraws)
      H0_gen <- H0_gen_list$new_samples
      sum_log_prop <- sum_log_prop + H0_gen_list$sum_log_prop

      #################### End prior ###########################

      if (dist == "Gaussian") {
        # param <- cbind(B_gen, A_gen, Sigma_gen)

        sum_log <- parallel::mclapply(1:ndraws,
                                      FUN = function(j) {
                                        B <- matrix(B_gen[j,], nrow = K)
                                        A <- a0toA(A_gen[j,], K)
                                        sigma_h <- Sigma_gen[j,]
                                        h0 <- H0_gen[j,]

                                        ytilde <- as.numeric(A %*% (yt - B %*% xt)) # from dim K * t_max to (Kt_max) * 1
                                        int_h_Gaussian(ytilde = ytilde, h0 = h0, sigma_h = sigma_h, t_max = t_max, K = K) +
                                          mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                          mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                          sum(dgamma(sigma_h, shape = 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                                          #sum(dnorm(sqrt(sigma_h), log = T)) +
                                          sum(dnorm(h0, mean = log(prior$sigma), sd = sqrt(4), log = T))
                                      },
                                      mc.cores = numCores)
      } else {
        Nu_gen_list <- Nu_Gamma_approx(Nu_mat, ndraws = ndraws, dist)
        Nu_gen <- Nu_gen_list$new_samples
        sum_log_prop <- sum_log_prop + Nu_gen_list$sum_log_prop
        # w_mean <- matrix(apply(W_mat, MARGIN = 2, mean), nrow = K)

        if (dist == "Student") {
          w_mean <- apply(W_mat, MARGIN = 2, mean)
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          sigma_h <- Sigma_gen[j,]
                                          h0 <- H0_gen[j,]
                                          nu <- Nu_gen[j]

                                          llw <- int_h_StudentSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                                h_mean = h_mean, w_mean = w_mean, t_max = t_max, K = K, R = 100)

                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(dgamma(sigma_h, shape = 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                                            #sum(dnorm(sqrt(sigma_h), log = T)) +
                                            sum(dnorm(h0, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                                            dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)
                                        },
                                        mc.cores = numCores)
          }


        if (dist == "Skew.Student") {
          w_mean <- apply(W_mat, MARGIN = 2, mean)
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          gamma <- Gamma_gen[j,]

                                          sigma_h <- Sigma_gen[j,]
                                          h0 <- H0_gen[j,]
                                          nu <- Nu_gen[j]


                                          llw <- int_h_SkewSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                               gamma = gamma, h_mean = h_mean, w_mean = w_mean, t_max = t_max, K = K, R = 100)

                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(dgamma(sigma_h, shape = 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                                            #sum(dnorm(sqrt(sigma_h), log = T)) +
                                            sum(dnorm(h0, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                                            dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T) +
                                            mvnfast::dmvn(X = gamma, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
                                        },
                                        mc.cores = numCores)
          }

        if (dist == "MT") {
          w_mean <- matrix(apply(W_mat, MARGIN = 2, mean), nrow = K)
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          sigma_h <- Sigma_gen[j,]
                                          h0 <- H0_gen[j,]
                                          nu <- Nu_gen[j,]


                                          llw <- int_h_MTSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                                    h_mean = h_mean, w_mean = w_mean, t_max = t_max, K = K, R = 100)

                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(dgamma(sigma_h, shape = 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                                            #sum(dnorm(sqrt(sigma_h), log = T)) +
                                            sum(dnorm(h0, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                                            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T))
                                        },
                                        mc.cores = numCores)


          }

        if (dist == "MST") {
          w_mean <- matrix(apply(W_mat, MARGIN = 2, mean), nrow = K)
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          gamma <- Gamma_gen[j,]

                                          sigma_h <- Sigma_gen[j,]
                                          h0 <- H0_gen[j,]
                                          nu <- Nu_gen[j,]

                                          llw <- int_h_MSTSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                                    gamma = gamma, h_mean = h_mean, w_mean = w_mean, t_max = t_max, K = K, R = 100)

                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(dgamma(sigma_h, shape = 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                                            #sum(dnorm(sqrt(sigma_h), log = T)) +
                                            sum(dnorm(h0, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                                            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
                                            mvnfast::dmvn(X = gamma, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
                                        },
                                        mc.cores = numCores)

        }

        if (dist == "OT") {
          w_mean <- matrix(apply(W_mat, MARGIN = 2, mean), nrow = K)
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          sigma_h <- Sigma_gen[j,]
                                          h0 <- H0_gen[j,]
                                          nu <- Nu_gen[j,]


                                          llw <- int_h_OTSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                                    h_mean = h_mean, w_mean = w_mean, t_max = t_max, K = K, R = 100)

                                          # ytilde <- A %*% (yt - B %*% xt) # dim K * t_max
                                          # lluni <- rep(0,K)
                                          # for (i in c(1:K)){
                                          #   lluni[i] <- PF_StudentSV(ytilde = ytilde[i,], h0 = h0[i], sigma_h = sigma_h[i], nu = nu[i], noParticles = 1000)$logLikelihood
                                          # }
                                          # sum(lluni) +
                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(dgamma(sigma_h, shape = 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                                            #sum(dnorm(sqrt(sigma_h), log = T)) +
                                            sum(dnorm(h0, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                                            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T))
                                        },
                                        mc.cores = numCores)


        }

        if (dist == "OST") {
          w_mean <- matrix(apply(W_mat, MARGIN = 2, mean), nrow = K)
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          gamma <- Gamma_gen[j,]

                                          sigma_h <- Sigma_gen[j,]
                                          h0 <- H0_gen[j,]
                                          nu <- Nu_gen[j,]


                                          llw <- int_h_OSTSV(yt = yt, xt = xt, B = B, A = A, h0 = h0, sigma_h = sigma_h, nu = nu,
                                                                    gamma = gamma, h_mean = h_mean, w_mean = w_mean, t_max = t_max, K = K, R = 100)

                                          # ytilde <- A %*% (yt - B %*% xt) # dim K * t_max
                                          # lluni <- rep(0,K)
                                          # for (i in c(1:K)){
                                          #   lluni[i] <- PF_uniSkewSV(ytilde = ytilde[i,], h0 = h0[i], sigma_h = sigma_h[i], nu = nu[i], gamma = gamma[i], noParticles = 1000)$logLikelihood
                                          #   #PF_StudentSV(ytilde = ytilde[i,], h0 = h0[i], sigma_h = sigma_h[i], nu = nu[i], noParticles = 1000)$logLikelihood
                                          # }
                                          # sum(lluni) +
                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(dgamma(sigma_h, shape = 0.5, rate = 0.5 * prior$sigma_S0, log = T)) +
                                            #sum(dnorm(sqrt(sigma_h), log = T)) +
                                            sum(dnorm(h0, mean = log(prior$sigma), sd = sqrt(4), log = T)) +
                                            sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
                                            mvnfast::dmvn(X = gamma, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
                                        },
                                        mc.cores = numCores)

        }


      } # end if student



    } else {


      if (dist == "Gaussian") {
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

        # sum_log_inv <- parallel::mclapply(1:ndraws,
        #                               FUN = function(j) {
        #                                 B <- matrix(B_gen[j,], nrow = K)
        #                                 A <- a0toA(A_gen[j,], K)
        #                                 sigma <- Sigma_gen_inv[j,]
        #                                 sum(mvnfast::dmvn(X = (y - t(xt) %*% t(B)),
        #                                                   mu = rep(0,K),
        #                                                   sigma = t(solve(A) %*% diag(sigma, nrow = K)), log = T, isChol = T)) +
        #                                   mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
        #                                   mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
        #                                   sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5, rate = prior$sigma_S0 * 0.5, log = T))
        #                               },
        #                               mc.cores = numCores)
      } else {
        # Fat tail
        Nu_gen_list <- Nu_Gamma_approx(Nu_mat, ndraws = ndraws, dist)
        Nu_gen <- Nu_gen_list$new_samples
        sum_log_prop <- sum_log_prop + Nu_gen_list$sum_log_prop

        if (dist == "Student") {
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          sigma <- Sigma_gen[j,]
                                          nu <- Nu_gen[j]

                                          # llw <- int_w_Student(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu,
                                          #                     t_max = t_max, K = K, R = 1)
                                          llw <- int_w_TnonSV(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu, t_max = t_max, K = K)

                                          llw +
                                            mvnfast::dmvn(X = B_gen[j,], mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
                                            mvnfast::dmvn(X = A_gen[j,], mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
                                            sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5,
                                                                    rate = prior$sigma_S0 * 0.5, log = T)) +
                                            dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)
                                        },
                                        mc.cores = numCores)
          }

        if (dist == "Skew.Student") {
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          gamma <- Gamma_gen[j,]

                                          sigma <- Sigma_gen[j,]
                                          nu <- Nu_gen[j]

                                          llw <- int_w_STnonSV(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu, gamma = gamma,
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

        if (dist == "MT") {
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          sigma <- Sigma_gen[j,]
                                          nu <- Nu_gen[j,]

                                          llw <- int_w_MTnonSV(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu,
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

        if (dist == "MST") {
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          gamma <- Gamma_gen[j,]
                                          sigma <- Sigma_gen[j,]
                                          nu <- Nu_gen[j,]

                                          llw <- int_w_MSTnonSV(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu, gamma = gamma,
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

        if (dist == "OT") {
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          sigma <- Sigma_gen[j,]
                                          nu <- Nu_gen[j,]

                                          llw <- int_w_OSTnonSV(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu,
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

        if (dist == "OST") {
          sum_log <- parallel::mclapply(1:ndraws,
                                        FUN = function(j) {
                                          B <- matrix(B_gen[j,], nrow = K)
                                          A <- a0toA(A_gen[j,], K)
                                          gamma <- Gamma_gen[j,]
                                          sigma <- Sigma_gen[j,]
                                          nu <- Nu_gen[j,]

                                          llw <- int_w_OSTnonSV(y = y, xt = xt, A = A, B = B, sigma = sigma, nu = nu, gamma = gamma,
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
    sum_log <- sum_log - sum_log_prop # log posterior - log proposal

    short_sumlog = matrix(sum_log, nrow = ndraws/20, ncol = 20)
    max_sumlog = apply(short_sumlog, MARGIN = 2 , FUN = max)
    bigml = log( apply(exp(short_sumlog-reprow(max_sumlog,ndraws/20)), MARGIN = 2, FUN = mean )) + max_sumlog
    ml = mean(bigml)
    mlstd = sd(bigml)/sqrt(20)
    # impt <- importance_check(sum_log)
    # ml_est <- TrimmedML_est(sum_log)
  return( list( LL = ml,
                std = mlstd,
                sum_log = sum_log))
  } else {
    return( marginalLLSing(Chain, ndraws) )
  }
}

#' @export
TrimmedML_est <- function(sum_log, turnk = 50, nbatch = 10){
  # nbatch <- 1
  # turnk <- 50
  ml_batch <- rep(NA, nbatch)
  ml_mean <- rep(NA, nbatch)
  msum_log <- matrix(sum_log, ncol = nbatch)
  for (i in c(1:nbatch)){
    log_w <- msum_log[,i]
    maxlogw <- max(log_w)
    log_w <- log_w - maxlogw
    impt_w <- exp(log_w)

    N <- length(impt_w)
    m_n <- round((N)^(turnk/100) * log(N), digits = 0)
    if (m_n > N) cat("Error m_n > N")
    k_n <- round((N)^(turnk/100), digits = 0)
    sort_weights <- sort(impt_w, decreasing = TRUE)
    l_n <- sort_weights[k_n]
    impt_w_star <- impt_w * (impt_w <= l_n)
    phi_star <- mean(impt_w_star)
    log_weight_s <- log(sort_weights[c(1:m_n)])
    alpha_inv <- 1/ m_n * sum(log_weight_s - log(sort_weights[m_n]))
    if (alpha_inv > 1) {
      alpha_inv <- min(sapply(10:(min(5000,m_n*0.5)), FUN = function(k) gPdtest::gpd.fit(sort_weights[c(1:k)],"amle")[1] ))
      B_n <- 1 / (1 - alpha_inv) * k_n / N * l_n
    } else {
      B_n <- 1 / (1 - alpha_inv) * k_n / N * l_n
    }

    phi_ubias <- phi_star + B_n
    ml_batch[i] <- maxlogw + log(phi_ubias)
    ml_mean[i] <- maxlogw + log(mean(impt_w))

  }
  ml_batch <- ml_batch[!is.na(ml_batch)]
  if (length(ml_batch) > 1){
    ml = mean(ml_batch)
    mlstd = sd(ml_batch)/sqrt(length(ml_batch))
    mlm = mean(ml_mean)
    mlmstd = sd(ml_mean)/sqrt(length(ml_mean))
  } else {
    ml = mean(ml_batch)
    T_n = matrix(c(1, - 1/(1/alpha_inv -1) * sqrt(k_n) / sqrt(N) * l_n ), ncol = 1)
    Yta_n = cbind(impt_w_star - mean(impt_w_star), sqrt(N)/ sqrt(k_n) * ((impt_w > l_n) - k_n / N)  )
    Omega_n <- t(Yta_n) %*% Yta_n # / N # ???
    mlstd = as.numeric( sqrt(t(T_n) %*% Omega_n %*% T_n) ) # / sqrt(N) # ???
    mlm = mean(ml_mean)
    if (length(ml_mean) > 1) {
      mlmstd = sd(ml_mean)
    } else {
      mlmstd = NULL
    }
  }
  return(list(ml = ml, mlstd = mlstd, mlm = mlm, mlmstd = mlmstd))
}

# No use
int_w_Student <- function(y, xt, A, B, sigma, nu, gamma = rep(0,K), t_max, K, R = 100){
  u <- (y - t(xt) %*% t(B))
  ytilde <- (A %*% t(u))
  shape <- nu*0.5 + K*0.5
  rate <- nu*0.5 + 0.5 * apply(ytilde^2/sigma^2, 2, sum)

  store_llw <- matrix(NA, nrow = t_max, ncol = R)
  for (i in c(1:R)){
    w <- matrix(mapply(FUN = rinvgamma, n = 1, shape = shape, rate = rate), ncol = 1)
    dens_w <- mapply(FUN = dinvgamma,x = w, shape = shape, rate = rate, log = T)
    Prior_dens_w <- mapply(FUN = dinvgamma,x = w, shape = nu * 0.5, rate = nu * 0.5, log = T)

    store_llw[,i] <- mvnfast::dmvn(X = (u - repcol(w,K) * reprow(gamma, t_max))/repcol(sqrt(w),K),
                                   mu = rep(0,K),
                                   sigma = t(solve(A) %*% diag(sigma, nrow = K)), log = T, isChol = T) -
      as.numeric(0.5 * K * log(w)) + Prior_dens_w - dens_w # proposal
  }
  maxllike = apply(store_llw, MARGIN = 1, max)
  llw = log( apply(exp(store_llw-maxllike), MARGIN = 1, FUN = mean) ) + maxllike
  return(sum(llw))
}

#' @export
int_w_TnonSV <- function(y, xt, A, B, sigma, nu, t_max, K, R = 100){
  u <- (y - t(xt) %*% t(B))
  allw <- mvnfast::dmvt(X = u,
                        mu = rep(0,K),
                        sigma = t(solve(A) %*% diag(sigma, nrow = K)), df = nu, log = T, isChol = T)

  return(sum(allw))
}

#' @export
int_w_STnonSV <- function(y, xt, A, B, sigma, nu, gamma, t_max, K, R = 100){
  u <- (y - t(xt) %*% t(B))
  if (K == 1){
    allw <- ghyp::dghyp(u, object = ghyp::ghyp(lambda = -0.5*nu, chi = nu, psi = 0, mu = -nu/(nu-2)*gamma,
                                               sigma = diag(sigma, nrow = K),
                                               gamma = gamma), logvalue = T)
  } else {
    allw <- ghyp::dghyp(u, object = ghyp::ghyp(lambda = -0.5*nu, chi = nu, psi = 0, mu = -nu/(nu-2)*gamma,
                                             sigma = solve(A) %*% diag(sigma^2, nrow = K) %*% t(solve(A)),
                                             gamma = gamma), logvalue = T)
  }
  return(sum(allw))
}

#' @export
int_w_MTnonSV <- function(y, xt, A, B, sigma, nu, t_max, K, R = 100){
  u <- (t(y) - B %*%xt)
  u_proposal <- A %*% (u)

  # A.tmp <- diag(1/sigma) %*% A
  a_target <- (nu*0.5 + 1*0.5)
  b_target <- (nu*0.5 + 0.5 * u_proposal^2 / sigma^2)

  store_llw <- array(NA, dim = c(K,t_max,R))
  for (i in c(1:R)){

    w <- matrix(rinvgamma(K * t_max, shape = a_target, rate = b_target), nrow = K)
    dens_w <- dinvgamma(x = w, shape = a_target, rate = b_target, log = T)
    Prior_dens_w <- dinvgamma(x = w, shape = nu*0.5, rate = nu*0.5, log = T)

    w_sqrt <- sqrt(w)
    u_proposal <- A %*% ((u/w_sqrt)) / sigma
    store_llw[,,i] <- - 0.5 * log(2*pi) - 0.5 * u_proposal^2 -
      0.5 * log(w) - log(sigma) + Prior_dens_w - dens_w
  }
  maxllike = apply(store_llw, MARGIN = c(1,2), max)
  llw = log( apply(exp(store_llw-c(maxllike)), MARGIN = c(1,2), FUN = mean) ) + maxllike
  sum(llw)
  return(sum(llw))
}

#' @export
int_w_MSTnonSV <- function(y, xt, A, B, sigma, nu, gamma = rep(0,K), t_max, K, R = 100){
  u <- t(y) - B %*%xt + nu/(nu-2)*gamma
  u_proposal <- A %*% (t(y) - B %*%xt) + nu/(nu-2)*gamma

  a_target <- (nu*0.5 + 1*0.5)
  b_target <- (nu*0.5 + 0.5 * u_proposal^2 / sigma^2)

  # lambda = -(nu+1)*0.5
  # chi = nu + (u_proposal^2)/ sigma^2
  # psi = gamma^2/sigma^2

  store_llw <- array(NA, dim = c(K,t_max,R))
  for (i in c(1:R)){

    w <- matrix(rinvgamma(K * t_max, shape = a_target, rate = b_target), nrow = K)
    dens_w <- dinvgamma(x = w, shape = a_target, rate = b_target, log = T)
    Prior_dens_w <- dinvgamma(x = w, shape = nu*0.5, rate = nu*0.5, log = T)

    w_sqrt <- sqrt(w)
    u_proposal <- A %*% ((u/w_sqrt - w_sqrt * gamma)) / sigma

    store_llw[,,i] <- - 0.5 * log(2*pi) - 0.5 * u_proposal^2 -
                        0.5 * log(w) - log(sigma) + Prior_dens_w - dens_w
  }
  maxllike = apply(store_llw, MARGIN = c(1,2), max)
  llw = log( apply(exp(store_llw-c(maxllike)), MARGIN = c(1,2), FUN = mean) ) + maxllike
  sum(llw)
  return(sum(llw))
}

#' @export
int_w_OSTnonSV <- function(y, xt, A, B, sigma, nu, gamma = rep(0,K), t_max, K, R = 100){
  u <- A %*% (t(y) - B %*%xt) + nu/(nu-2) * gamma
  llw <- rep(0,K)
  for (j in c(1:K)){
    llw[j] <- sum(ghyp::dghyp(u[j,], object = ghyp::ghyp(lambda = -0.5*nu[j], chi = nu[j], psi = 0,
                                                         mu = 0, sigma = sigma[j],
                                                        gamma = gamma[j]), logvalue = TRUE)) # input here sigma, not sigma^2
  }
  return(sum(llw))
}


# This function evaluates the integrated likelihood of the VAR-SV model in
# Chan and Eisenstat (2018)
#' @export
int_h_Gaussian <- function(ytilde, h0, sigma_h, t_max, K, R = 100){
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

    c_pri = -t_max*K*0.5*log(2*pi) -.5*t_max*sum(log(sigma_h))
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
int_h_StudentSV <- function(yt, xt, B, A, h0, sigma_h, nu, gamma = rep(0,K), h_mean, w_mean, t_max, K, R = 100){
  y_tilde = A %*% ( (yt - B %*% xt) )  # Student dist with Identity correlation matrix
  s2 = as.numeric(y_tilde)^2
  max_loop = 10000
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


  c_pri = -t_max*K*0.5*log(2*pi) -.5*t_max*sum(log(sigma_h))
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
int_h_SkewSV <- function(yt, xt, B, A, h0, sigma_h, nu, gamma = rep(0,K), h_mean, w_mean, t_max, K, R = 100){
  u <- yt - B %*% xt + nu/(nu-2) * gamma
  y_tilde = A %*% (( u - reprow(w_mean,K)*gamma ) ) # Approximation multivariate Student dist
  s2 = as.numeric(y_tilde)^2
  max_loop = 10000
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

  c_pri = -t_max*K*0.5*log(2*pi) -.5*t_max*sum(log(sigma_h))
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
int_h_MTSV <- function(yt, xt, B, A, h0, sigma_h, nu, h_mean, w_mean, t_max, K, R = 100){
  u <- (yt - B %*% xt)
  y_tilde <- (A %*% u) # Mimic univariate Student
  s2 = as.numeric(y_tilde)^2
  max_loop = 10000
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

  c_pri = -t_max*K*0.5*log(2*pi) -.5*t_max*sum(log(sigma_h))
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
int_h_MSTSV <- function(yt, xt, B, A, h0, sigma_h, nu, gamma = rep(0,K), h_mean, w_mean, t_max, K, R = 100){
  u <- (yt - B %*% xt - w_mean * gamma + nu/(nu-2)*gamma)
  y_tilde <- (A %*% u) # Mimic Univariate Student
  s2 = as.numeric(y_tilde)^2
  max_loop = 10000
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

  c_pri = -t_max*K*0.5*log(2*pi) -.5*t_max*sum(log(sigma_h))
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
int_h_OTSV <- function(yt, xt, B, A, h0, sigma_h, nu, h_mean, w_mean, t_max, K, R = 100){
  u <- (yt - B %*% xt)
  y_tilde <- (A %*% u) # Univarite Student
  s2 = as.numeric(y_tilde)^2
  max_loop = 10000
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


  c_pri = -t_max*K*0.5*log(2*pi) -.5*t_max*sum(log(sigma_h))
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
int_h_OSTSV <- function(yt, xt, B, A, h0, sigma_h, nu, gamma = rep(0,K), h_mean, w_mean, t_max, K, R = 100){
  u <- (yt - B %*% xt )
  y_tilde <- (A %*% u - w_mean * gamma + nu/(nu-2)*gamma) # Mimic univariate Student distribution
  s2 = as.numeric(y_tilde)^2
  max_loop = 10000
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

  c_pri = -t_max*K*0.5*log(2*pi) -.5*t_max*sum(log(sigma_h))
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

#' @export
Normal_approx <- function(mcmc_sample, ndraws){
  mcmc_mean <- apply(mcmc_sample, 2, mean)
  mcmc_Sigma <- cov(mcmc_sample)
  nElements <- length(mcmc_mean)

  new_samples <- mvnfast::rmvn(ndraws, mu = mcmc_mean, sigma = mcmc_Sigma)
  colnames(new_samples) <- colnames(mcmc_sample)
  sum_log_prop <- mvnfast::dmvn(X = new_samples,
                                mu = as.numeric(mcmc_mean), sigma = mcmc_Sigma, log = T)
  return(list(new_samples = new_samples,
              sum_log_prop = sum_log_prop))
}


#' @export
Gamma_approx <- function(mcmc_sample, ndraws){

  # MASS::fitdistr(1/Sigma_mat[,1], "gamma")
  # MASS::fitdistr(1/Sigma_mat[,3]/1000, "gamma")
  # MASS::fitdistr(Sigma_mat[,3], "exponential")

  nElements <- ncol(mcmc_sample)
  new_samples <- matrix(NA, ncol = nElements, nrow = ndraws, dimnames = list(c(), colnames(mcmc_sample)))
  Density_prop <-  matrix(NA, ncol = nElements, nrow = ndraws)
  shape_param <- rep(0, nElements)
  rate_param <- rep(0, nElements)
  mcmc_mean <- apply(mcmc_sample, 2, mean)
  mcmc_sd <- apply(mcmc_sample, 2, sd)

  for (i in c(1:nElements)){
    if (mcmc_sd[i] > 0){
        fit.gamma <- tryCatch({
          fitdistrplus::fitdist(as.numeric(mcmc_sample[,i]), distr = "gamma", method = "mle")
        }, error = function(e) {
          fitdistrplus::fitdist(as.numeric(mcmc_sample[,i]), distr = "exp", method = "mle")
        })
      if (fit.gamma$distname == "gamma"){
        shape_param[i] <- fit.gamma$estimate[1]
        rate_param[i] <- fit.gamma$estimate[2]
        new_samples[,i] <- rgamma(ndraws, shape = shape_param[i], rate_param[i])
        Density_prop[,i] <- dgamma(new_samples[,i],
                                                shape = shape_param[i],
                                                rate = rate_param[i], log = T)
      }
      if (fit.gamma$distname == "exp"){
        rate_param[i] <- fit.gamma$estimate[1]
        new_samples[,i] <- rexp(ndraws, rate =  rate_param[i])
        Density_prop[,i] <- dexp(new_samples[,i], rate = rate_param[i], log = T)

      }

    } else {
      new_samples[,i] <- as.numeric(mcmc_sample[,i])
      Density_prop[,i] <- 0
    }

  }

  return(list(new_samples = new_samples,
              sum_log_prop = apply(Density_prop, 1, sum)))
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

  for (i in c(1:nElements)){
      fit.gamma <- tryCatch({
        fitdistrplus::fitdist(1/as.numeric(mcmc_sample[,i]), distr = "gamma", method = "mle")
      }, error = function(e) {
        fitdistrplus::fitdist(1/as.numeric(mcmc_sample[,i]), distr = "gamma", method = "qme", probs = c(1/3, 2/3))
      })
      # cat(i, " = Invgamma ")
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
Nu_Gamma_approx <- function(mcmc_sample, ndraws, dist){
  nElements <- ncol(mcmc_sample)
  new_samples <- matrix(NA, ncol = nElements, nrow = ndraws, dimnames = list(c(), colnames(mcmc_sample)))
  Density_prop <-  matrix(NA, ncol = nElements, nrow = ndraws)
  shape_param <- rep(0, nElements)
  rate_param <- rep(0, nElements)
  for (i in c(1:nElements)){
    fit.gamma <- fitdistrplus::fitdist(as.numeric(mcmc_sample[,i]), distr = "gamma", method = "mle")
    shape_param[i] <- fit.gamma$estimate[1]
    rate_param[i] <- fit.gamma$estimate[2]
    tmp <- rgamma(ndraws*2, shape = shape_param[i], rate_param[i])
    tmp <- tmp[tmp > 4]
    tmp <- tmp[tmp < 100]
    new_samples[,i] <- head(tmp,ndraws)
    Density_prop[,i] <- dgamma(new_samples[,i], shape = shape_param[i], rate = rate_param[i], log = T)
  }
  return(list(new_samples = new_samples,
              sum_log_prop = apply(Density_prop, 1, sum)))
}

#' @export
LL_tnorm <- function(par, data){
  - sum(log(dtruncnorm(x = data, a=0, b=Inf, mean = par[1], sd = par[2])))
}

#' @export
Normal_trunc_approx <- function(mcmc_sample, ndraws){

  nElements <- ncol(mcmc_sample)
  mcmc_mean <- rep(0, nElements)
  mcmc_Sigma <- rep(0, nElements)
  new_samples <- matrix(NA, ncol = nElements, nrow = ndraws)
  sum_log_prop <- rep(0, ndraws)

  for (i in c(1:nElements)){
    result <- optim(par = c(0, 1), fn = LL_tnorm, data = mcmc_sample[,i])
    new_samples[,i] <- rtruncnorm(ndraws, mean = result$par[1], sd = result$par[2], a = 0)
    sum_log_prop <- sum_log_prop +
      log(dtruncnorm(x = new_samples[,i], a=0, b=Inf, mean = result$par[1], sd = result$par[2])) - log(2) - log(new_samples[,i]) #Jacobian trans
  }

  colnames(new_samples) <- colnames(mcmc_sample)
  return(list(new_samples = new_samples,
              sum_log_prop = sum_log_prop))
}

#' @export
a0toA <- function(a0, K){
  A <- matrix(0, nrow = K, ncol = K)
  A[upper.tri(A)] <- a0
  A <- t(A)
  diag(A) <- 1
  return(A)
}

