#' @export
Chib.LML.nonSV <- function(Chain,
                           samples = 3000, burnin = 0, thin = 1, numCores = 4){
  # samples = 3000; burnin = 0; thin = 1;
  # numCores <- 16;

  P1 <- P2 <- P3 <- P4 <- P5 <- 0

  prior <- Chain$prior
  inits <- Chain$inits
  dist <- prior$dist

  B_mat <- get_post(Chain, element = "B")
  A_mat <- get_post(Chain, element = "a")
  Gamma_mat <- get_post(Chain, element = "gamma")
  Sigma_mat <- get_post(Chain, element = "sigma")
  Nu_mat <- get_post(Chain, element = "nu")
  W_mat <- get_post(Chain, element = "w")
  H_mat <- get_post(Chain, element = "h")

  B_med <- apply(B_mat, MARGIN = 2, FUN = mean)
  A_med <- apply(A_mat, MARGIN = 2, FUN = mean)
  Gamma_med <- apply(Gamma_mat, MARGIN = 2, FUN = mean)
  Sigma_med <- apply(Sigma_mat, MARGIN = 2, FUN = mean)
  Nu_med <- apply(Nu_mat, MARGIN = 2, FUN = mean)
  W_med <- apply(W_mat, MARGIN = 2, FUN = mean)
  H_med <- apply(H_mat, MARGIN = 2, FUN = mean)

  inits$samples <- samples # Reset
  inits$burnin <- burnin # Reset
  inits$thin <- thin # Reset

  # Commom parameters
  inits$A0 <- a0toA(A_med, K = K)
  inits$B0 <- matrix(B_med, nrow = K)
  inits$sigma <- Sigma_med

  str(inits)

  if (prior$SV == FALSE){
    Start = Sys.time()
    if (dist == "Gaussian") {

      cat("Pi( B* | y, A*, Sigma*) \n ")
      ChibLLP_chain <- ChibLLP(Chain, ndraws = 100, numCores = numCores)
      RhpcBLASctl::blas_set_num_threads(RhpcBLASctl::get_num_cores())

      P1 <- Chib.Gaussian.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                                samples = samples, burnin = burnin, thin = thin,
                                fix_B = FALSE, fix_A = TRUE, fix_Sigma = TRUE,
                                cal_B = TRUE, cal_A = FALSE, cal_Sigma = FALSE)
      cat("Pi( A* | y, Sigma*) \n ")
      P2 <- Chib.Gaussian.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                                samples = samples, burnin = burnin, thin = thin,
                                fix_B = FALSE, fix_A = FALSE, fix_Sigma = TRUE,
                                cal_B = FALSE, cal_A = TRUE, cal_Sigma = FALSE)
      cat("Pi( Sigma* | y) \n ")
      P3 <- Chib.Gaussian.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                                samples = samples, burnin = burnin, thin = thin,
                                fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE,
                                cal_B = FALSE, cal_A = FALSE, cal_Sigma = TRUE)

      CML_Chain <- ChibLLP_chain$LL - logmeanexp(P1) - logmeanexp(P2) - logmeanexp(P3)
    } else {
      inits$nu <- Nu_med
      inits$logsigma_nu <- log(apply(Nu_mat, MARGIN = 2, sd))
      inits$gamma <- Gamma_med
    }
    if (dist == "Student") {
      inits$w <- W_med
      ChibLLP_chain <- ChibLLP(Chain, ndraws = 100, numCores = numCores)
      RhpcBLASctl::blas_set_num_threads(RhpcBLASctl::get_num_cores())

      cat("Pi( B* | y, A*, Sigma*, Nu*) \n ")
      P1 <- Chib.Student.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                               samples = samples, burnin = burnin, thin = thin,
                               fix_B = FALSE, fix_A = TRUE, fix_Sigma = TRUE, fix_Nu = TRUE,
                               cal_B = TRUE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE)


      cat("Pi( A* | y, Sigma*, Nu*) \n ")
      P2 <- Chib.Student.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                               samples = samples, burnin = burnin, thin = thin,
                               fix_B = FALSE, fix_A = FALSE, fix_Sigma = TRUE, fix_Nu = TRUE,
                               cal_B = FALSE, cal_A = TRUE, cal_Sigma = FALSE, cal_Nu = FALSE)
      cat(" Pi( Sigma* | y, Nu*) \n ")
      P3 <- Chib.Student.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                               samples = samples, burnin = burnin, thin = thin,
                               fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                               cal_B = FALSE, cal_A = FALSE, cal_Sigma = TRUE, cal_Nu = FALSE)

      cat(" Pi( Nu* | y) --> Nominator \n ")
      P4 <- Chib.Student.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                               samples = samples, burnin = burnin, thin = thin,
                               fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                               cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      cat("Pi( Nu* | y) --> Denominator \n ")
      P5 <- Chib.Student.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                               samples = samples, burnin = burnin, thin = thin,
                               fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                               cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      CML_Chain <- ChibLLP_chain$LL - logmeanexp(P1) - logmeanexp(P2) - logmeanexp(P3) - logmeanexp(P4) + logmeanexp(P5) # P5 Denominator

    }
    if (dist == "Skew.Student") {
      inits$w <- W_med
      ChibLLP_chain <- ChibLLP(Chain, ndraws = 100, numCores = numCores)
      RhpcBLASctl::blas_set_num_threads(RhpcBLASctl::get_num_cores())

      cat("Pi( B* | y, A*, Sigma*, Nu*) \n ")
      P1 <- Chib.Skew.Student.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                                    samples = samples, burnin = burnin, thin = thin,
                                    fix_B = FALSE, fix_A = TRUE, fix_Sigma = TRUE, fix_Nu = TRUE,
                                    cal_B = TRUE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE)


      cat("Pi( A* | y, Sigma*, Nu*) \n ")
      P2 <- Chib.Skew.Student.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                                    samples = samples, burnin = burnin, thin = thin,
                                    fix_B = FALSE, fix_A = FALSE, fix_Sigma = TRUE, fix_Nu = TRUE,
                                    cal_B = FALSE, cal_A = TRUE, cal_Sigma = FALSE, cal_Nu = FALSE)
      cat(" Pi( Sigma* | y, Nu*) \n ")
      P3 <- Chib.Skew.Student.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                                    samples = samples, burnin = burnin, thin = thin,
                                    fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                                    cal_B = FALSE, cal_A = FALSE, cal_Sigma = TRUE, cal_Nu = FALSE)

      cat(" Pi( Nu* | y) --> Nominator \n ")
      P4 <- Chib.Skew.Student.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                                    samples = samples, burnin = burnin, thin = thin,
                                    fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                                    cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      cat("Pi( Nu* | y) --> Denominator \n ")
      P5 <- Chib.Skew.Student.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                                    samples = samples, burnin = burnin, thin = thin,
                                    fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                                    cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      CML_Chain <- ChibLLP_chain$LL - logmeanexp(P1) - logmeanexp(P2) - logmeanexp(P3) - logmeanexp(P4) + logmeanexp(P5) # P5 Denominator


    }
    if (dist == "MT") {
      inits$w <- matrix(W_med, nrow = K)
      ChibLLP_chain <- ChibLLP(Chain, ndraws = 1000, numCores = numCores)
      RhpcBLASctl::blas_set_num_threads(RhpcBLASctl::get_num_cores())


      cat("Pi( B* | y, A*, Sigma*, Nu*) \n ")
      P1 <- Chib.MT.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = TRUE, fix_Sigma = TRUE, fix_Nu = TRUE,
                          cal_B = TRUE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE)


      cat("Pi( A* | y, Sigma*, Nu*) \n ")
      P2 <- Chib.MT.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = TRUE, fix_Nu = TRUE,
                          cal_B = FALSE, cal_A = TRUE, cal_Sigma = FALSE, cal_Nu = FALSE)
      cat("Pi( Sigma* | y, Nu*) \n ")
      P3 <- Chib.MT.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                          cal_B = FALSE, cal_A = FALSE, cal_Sigma = TRUE, cal_Nu = FALSE)

      cat("Pi( Nu* | y) --> Nominator \n ")#
      P4 <- Chib.MT.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                          cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      cat("Pi( Nu* | y) --> Denominator \n ")
      P5 <- Chib.MT.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                          cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      CML_Chain <- ChibLLP_chain$LL - logmeanexp(P1) - logmeanexp(P2) - logmeanexp(P3) - sum(logmeanexp(P4)) + sum(logmeanexp(P5) ) # P5 Denominator

    }
    if (dist == "MST") {
      inits$w <- matrix(W_med, nrow = K)
      ChibLLP_chain <- ChibLLP(Chain, ndraws = 1000, numCores = numCores)
      RhpcBLASctl::blas_set_num_threads(RhpcBLASctl::get_num_cores())


      cat("Pi( B* | y, A*, Sigma*, Nu*) \n ")
      P1 <- Chib.MST.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = TRUE, fix_Sigma = TRUE, fix_Nu = TRUE,
                          cal_B = TRUE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE)


      cat("Pi( A* | y, Sigma*, Nu*) \n ")
      P2 <- Chib.MST.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = TRUE, fix_Nu = TRUE,
                          cal_B = FALSE, cal_A = TRUE, cal_Sigma = FALSE, cal_Nu = FALSE)
      cat("Pi( Sigma* | y, Nu*) \n ")
      P3 <- Chib.MST.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                          cal_B = FALSE, cal_A = FALSE, cal_Sigma = TRUE, cal_Nu = FALSE)

      cat("Pi( Nu* | y) --> Nominator \n ")#
      P4 <- Chib.MST.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                          cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      cat("Pi( Nu* | y) --> Denominator \n ")
      P5 <- Chib.MST.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                          cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      CML_Chain <- ChibLLP_chain$LL - logmeanexp(P1) - logmeanexp(P2) - logmeanexp(P3) - sum(logmeanexp(P4)) + sum(logmeanexp(P5) ) # P5 Denominator

    }
    if (dist == "OT") {
      inits$w <- matrix(W_med, nrow = K)
      ChibLLP_chain <- ChibLLP(Chain, ndraws = 1000, numCores = numCores)
      RhpcBLASctl::blas_set_num_threads(RhpcBLASctl::get_num_cores())


      cat("Pi( B* | y, A*, Sigma*, Nu*) \n ")
      P1 <- Chib.OT.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                           samples = samples, burnin = burnin, thin = thin,
                           fix_B = FALSE, fix_A = TRUE, fix_Sigma = TRUE, fix_Nu = TRUE,
                           cal_B = TRUE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE)


      cat("Pi( A* | y, Sigma*, Nu*) \n ")
      P2 <- Chib.OT.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                           samples = samples, burnin = burnin, thin = thin,
                           fix_B = FALSE, fix_A = FALSE, fix_Sigma = TRUE, fix_Nu = TRUE,
                           cal_B = FALSE, cal_A = TRUE, cal_Sigma = FALSE, cal_Nu = FALSE)
      cat("Pi( Sigma* | y, Nu*) \n ")
      P3 <- Chib.OT.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                           samples = samples, burnin = burnin, thin = thin,
                           fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                           cal_B = FALSE, cal_A = FALSE, cal_Sigma = TRUE, cal_Nu = FALSE)

      cat("Pi( Nu* | y) --> Nominator \n ")#
      P4 <- Chib.OT.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                           samples = samples, burnin = burnin, thin = thin,
                           fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                           cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      cat("Pi( Nu* | y) --> Denominator \n ")
      P5 <- Chib.OT.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                           samples = samples, burnin = burnin, thin = thin,
                           fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                           cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      CML_Chain <- ChibLLP_chain$LL - logmeanexp(P1) - logmeanexp(P2) - logmeanexp(P3) - sum(logmeanexp(P4)) + sum(logmeanexp(P5) ) # P5 Denominator

    }
    if (dist == "OST") {
      inits$w <- matrix(W_med, nrow = K)
      ChibLLP_chain <- ChibLLP(Chain, ndraws = 1000, numCores = numCores)
      RhpcBLASctl::blas_set_num_threads(RhpcBLASctl::get_num_cores())


      cat("Pi( B* | y, A*, Sigma*, Nu*) \n ")
      P1 <- Chib.OST.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = TRUE, fix_Sigma = TRUE, fix_Nu = TRUE,
                          cal_B = TRUE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE)


      cat("Pi( A* | y, Sigma*, Nu*) \n ")
      P2 <- Chib.OST.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = TRUE, fix_Nu = TRUE,
                          cal_B = FALSE, cal_A = TRUE, cal_Sigma = FALSE, cal_Nu = FALSE)
      cat("Pi( Sigma* | y, Nu*) \n ")
      P3 <- Chib.OST.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                          cal_B = FALSE, cal_A = FALSE, cal_Sigma = TRUE, cal_Nu = FALSE)

      cat("Pi( Nu* | y) --> Nominator \n ")#
      P4 <- Chib.OST.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                          cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      cat("Pi( Nu* | y) --> Denominator \n ")
      P5 <- Chib.OST.novol(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                          cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      CML_Chain <- ChibLLP_chain$LL - logmeanexp(P1) - logmeanexp(P2) - logmeanexp(P3) - sum(logmeanexp(P4)) + sum(logmeanexp(P5) ) # P5 Denominator

    }

    elapsedTime = Sys.time() - Start
    print(elapsedTime)
    out <- list(CML = CML_Chain,
                LLP = ChibLLP_chain$LL,
                lprior = ChibLLP_chain$lprior,
                lc = ChibLLP_chain$lc,
                LLP_chain = ChibLLP_chain,
                aP1 = logmeanexp(P1), aP2 = logmeanexp(P2), aP3 = logmeanexp(P3),
                aP4 = sum(logmeanexp(P4)), aP5 = sum(logmeanexp(P5)),
                P1 = P1, P2 = P2, P3 = P3, P4 = P4, P5 = P5,
                ChibLLP_chain = ChibLLP_chain,
                esttime = elapsedTime)
    class(out) <- c("ChibMLM")

  } else {
    warning("prior$SV is TRUE")
  }
  return(out)
}

###########################################################################
#' @export
Chib.Gaussian.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL,
                                samples = 1000, burnin = 0, thin = 1,
                                fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE,
                                cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)

  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Gaussian", SV = FALSE)
  }
  # prior B
  b_prior = prior$b_prior
  V_b_prior = prior$V_b_prior
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior A
  a_prior = prior$a_prior
  V_a_prior = prior$V_a_prior


  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  #samples <- inits$samples
  # burnin <- inits$burnin
  # thin <- inits$thin

  A <- inits$A0
  B <- inits$B0
  b_sample <- as.numeric(B)
  a_sample <- as.numeric(A[lower.tri(A)])
  sigma <- inits$sigma
  # Sigma <- solve(inits$A0) %*% diag(sigma, nrow = length(sigma))
  # Sigma2 <- Sigma %*% t(Sigma)
  # Sigma2_inv <- solve(Sigma2)
  V_b_prior_inv <- solve(V_b_prior)
  # new precompute
  theta.prior.precmean <- V_b_prior_inv %*% b_prior

  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K,
                 ncol = (samples - burnin)%/% thin)
  lpost <- rep(0, (samples - burnin)%/% thin)
  for (j in c(1:samples)){

    # sample B
    if(!fix_B) {

      A.tmp <- diag(1/sigma) %*% A
      y.tilde <- as.vector( A.tmp %*% yt )
      x.tilde <- kronecker( t(xt), A.tmp )
      theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
      b_sample <- backsolve( theta.prec.chol,
                             backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                        upper.tri = T, transpose = T )
                           + rnorm(K*m) )
      B <- matrix(b_sample,K,m)
    }
    if (cal_B & (j > burnin) & (j %% thin == 0)){
      B_star <- vec(inits$B)
      A.tmp <- diag(1/sigma) %*% A
      y.tilde <- as.vector( A.tmp %*% yt )
      x.tilde <- kronecker( t(xt), A.tmp )

      theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
      theta.prior.precmean <- V_b_prior_inv %*% prior$b_prior
      b_star <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T ))
      #mvnfast::dmvn(B_star, mu = b_star, sigma = (solve(V_b_prior_inv + crossprod(x.tilde) )), log = T)
      lpost[(j - burnin) %/% thin] <- - length(B_star) * 0.5 * log(2*pi) + sum(log(diag(theta.prec.chol))) -
                                                  0.5 * t(B_star - b_star) %*% (V_b_prior_inv + crossprod(x.tilde)) %*% (B_star - b_star)
    }

    # Sample sigma
    if(!fix_Sigma) {

      sigma2 <- rep(0,K)
      u_ort <- A %*% (yt - B %*%xt)
      sigma_post_a <- sigma0_T0 + rep(t_max,K)
      sigma_post_b <- sigma0_S0 + rowSums(u_ort^2)
      for (i in c(1:K)){
        sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
      }
      sigma <- sqrt(sigma2)
    }

    if (cal_Sigma & (j > burnin) & (j %% thin == 0)){
        u_ort <- A %*% (yt - B %*%xt)
        sigma_post_a <- prior$sigma_T0 + rep(t_max,K)
        sigma_post_b <- prior$sigma_S0 + rowSums(u_ort^2)
        lpost[(j - burnin) %/% thin] <- sum(dinvgamma(inits$sigma^2, shape = sigma_post_a * 0.5, rate = sigma_post_b * 0.5, log = T))
      }


    # Sample A0
    if(!fix_A) {

      u_std <- (yt - B %*%xt)
      u_neg <- - u_std
      a_sample <- rep(0, K * (K - 1) /2)
      if (K > 1) {
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- a_prior[id_start:id_end]
          V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
          a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sigma[i],
                                                       xsub = matrix(u_neg[1:(i-1),] / sigma[i], nrow = i-1),
                                                       a_sub = a_sub,
                                                       V_a_sub = V_a_sub)
        }
      }
      A_post <- matrix(0, nrow = K, ncol = K)
      A_post[upper.tri(A)] <- a_sample
      A <- t(A_post)
      diag(A) <- 1
    }

    if (cal_A & (j > burnin) & (j %% thin == 0)){
      u_std <- (yt - B %*%xt)
      u_neg <- - u_std
      A_star <- inits$A0[lower.tri(A)]
      if (K > 1) {
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- prior$a_prior[id_start:id_end]
          V_a_sub <- prior$V_a_prior[id_start:id_end, id_start:id_end]
          sub_A_star <- A_star[id_start:id_end]

          ysub = u_std[i,] / sigma[i]
          xsub = matrix(u_neg[1:(i-1),] / sigma[i], nrow = i-1)

          n = nrow(xsub)

          V_a_sub_inv = solve(V_a_sub)

          V_a_post_inv <- V_a_sub_inv + xsub %*% t(xsub)
          V_a_post <- solve(V_a_post_inv)
          a_post <- V_a_post %*% ( V_a_sub_inv %*% a_sub + xsub %*% ysub)

          lpost[(j - burnin) %/% thin] <- lpost[(j - burnin) %/% thin] + mvnfast::dmvn(sub_A_star, mu = a_post, sigma = V_a_post, log = T)
        }
      }
    }

    if ((j > burnin) & (j %% thin == 0))
      mcmc[, (j - burnin) %/% thin] <- c(b_sample, a_sample, sigma)
    if (j %% 1000 == 0) { cat(" Iteration ", j, " \n")}
  }
  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("sigma",c(1:K), sep = ""))
  return(lpost)
}

###########################################################################
#' @export
Chib.Student.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL,
                               samples = 1000, burnin = 0, thin = 1,
                               fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                               cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)

  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Student", SV = FALSE)
  }
  # prior B
  b_prior = prior$b_prior
  V_b_prior = prior$V_b_prior
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior A
  a_prior = prior$a_prior
  V_a_prior = prior$V_a_prior
  # prior nu
  nu_gam_a = prior$nu_gam_a
  nu_gam_b = prior$nu_gam_b

  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }

  #samples <- inits$samples
  # burnin <- inits$burnin
  # thin <- inits$thin
  A <- inits$A0
  B <- Vec_to_Mat(inits$B0, K, p)
  b_sample <- as.numeric(B)
  a_sample <- as.numeric(A[lower.tri(A)])
  sigma <- inits$sigma
  # Sigma <- solve(inits$A0) %*% diag(inits$sigma, nrow = K)
  # Sigma2 <- Sigma %*% t(Sigma)
  # Sigma2_inv <- solve(Sigma2)
  V_b_prior_inv <- solve(V_b_prior)

  nu <- inits$nu
  logsigma_nu <- inits$logsigma_nu
  acount_nu <- 0
  acount_w <- rep(0, t_max)
  # Init w as Gaussian
  w_sample <- inits$w
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)

  # new precompute
  theta.prior.precmean <- V_b_prior_inv %*% b_prior

  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + 1 + t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  lpost <- rep(0, (samples - burnin)%/% thin)
  for (j in c(1:samples)){

    # sample B
    if(!fix_B) {
      A.tmp <- diag(1/sigma) %*% A
      wt <- as.vector(1/w_sqrt)
      y.tilde <- as.vector( A.tmp %*% yt ) * wt
      x.tilde <- kronecker( t(xt), A.tmp ) * wt
      theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
      b_sample <- backsolve( theta.prec.chol,
                             backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                        upper.tri = T, transpose = T )
                             + rnorm(K*m) )
      B <- matrix(b_sample,K,m)
    }

    if (cal_B & (j > burnin) & (j %% thin == 0)){
      B_star <- vec(inits$B)
      A.tmp <- diag(1/sigma) %*% A
      wt <- as.vector(1/w_sqrt)
      y.tilde <- as.vector( A.tmp %*% yt ) * wt
      x.tilde <- kronecker( t(xt), A.tmp ) * wt
      theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )

      theta.prior.precmean <- V_b_prior_inv %*% prior$b_prior
      b_star <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T ))
      #mvnfast::dmvn(B_star, mu = b_star, sigma = (solve(V_b_prior_inv + crossprod(x.tilde) )), log = T)
      lpost[(j - burnin) %/% thin] <- - length(B_star) * 0.5 * log(2*pi) + sum(log(diag(theta.prec.chol))) -
        0.5 * t(B_star - b_star) %*% (V_b_prior_inv + crossprod(x.tilde)) %*% (B_star - b_star)
    }


    # Sample sigma
    if(!fix_Sigma) {
    sigma2 <- rep(0,K)
    u_ort <- A %*% ((yt - B %*%xt)/ w_sqrt)    # change from Gaussian
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
    sigma_post_b <- sigma0_S0 + rowSums(u_ort^2)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)
    }

    if (cal_Sigma & (j > burnin) & (j %% thin == 0)){
      u_ort <- A %*% ((yt - B %*%xt)/ w_sqrt)    # change from Gaussian
      sigma_post_a <- prior$sigma_T0 + rep(t_max,K)
      sigma_post_b <- prior$sigma_S0 + rowSums(u_ort^2)
      lpost[(j - burnin) %/% thin] <- sum(dinvgamma(inits$sigma^2, shape = sigma_post_a * 0.5, rate = sigma_post_b * 0.5, log = T))
    }

    # Sample A0
    if(!fix_A) {
    u_std <- (yt - B %*%xt) / w_sqrt # change from Gaussian
    u_neg <- - u_std
    a_sample <- rep(0, K * (K - 1) /2)
    if (K > 1) {
      for (i in c(2:K)){
        id_end <- i*(i-1)/2
        id_start <- id_end - i + 2
        a_sub <- a_prior[id_start:id_end]
        V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
        a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sigma[i] ,
                                                     xsub = matrix(u_neg[1:(i-1),] / sigma[i], nrow = i-1),
                                                     a_sub = a_sub,
                                                     V_a_sub = V_a_sub)
      }
    }
    A_post <- matrix(0, nrow = K, ncol = K)
    A_post[upper.tri(A)] <- a_sample
    A <- t(A_post)
    diag(A) <- 1
    }

    if (cal_A & (j > burnin) & (j %% thin == 0)){
      u_std <- (yt - B %*%xt) / w_sqrt # change from Gaussian
      u_neg <- - u_std
      A_star <- inits$A0[lower.tri(A)]
      if (K > 1) {
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- prior$a_prior[id_start:id_end]
          V_a_sub <- prior$V_a_prior[id_start:id_end, id_start:id_end]
          sub_A_star <- A_star[id_start:id_end]

          ysub = u_std[i,] / sigma[i]
          xsub = matrix(u_neg[1:(i-1),] / sigma[i], nrow = i-1)

          n = nrow(xsub)

          V_a_sub_inv = solve(V_a_sub)

          V_a_post_inv <- V_a_sub_inv + xsub %*% t(xsub)
          V_a_post <- solve(V_a_post_inv)
          a_post <- V_a_post %*% ( V_a_sub_inv %*% a_sub + xsub %*% ysub)

          lpost[(j - burnin) %/% thin] <- lpost[(j - burnin) %/% thin] + mvnfast::dmvn(sub_A_star, mu = a_post, sigma = V_a_post, log = T)
        }
      }
    }
    # Sample w
    u <- (yt - B %*%xt)
    shape_w <- nu*0.5 + K*0.5
    rate_w <- as.numeric(nu*0.5 + 0.5 * colSums((A %*% u / sigma)^2))
    w_sample <- rinvgamma( n=t_max, shape = shape_w, rate = rate_w )
    w <- reprow(w_sample, K)
    w_sqrt <- sqrt(w)


    # Sample nu
    if(!fix_Nu) {
      nu_temp = nu + exp(logsigma_nu)*rnorm(1)
      if (nu_temp > 4 && nu_temp < 100){
        num_mh = dgamma(nu_temp, shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w_sample, shape = nu_temp*0.5, rate = nu_temp*0.5, log = T))
        denum_mh = dgamma(nu, shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w_sample, shape = nu*0.5, rate = nu*0.5, log = T))
        alpha = num_mh - denum_mh;
        temp = log(runif(1));
        if (alpha > temp){
          nu = nu_temp
          acount_nu = acount_nu + 1
        }

      }
    }
    if (cal_Nu & (j > burnin) & (j %% thin == 0)){
      # Denominator
      if (fix_Nu) {
        nu_temp = inits$nu + exp(logsigma_nu)*rnorm(1)
        if (nu_temp > 4 && nu_temp < 100){
          num_mh = dgamma(nu_temp, shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w_sample, shape = nu_temp*0.5, rate = nu_temp*0.5, log = T))

          denum_mh = dgamma(inits$nu, shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w_sample, shape = inits$nu*0.5, rate = inits$nu*0.5, log = T))
          alpha = min(1, exp(num_mh - denum_mh))
          temp = runif(1)
          if (alpha > temp){
            acount_nu = acount_nu + 1
          }
          lpost[(j - burnin) %/% thin] <- log(alpha) # Take a negative for denominator

        }
      } else {
        # Nominator
        q_nu_nustar <- dnorm(inits$nu, mean = nu, sd = exp(logsigma_nu) ) # Move nu to inits$nu
        num_mh = dgamma(inits$nu, shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w_sample, shape = inits$nu*0.5, rate = inits$nu*0.5, log = T))
        denum_mh = dgamma(nu, shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w_sample, shape = nu*0.5, rate = nu*0.5, log = T))
        alpha = min(1, exp(num_mh - denum_mh))
        lpost[(j - burnin) %/% thin] <- log( alpha * q_nu_nustar)
      }
    }



    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, nu, as.numeric(w_sample))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", round(acount_nu/j,2) ," ", round(nu,2)," \n")
      acount_w <- rep(0,t_max)
    }
    # if(j %% batchlength == 0 ){
    #   if (acount_nu > batchlength * TARGACCEPT){
    #     logsigma_nu = logsigma_nu + adaptamount(j %/% batchlength);
    #   }
    #   if (acount_nu < batchlength * TARGACCEPT){
    #     logsigma_nu = logsigma_nu - adaptamount(j %/% batchlength);
    #   }
    #   acount_nu = 0
    # }
  }


  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("sigma",c(1:K), sep = ""),
                        paste("nu"),
                        paste("w",c(1:t_max), sep = ""))

  return(lpost)
}


#############################################################################################
#' @export
Chib.Skew.Student.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL,
                                    samples = 1000, burnin = 0, thin = 1,
                                    fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                                    cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)

  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Skew.Student", SV = FALSE)
  }
  # prior B
  b_prior = prior$b_prior
  V_b_prior = prior$V_b_prior
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior A
  a_prior = prior$a_prior
  V_a_prior = prior$V_a_prior
  # prior nu
  nu_gam_a = prior$nu_gam_a
  nu_gam_b = prior$nu_gam_b
  # prior gamma
  gamma_prior = prior$gamma_prior
  V_gamma_prior = prior$V_gamma_prior
  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  #samples <- inits$samples
  # burnin <- inits$burnin
  # thin <- inits$thin
  A <- inits$A0
  B <- Vec_to_Mat(inits$B0, K, p)
  b_sample <- as.numeric(B)
  a_sample <- as.numeric(A[lower.tri(A)])
  sigma <- inits$sigma
  V_b_prior_inv <- solve(V_b_prior)

  nu <- inits$nu
  mu.xi <- nu/( nu - 2 )
  logsigma_nu <- inits$logsigma_nu
  acount_nu <- 0
  acount_w <- rep(0, t_max)
  gamma <- inits$gamma
  D <- diag(gamma, nrow = length(gamma))

  # Init w as Gaussian
  w_sample <- inits$w
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)

  # new precompute
  i <- K*m
  theta.prior.prec = matrix(0,i+K,i+K)
  theta.prior.prec[1:i,1:i] <- V_b_prior_inv
  theta.prior.prec[i+(1:K),i+(1:K)] <- solve( V_gamma_prior )
  theta.prior.precmean <- theta.prior.prec %*% c( b_prior, gamma_prior )

  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + 1 + t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  lpost <- rep(0, (samples - burnin)%/% thin)

  for (j in c(1:samples)){

    # Sample B and gamma
    if(!fix_B) {
      wt <- as.vector(1/(sigma*w_sqrt))
      y.tilde <- as.vector( A %*% yt ) * wt
      x.tilde <- kronecker( cbind( t(xt), w_sample - mu.xi ), A ) * wt
      theta.prec.chol <- chol( theta.prior.prec + crossprod(x.tilde) )
      theta <- backsolve( theta.prec.chol,
                          backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                     upper.tri = T, transpose = T )
                          + rnorm(K*(m+1)) )
      B <- matrix(theta[1:(m*K)],K,m)
      b_sample <- as.vector(B)
      gamma <- theta[(m*K+1):((m+1)*K)]
      D <- diag(gamma)
    }
    if (cal_B & (j > burnin) & (j %% thin == 0)){
      B_star <- c(vec(inits$B), inits$gamma)
      wt <- as.vector(1/(sigma*w_sqrt))
      y.tilde <- as.vector( A %*% yt ) * wt
      x.tilde <- kronecker( cbind( t(xt), w_sample - mu.xi ), A ) * wt
      theta.prec.chol <- chol( theta.prior.prec + crossprod(x.tilde) )
      b_star <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T ))
      #mvnfast::dmvn(B_star, mu = b_star, sigma = (solve(V_b_prior_inv + crossprod(x.tilde) )), log = T)
      lpost[(j - burnin) %/% thin] <- - length(B_star) * 0.5 * log(2*pi) + sum(log(diag(theta.prec.chol))) -
        0.5 * t(B_star - b_star) %*% (theta.prior.prec + crossprod(x.tilde)) %*% (B_star - b_star)
    }

    # Sample sigma
    if(!fix_Sigma) {
      sigma2 <- rep(0,K)
      u_ort <- A %*% ((yt - B %*% xt - D %*% (w-mu.xi))/ w_sqrt)    # change from Gaussian
      sigma_post_a <- sigma0_T0 + rep(t_max,K)
      sigma_post_b <- sigma0_S0 + rowSums(u_ort^2)
      for (i in c(1:K)){
        sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
      }
      sigma <- sqrt(sigma2)
    }

    if (cal_Sigma & (j > burnin) & (j %% thin == 0)){
      u_ort <- A %*% ((yt - B %*% xt - D %*% (w-mu.xi))/ w_sqrt)    # change from Gaussian
      sigma_post_a <- prior$sigma_T0 + rep(t_max,K)
      sigma_post_b <- prior$sigma_S0 + rowSums(u_ort^2)
      lpost[(j - burnin) %/% thin] <- sum(dinvgamma(inits$sigma^2, shape = sigma_post_a * 0.5, rate = sigma_post_b * 0.5, log = T))
    }

    # Sample A0
    if(!fix_A) {
      u_std <- (yt - B %*% xt - D%*% (w-mu.xi))/ w_sqrt # change from Gaussian
      u_neg <- - u_std
      a_sample <- rep(0, K * (K - 1) /2)
      if (K > 1) {
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- a_prior[id_start:id_end]
          V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
          a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sigma[i],
                                                       xsub = matrix(u_neg[1:(i-1),] / sigma[i], nrow = i-1),
                                                       a_sub = a_sub,
                                                       V_a_sub = V_a_sub)
        }
      }
      A_post <- matrix(0, nrow = K, ncol = K)
      A_post[upper.tri(A)] <- a_sample
      A <- t(A_post)
      diag(A) <- 1
    }

    if (cal_A & (j > burnin) & (j %% thin == 0)){
      u_std <- (yt - B %*% xt - D%*% (w-mu.xi))/ w_sqrt # change from Gaussian
      u_neg <- - u_std
      A_star <- inits$A0[lower.tri(A)]
      if (K > 1) {
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- prior$a_prior[id_start:id_end]
          V_a_sub <- prior$V_a_prior[id_start:id_end, id_start:id_end]
          sub_A_star <- A_star[id_start:id_end]

          ysub = u_std[i,] / sigma[i]
          xsub = matrix(u_neg[1:(i-1),] / sigma[i], nrow = i-1)

          n = nrow(xsub)

          V_a_sub_inv = solve(V_a_sub)

          V_a_post_inv <- V_a_sub_inv + xsub %*% t(xsub)
          V_a_post <- solve(V_a_post_inv)
          a_post <- V_a_post %*% ( V_a_sub_inv %*% a_sub + xsub %*% ysub)

          lpost[(j - burnin) %/% thin] <- lpost[(j - burnin) %/% thin] + mvnfast::dmvn(sub_A_star, mu = a_post, sigma = V_a_post, log = T)
        }
      }
    }
    # Sample w
    q2 <- colSums( ( ( A %*% ( yt - B%*%xt + mu.xi * gamma) ) / sigma )^2 )
    p2 <- sum( ( c( A %*% gamma ) / sigma )^2 )
    w_sample <- mapply( GIGrvg::rgig, n = 1, lambda = -(nu+K)*0.5, chi = nu + q2,
                        psi = p2 )
    w <- reprow(w_sample, K)
    w_sqrt <- sqrt(w)

    # Sample nu
    if(!fix_Nu) {
      nu_temp = nu + exp(logsigma_nu)*rnorm(1)
      if (nu_temp > 4 && nu_temp < 100){
        num_mh = dgamma(nu_temp, shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w_sample, shape = nu_temp*0.5, rate = nu_temp*0.5, log = T))
        denum_mh = dgamma(nu, shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w_sample, shape = nu*0.5, rate = nu*0.5, log = T))
        alpha = num_mh - denum_mh;
        temp = log(runif(1));
        if (alpha > temp){
          nu = nu_temp
          acount_nu = acount_nu + 1
        }

      }
      mu.xi <- nu/( nu - 2 )
    }

    if (cal_Nu & (j > burnin) & (j %% thin == 0)){
      # Denominator
      if (fix_Nu) {
        nu_temp = inits$nu + exp(logsigma_nu)*rnorm(1)
        if (nu_temp > 4 && nu_temp < 100){
          num_mh = dgamma(nu_temp, shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w_sample, shape = nu_temp*0.5, rate = nu_temp*0.5, log = T))

          denum_mh = dgamma(inits$nu, shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w_sample, shape = inits$nu*0.5, rate = inits$nu*0.5, log = T))
          alpha = min(1, exp(num_mh - denum_mh))
          temp = runif(1)
          if (alpha > temp){
            acount_nu = acount_nu + 1
          }
          lpost[(j - burnin) %/% thin] <- log(alpha) # Take a negative for denominator

        }
      } else {
        # Nominator
        q_nu_nustar <- dnorm(inits$nu, mean = nu, sd = exp(logsigma_nu) ) # Move nu to inits$nu
        num_mh = dgamma(inits$nu, shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w_sample, shape = inits$nu*0.5, rate = inits$nu*0.5, log = T))
        denum_mh = dgamma(nu, shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w_sample, shape = nu*0.5, rate = nu*0.5, log = T))
        alpha = min(1, exp(num_mh - denum_mh))
        lpost[(j - burnin) %/% thin] <- log( alpha * q_nu_nustar)
      }
    }


    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, gamma, nu, as.numeric(w_sample))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", round(acount_nu/j,2)," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2), " ", round(gamma,2), " \n")
      acount_w <- rep(0,t_max)
    }

    # if(j %% batchlength == 0 ){
    #   if (acount_nu > batchlength * TARGACCEPT){
    #     logsigma_nu = logsigma_nu + adaptamount(j %/% batchlength);
    #   }
    #   if (acount_nu < batchlength * TARGACCEPT){
    #     logsigma_nu = logsigma_nu - adaptamount(j %/% batchlength);
    #   }
    #   acount_nu = 0
    # }
  }
  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("sigma",c(1:K), sep = ""),
                        paste("gamma",c(1:K), sep = ""),
                        paste("nu"),
                        paste("w",c(1:t_max), sep = ""))

  return(lpost)
}

###########################################################################
#' @export
Chib.MT.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL,
                          samples = 1000, burnin = 0, thin = 1,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                          cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)

  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "MT", SV = FALSE)
  }
  # prior B
  b_prior = prior$b_prior
  V_b_prior = prior$V_b_prior
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior A
  a_prior = prior$a_prior
  V_a_prior = prior$V_a_prior
  # prior nu
  nu_gam_a = prior$nu_gam_a
  nu_gam_b = prior$nu_gam_b

  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  #samples <- inits$samples
  # burnin <- inits$burnin
  # thin <- inits$thin
  A <- inits$A0
  B <- Vec_to_Mat(inits$B0, K, p)
  b_sample <- as.numeric(B)
  a_sample <- as.numeric(A[lower.tri(A)])
  sigma <- inits$sigma
  # Sigma <- solve(inits$A0) %*% diag(sigma, nrow = length(sigma))
  # Sigma2 <- Sigma %*% t(Sigma)
  # Sigma2_inv <- solve(Sigma2)
  V_b_prior_inv <- solve(V_b_prior)

  # Multi degrees of freedom
  nu <- inits$nu
  logsigma_nu <- inits$logsigma_nu
  acount_nu <- rep(0,K)
  acount_w <- rep(0, t_max)
  # Init w as Gaussian
  w <- inits$w
  w_sqrt <- sqrt(w)

  # new precompute
  theta.prior.precmean <- V_b_prior_inv %*% b_prior

  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K*t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  lpost <- rep(0, (samples - burnin)%/% thin)
  if (cal_Nu) lpost <- matrix(0,nrow = (samples - burnin)%/% thin, ncol = K)
  for (j in c(1:samples)){

    # Sample B
    if(!fix_B) {
      wt <- rep(1/sigma, t_max)
      y.tilde <- as.vector( A %*% ( w^(-0.5) * yt ) ) * wt
      # make and stack X matrices
      x.tilde <- kronecker( t(xt), A )
      tmp <- kronecker( t(w^(-1/2)), matrix(1,K,1) ) * wt # W^(-1/2)
      ii <- 0
      for (i in 1:m) {
        x.tilde[,(ii+1):(ii+K)] <- x.tilde[,(ii+1):(ii+K)] * tmp
        ii <- ii+K
      }
      theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
      b_sample <- backsolve( theta.prec.chol,
                             backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                        upper.tri = T, transpose = T )
                             + rnorm(K*m) )
      B <- matrix(b_sample,K,m)
    }
    if (cal_B & (j > burnin) & (j %% thin == 0)){
      B_star <- vec(inits$B)
      wt <- rep(1/sigma, t_max)
      y.tilde <- as.vector( A %*% ( w^(-0.5) * yt ) ) * wt
      # make and stack X matrices
      x.tilde <- kronecker( t(xt), A )
      tmp <- kronecker( t(w^(-1/2)), matrix(1,K,1) ) * wt # W^(-1/2)
      ii <- 0
      for (i in 1:m) {
        x.tilde[,(ii+1):(ii+K)] <- x.tilde[,(ii+1):(ii+K)] * tmp
        ii <- ii+K
      }
      theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )

      theta.prior.precmean <- V_b_prior_inv %*% prior$b_prior
      b_star <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T ))
      #mvnfast::dmvn(B_star, mu = b_star, sigma = (solve(V_b_prior_inv + crossprod(x.tilde) )), log = T)
      lpost[(j - burnin) %/% thin] <- - length(B_star) * 0.5 * log(2*pi) + sum(log(diag(theta.prec.chol))) -
                                      0.5 * t(B_star - b_star) %*% (V_b_prior_inv + crossprod(x.tilde)) %*% (B_star - b_star)
    }
    # Sample sigma
    if(!fix_Sigma) {
      sigma2 <- rep(0,K)
      u_ort <- A %*% ((yt - B %*%xt)/ w_sqrt)    # change from Gaussian
      sigma_post_a <- sigma0_T0 + rep(t_max,K)
      # sigma_post_b <- sigma0_S0 + apply(u_ort^2, 1, sum)
      sigma_post_b <- sigma0_S0 + rowSums(u_ort^2)
      for (i in c(1:K)){
        sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
      }
      sigma <- sqrt(sigma2)
    }

    if (cal_Sigma & (j > burnin) & (j %% thin == 0)){
      u_ort <- A %*% ((yt - B %*%xt)/ w_sqrt)    # change from Gaussian
      sigma_post_a <- prior$sigma_T0 + rep(t_max,K)
      sigma_post_b <- prior$sigma_S0 + rowSums(u_ort^2)
      lpost[(j - burnin) %/% thin] <- sum(dinvgamma(inits$sigma^2, shape = sigma_post_a * 0.5, rate = sigma_post_b * 0.5, log = T))
    }

    # Sample A0
    if(!fix_A) {
      u_std <- (yt - B %*%xt) / w_sqrt # change from Gaussian
      u_neg <- - u_std
      a_sample <- rep(0, K * (K - 1) /2)
      if (K > 1) {
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- a_prior[id_start:id_end]
          V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
          a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sigma[i],
                                                       xsub = matrix(u_neg[1:(i-1),] / sigma[i], nrow = i-1),
                                                       a_sub = a_sub,
                                                       V_a_sub = V_a_sub)
        }
      }
      A_post <- matrix(0, nrow = K, ncol = K)
      A_post[upper.tri(A)] <- a_sample
      A <- t(A_post)
      diag(A) <- 1
    }

    if (cal_A & (j > burnin) & (j %% thin == 0)){
      u_std <- (yt - B %*%xt) / w_sqrt # change from Gaussian
      u_neg <- - u_std
      A_star <- inits$A0[lower.tri(A)]
      if (K > 1) {
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- prior$a_prior[id_start:id_end]
          V_a_sub <- prior$V_a_prior[id_start:id_end, id_start:id_end]
          sub_A_star <- A_star[id_start:id_end]

          ysub = u_std[i,] / sigma[i]
          xsub = matrix(u_neg[1:(i-1),] / sigma[i], nrow = i-1)

          n = nrow(xsub)

          V_a_sub_inv = solve(V_a_sub)

          V_a_post_inv <- V_a_sub_inv + xsub %*% t(xsub)
          V_a_post <- solve(V_a_post_inv)
          a_post <- V_a_post %*% ( V_a_sub_inv %*% a_sub + xsub %*% ysub)

          lpost[(j - burnin) %/% thin] <- lpost[(j - burnin) %/% thin] + mvnfast::dmvn(sub_A_star, mu = a_post, sigma = V_a_post, log = T)
        }
      }
    }

    Sigma <- solve(A) %*% diag(sigma, nrow = length(sigma))
    Sigma2 <- Sigma %*% t(Sigma)
    Sigma2_inv <- solve(Sigma2)
    Sigma2_inv <- (Sigma2_inv + t(Sigma2_inv))*0.5

    # Sample w
    u <- (yt - B %*%xt)
    u_proposal <- A %*% (yt - B %*%xt)

    a_target <- (nu*0.5 + 1*0.5) * 0.75 # adjust by 0.75
    b_target <- (nu*0.5 + 0.5 * u_proposal^2 / sigma^2) * 0.75 # adjust by 0.75
    w_temp <- matrix(rinvgamma( n = K*t_max, shape = a_target, rate = b_target), nrow = K)
    w_temp_sqrt <- sqrt(w_temp)
    prior_num <- colSums(dinvgamma(w_temp, shape = nu*0.5, rate = nu*0.5, log = T))
    prior_denum <- colSums(dinvgamma(w, shape = nu*0.5, rate = nu*0.5, log = T))
    proposal_num <- colSums(dinvgamma(w_temp, shape = a_target, rate = b_target, log = T))
    proposal_denum <- colSums(dinvgamma(w, shape = a_target, rate = b_target, log = T))

    A.w.u.s.prop <- ( A %*% ( u/w_temp_sqrt ) )*( 1/ sigma )
    num_mh <- prior_num - proposal_num - 0.5 * colSums(log(w_temp)) - 0.5 * colSums(A.w.u.s.prop^2)
    A.w.u.s.curr <- ( A %*% ( u/w_sqrt ) )*( 1/ sigma )
    denum_mh <- prior_denum - proposal_denum - 0.5 * colSums(log(w)) - 0.5 * colSums(A.w.u.s.curr^2)
    alpha = num_mh - denum_mh
    temp = log(runif(t_max))
    acre = alpha > temp
    w[,acre] <- w_temp[,acre]
    w_sqrt[,acre] <- w_temp_sqrt[,acre]
    acount_w[acre] <- acount_w[acre] + 1

    # Sample nu
    if(!fix_Nu) {
      nu_temp = nu + exp(logsigma_nu)*rnorm(K)
      for (k in c(1:K)){
        if (nu_temp[k] > 4 && nu_temp[k] < 100){
          num_mh = dgamma(nu_temp[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = nu_temp[k]*0.5, rate = nu_temp[k]*0.5, log = T))
          denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))
          alpha = num_mh - denum_mh;
          temp = log(runif(1));
          if (alpha > temp){
            nu[k] = nu_temp[k]
            acount_nu[k] = acount_nu[k] + 1
          }

        }
      }
    }
    if (cal_Nu & (j > burnin) & (j %% thin == 0)){
      # Denominator
      if (fix_Nu) {

        for (k in c(1:K)){
          nu_temp = inits$nu[k] + exp(logsigma_nu[k])*rnorm(1)
          if (nu_temp > 4 && nu_temp < 100){
            num_mh = dgamma(nu_temp, shape = nu_gam_a, rate = nu_gam_b, log = T) +
              sum(dinvgamma(w[k,], shape = nu_temp*0.5, rate = nu_temp*0.5, log = T))
            denum_mh = dgamma(inits$nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
              sum(dinvgamma(w[k,], shape = inits$nu[k]*0.5, rate = inits$nu[k]*0.5, log = T))
            alpha = min(1, exp(num_mh - denum_mh))
            temp = runif(1)
            if (alpha > temp){
              acount_nu[k] = acount_nu[k] + 1
            }
            lpost[(j - burnin) %/% thin,k] <- log(alpha) # Take a negative for denominator

          }
        }


      } else {
        # Nominator
        for (k in c(1:K)){
          q_nu_nustar <- dnorm(inits$nu[k], mean = nu[k], sd = exp(logsigma_nu[k]) ) # Move nu to inits$nu
          num_mh = dgamma(inits$nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = inits$nu[k]*0.5, rate = inits$nu[k]*0.5, log = T))
          denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))
          alpha = min(1, exp(num_mh - denum_mh))
          lpost[(j - burnin) %/% thin,k] <- log( alpha * q_nu_nustar)
        }
      }
    }
    # if(j %% batchlength == 0 ){
    #   for (jj in c(1:K)) {
    #     if (acount_nu[jj] > batchlength * TARGACCEPT){
    #       logsigma_nu[jj] = logsigma_nu[jj] + adaptamount(j %/% batchlength);
    #     }
    #     if (acount_nu[jj] < batchlength * TARGACCEPT){
    #       logsigma_nu[jj] = logsigma_nu[jj] - adaptamount(j %/% batchlength);
    #     }
    #     acount_nu[jj] = 0
    #   }
    # }
    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, nu, as.numeric(w))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", round(acount_nu/j,2)," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2), " \n")
      acount_w <- rep(0,t_max)
    }
  }
  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("sigma",c(1:K), sep = ""),
                        paste("nu",c(1:K), sep = ""),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K))
  )

  return(lpost)
}

###########################################################################

#' @export
Chib.MST.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL,
                           samples = 1000, burnin = 0, thin = 1,
                           fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                           cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)

  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "MST", SV = FALSE)
  }
  # prior B
  b_prior = prior$b_prior
  V_b_prior = prior$V_b_prior
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior A
  a_prior = prior$a_prior
  V_a_prior = prior$V_a_prior
  # prior nu
  nu_gam_a = prior$nu_gam_a
  nu_gam_b = prior$nu_gam_b
  # prior gamma
  gamma_prior = prior$gamma_prior
  V_gamma_prior = prior$V_gamma_prior
  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  #samples <- inits$samples
  # burnin <- inits$burnin
  # thin <- inits$thin
  A <- inits$A0
  B <- Vec_to_Mat(inits$B0, K, p)
  b_sample <- as.numeric(B)
  a_sample <- as.numeric(A[lower.tri(A)])
  sigma <- inits$sigma
  # Sigma <- solve(inits$A0) %*% diag(sigma, nrow = length(sigma))
  # Sigma2 <- Sigma %*% t(Sigma)
  # Sigma2_inv <- solve(Sigma2)
  V_b_prior_inv <- solve(V_b_prior)

  # Multi degrees of freedom
  nu <- inits$nu
  mu.xi <- nu/( nu - 2 )
  logsigma_nu <- inits$logsigma_nu
  acount_nu <- rep(0,K)
  acount_w <- rep(0, t_max)
  gamma <- inits$gamma
  D <- diag(gamma)

  # Init w as Gaussian
  w <- inits$w
  w_sqrt <- sqrt(w)
  #w_sqrt_inv <- 1/w_sqrt

  # new precompute
  i <- K*m
  theta.prior.prec = matrix(0,i+K,i+K)
  theta.prior.prec[1:i,1:i] <- V_b_prior_inv
  theta.prior.prec[i+(1:K),i+(1:K)] <- solve( V_gamma_prior )
  theta.prior.precmean <- theta.prior.prec %*% c( b_prior, gamma_prior )

  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K*t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  lpost <- rep(0, (samples - burnin)%/% thin)
  if (cal_Nu) lpost <- matrix(0,nrow = (samples - burnin)%/% thin, ncol = K)
  for (j in c(1:samples)){
    if(!fix_B) {
      # Sample B and gamma
      wt <- rep(1/sigma,t_max)
      y.tilde <- as.vector( A %*% ( w^(-0.5) * yt ) ) * wt
      # make and stack X matrices
      x1 <- kronecker( t(xt), A )
      tmp <- kronecker( t(w^(-1/2)), matrix(1,K,1) ) * wt # W^(-1/2)
      ii <- 0
      for (i in 1:m) {
        x1[,(ii+1):(ii+K)] <- x1[,(ii+1):(ii+K)] * tmp
        ii <- ii+K
      }
      W.mat <- matrix(0,K*t_max,K)
      idx <- (0:(t_max-1))*K
      for ( i in 1:K ) {
        W.mat[idx + (i-1)*t_max*K + i] <- (w[i,] - mu.xi[i])/sqrt(w[i,]) # W^(-1/2) * (W - bar(W))
      }
      x2 <- matrix(0,K*t_max,K)
      for ( i in 1:K ) {
        x2[,i] <- as.vector( A%*%matrix(W.mat[,i],K,t_max) ) * wt
      }
      x.tilde <- cbind( x1, x2 )
      theta.prec.chol <- chol( theta.prior.prec + crossprod(x.tilde) )
      theta <- backsolve( theta.prec.chol,
                          backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                     upper.tri = T, transpose = T )
                          + rnorm(K*(m+1)) )
      B <- matrix(theta[1:(m*K)],K,m)
      b_sample <- as.vector(B)
      gamma <- theta[(m*K+1):((m+1)*K)]
      D <- diag(gamma)
    }

    if (cal_B & (j > burnin) & (j %% thin == 0)){
      B_star <- c(vec(inits$B), inits$gamma)

      # Sample B and gamma
      wt <- rep(1/sigma,t_max)
      y.tilde <- as.vector( A %*% ( w^(-0.5) * yt ) ) * wt
      # make and stack X matrices
      x1 <- kronecker( t(xt), A )
      tmp <- kronecker( t(w^(-1/2)), matrix(1,K,1) ) * wt # W^(-1/2)
      ii <- 0
      for (i in 1:m) {
        x1[,(ii+1):(ii+K)] <- x1[,(ii+1):(ii+K)] * tmp
        ii <- ii+K
      }
      W.mat <- matrix(0,K*t_max,K)
      idx <- (0:(t_max-1))*K
      for ( i in 1:K ) {
        W.mat[idx + (i-1)*t_max*K + i] <- (w[i,] - mu.xi[i])/sqrt(w[i,]) # W^(-1/2) * (W - bar(W))
      }
      x2 <- matrix(0,K*t_max,K)
      for ( i in 1:K ) {
        x2[,i] <- as.vector( A%*%matrix(W.mat[,i],K,t_max) ) * wt
      }
      x.tilde <- cbind( x1, x2 )

      theta.prec.chol <- chol( theta.prior.prec + crossprod(x.tilde) )
      b_star <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T ))
      #mvnfast::dmvn(B_star, mu = b_star, sigma = (solve(V_b_prior_inv + crossprod(x.tilde) )), log = T)
      lpost[(j - burnin) %/% thin] <- - length(B_star) * 0.5 * log(2*pi) + sum(log(diag(theta.prec.chol))) -
        0.5 * t(B_star - b_star) %*% (theta.prior.prec + crossprod(x.tilde)) %*% (B_star - b_star)
    }

    # Sample sigma
    if(!fix_Sigma) {
      sigma2 <- rep(0,K)
      u_ort <- A %*% ((yt - B %*% xt - D %*% ( w - mu.xi ))/ w_sqrt)    # change from Gaussian
      sigma_post_a <- sigma0_T0 + rep(t_max,K)
      sigma_post_b <- sigma0_S0 + rowSums(u_ort^2)
      for (i in c(1:K)){
        sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
      }
      sigma <- sqrt(sigma2)
    }

    if (cal_Sigma & (j > burnin) & (j %% thin == 0)){
      u_ort <- A %*% ((yt - B %*% xt - D %*% ( w - mu.xi ))/ w_sqrt)    # change from Gaussian
      sigma_post_a <- prior$sigma_T0 + rep(t_max,K)
      sigma_post_b <- prior$sigma_S0 + rowSums(u_ort^2)
      lpost[(j - burnin) %/% thin] <- sum(dinvgamma(inits$sigma^2, shape = sigma_post_a * 0.5, rate = sigma_post_b * 0.5, log = T))
    }

    # Sample A0
    if(!fix_A) {
      u_std <- (yt - B %*% xt - D%*% ( w - mu.xi )) / w_sqrt # change from Gaussian
      u_neg <- - u_std
      a_sample <- rep(0, K * (K - 1) /2)
      if (K > 1) {
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- a_prior[id_start:id_end]
          V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
          a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sigma[i],
                                                       xsub = matrix(u_neg[1:(i-1),] / sigma[i], nrow = i-1),
                                                       a_sub = a_sub,
                                                       V_a_sub = V_a_sub)
        }
      }

      A_post <- matrix(0, nrow = K, ncol = K)
      A_post[upper.tri(A)] <- a_sample
      A <- t(A_post)
      diag(A) <- 1
    }

    if (cal_A & (j > burnin) & (j %% thin == 0)){
      u_std <- (yt - B %*% xt - D%*% ( w - mu.xi )) / w_sqrt # change from Gaussian
      u_neg <- - u_std
      A_star <- inits$A0[lower.tri(A)]
      if (K > 1) {
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- prior$a_prior[id_start:id_end]
          V_a_sub <- prior$V_a_prior[id_start:id_end, id_start:id_end]
          sub_A_star <- A_star[id_start:id_end]

          ysub = u_std[i,] / sigma[i]
          xsub = matrix(u_neg[1:(i-1),] / sigma[i], nrow = i-1)

          n = nrow(xsub)

          V_a_sub_inv = solve(V_a_sub)

          V_a_post_inv <- V_a_sub_inv + xsub %*% t(xsub)
          V_a_post <- solve(V_a_post_inv)
          a_post <- V_a_post %*% ( V_a_sub_inv %*% a_sub + xsub %*% ysub)

          lpost[(j - burnin) %/% thin] <- lpost[(j - burnin) %/% thin] + mvnfast::dmvn(sub_A_star, mu = a_post, sigma = V_a_post, log = T)
        }
      }
    }

    Sigma2_inv <- t(A)%*%diag(1/sigma^2)%*%A

    # Sample w
    u <- (yt - B %*%xt + mu.xi * gamma)
    u_proposal <- A %*% (yt - B %*%xt) + mu.xi * gamma

    a_target <- (nu*0.5 + 1*0.5) * 0.75 # adjust by 0.75
    b_target <- (nu*0.5 + 0.5 * u_proposal^2 / sigma^2) * 0.75 # adjust by 0.75
    w_temp <- matrix(rinvgamma( n = K*t_max, shape = a_target, rate = b_target), nrow = K)
    log_w_temp <- log(w_temp)
    w_temp_sqrt <- sqrt(w_temp)
    A.w.u.s.prop <- ( A %*% ( u/w_temp_sqrt - w_temp_sqrt*gamma ) )/sigma
    num_mh <- colSums(dinvgamma(w_temp, shape = nu*0.5, rate = nu*0.5, log = T)) - #prior
      0.5*( colSums(log(w_temp)) + colSums(A.w.u.s.prop^2) ) - # posterior
      colSums( dinvgamma(w_temp, shape = a_target, rate = b_target, log = T)) # proposal
    A.w.u.s.curr <- ( A %*% ( u/w_sqrt - w_sqrt*gamma ) )/sigma
    denum_mh <- colSums(dinvgamma(w, shape = nu*0.5, rate = nu*0.5, log = T)) - #prior
      0.5*( colSums(log(w)) + colSums(A.w.u.s.curr^2) ) - # posterior
      colSums( dinvgamma(w, shape = a_target, rate = b_target, log = T)) # proposal
    acc <- (num_mh-denum_mh) > log(runif(t_max))
    w[,acc] <- w_temp[,acc]
    w_sqrt[,acc] <- w_temp_sqrt[,acc]
    acount_w[acc] <- acount_w[acc] + 1



    # Sample nu
    if(!fix_Nu) {
      nu_temp = nu + exp(logsigma_nu)*rnorm(K)
      for (k in c(1:K)){
        if (nu_temp[k] > 4 && nu_temp[k] < 100){
          num_mh = dgamma(nu_temp[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = nu_temp[k]*0.5, rate = nu_temp[k]*0.5, log = T))
          denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))
          alpha = num_mh - denum_mh;
          temp = log(runif(1));
          if (alpha > temp){
            nu[k] = nu_temp[k]
            acount_nu[k] = acount_nu[k] + 1
          }

        }

      }
      mu.xi <- nu/( nu - 2 )
    }
    if (cal_Nu & (j > burnin) & (j %% thin == 0)){
      # Denominator
      if (fix_Nu) {

        for (k in c(1:K)){
          nu_temp = inits$nu[k] + exp(logsigma_nu[k])*rnorm(1)
          if (nu_temp > 4 && nu_temp < 100){
            num_mh = dgamma(nu_temp, shape = nu_gam_a, rate = nu_gam_b, log = T) +
              sum(dinvgamma(w[k,], shape = nu_temp*0.5, rate = nu_temp*0.5, log = T))
            denum_mh = dgamma(inits$nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
              sum(dinvgamma(w[k,], shape = inits$nu[k]*0.5, rate = inits$nu[k]*0.5, log = T))
            alpha = min(1, exp(num_mh - denum_mh))
            temp = runif(1)
            if (alpha > temp){
              acount_nu[k] = acount_nu[k] + 1
            }
            lpost[(j - burnin) %/% thin,k] <- log(alpha) # Take a negative for denominator

          }
        }


      } else {
        # Nominator
        for (k in c(1:K)){
          q_nu_nustar <- dnorm(inits$nu[k], mean = nu[k], sd = exp(logsigma_nu[k]) ) # Move nu to inits$nu
          num_mh = dgamma(inits$nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = inits$nu[k]*0.5, rate = inits$nu[k]*0.5, log = T))
          denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))
          alpha = min(1, exp(num_mh - denum_mh))
          lpost[(j - burnin) %/% thin,k] <- log( alpha * q_nu_nustar)
        }
      }
    }
    # if(j %% batchlength == 0 ){
    #   for (jj in c(1:K)) {
    #     if (acount_nu[jj] > batchlength * TARGACCEPT){
    #       logsigma_nu[jj] = logsigma_nu[jj] + adaptamount(j %/% batchlength);
    #     }
    #     if (acount_nu[jj] < batchlength * TARGACCEPT){
    #       logsigma_nu[jj] = logsigma_nu[jj] - adaptamount(j %/% batchlength);
    #     }
    #     acount_nu[jj] = 0
    #   }
    # }

    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, gamma, nu, as.numeric(w))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", round(acount_nu/j,2)," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2) , " ", round(gamma,2),  " \n")
      acount_w <- rep(0,t_max)
    }
  }
  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("sigma",c(1:K), sep = ""),
                        paste("gamma",c(1:K), sep = ""),
                        paste("nu",c(1:K), sep = ""),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K))
  )

  return(lpost)
}

###########################################################################
#' @export
Chib.OT.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL,
                          samples = 1000, burnin = 0, thin = 1,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                          cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)

  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "OT", SV = FALSE)
  }
  # prior B
  b_prior = prior$b_prior
  V_b_prior = prior$V_b_prior
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior A
  a_prior = prior$a_prior
  V_a_prior = prior$V_a_prior
  # prior nu
  nu_gam_a = prior$nu_gam_a
  nu_gam_b = prior$nu_gam_b

  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  #samples <- inits$samples
  # burnin <- inits$burnin
  # thin <- inits$thin
  A <- inits$A0
  B <- Vec_to_Mat(inits$B0, K, p)
  b_sample <- as.numeric(B)
  a_sample <- as.numeric(A[lower.tri(A)])
  sigma <- inits$sigma
  sigma2 <- sigma^2

  V_b_prior_inv <- solve(V_b_prior)

  # Multi degrees of freedom
  nu <- inits$nu
  logsigma_nu <- inits$logsigma_nu
  acount_nu <- rep(0,K)
  acount_w <- rep(0, t_max)
  # Init w as Gaussian
  w <- inits$w
  w_sqrt <- sqrt(w)
  w_sqrt_inv <- 1/w_sqrt

  # new precompute
  theta.prior.precmean <- V_b_prior_inv %*% b_prior

  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K*t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  lpost <- rep(0, (samples - burnin)%/% thin)
  if (cal_Nu) lpost <- matrix(0,nrow = (samples - burnin)%/% thin, ncol = K)
  for (j in c(1:samples)){

    # Sample B
    if(!fix_B) {
      wt <- rep( 1/sigma, t_max ) / as.numeric(w_sqrt)
      y.tilde <- as.vector( A %*% yt ) * wt
      x.tilde <- kronecker( t(xt), A ) * wt
      theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
      b_sample <- backsolve( theta.prec.chol,
                             backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                        upper.tri = T, transpose = T )
                             + rnorm(K*m) )
      B <- matrix(b_sample,K,m)
    }
    if (cal_B & (j > burnin) & (j %% thin == 0)){
      B_star <- vec(inits$B)
      y.tilde <- as.vector( A %*% yt ) * wt
      x.tilde <- kronecker( t(xt), A ) * wt
      theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )

      theta.prior.precmean <- V_b_prior_inv %*% prior$b_prior
      b_star <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T ))
      #mvnfast::dmvn(B_star, mu = b_star, sigma = (solve(V_b_prior_inv + crossprod(x.tilde) )), log = T)
      lpost[(j - burnin) %/% thin] <- - length(B_star) * 0.5 * log(2*pi) + sum(log(diag(theta.prec.chol))) -
        0.5 * t(B_star - b_star) %*% (V_b_prior_inv + crossprod(x.tilde)) %*% (B_star - b_star)
    }

    # Sample sigma
    if(!fix_Sigma) {
      u_ort <- (A %*% (yt - B %*%xt))/ w_sqrt    # change from Gaussian
      sigma_post_a <- sigma0_T0 + rep(t_max,K)
      sigma_post_b <- sigma0_S0 + apply(u_ort^2, 1, sum)
      for (i in c(1:K)){
        sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
      }
      sigma <- sqrt(sigma2)
    }

    if (cal_Sigma & (j > burnin) & (j %% thin == 0)){
      u_ort <- (A %*% (yt - B %*%xt))/ w_sqrt    # change from Gaussian
      sigma_post_a <- prior$sigma_T0 + rep(t_max,K)
      sigma_post_b <- prior$sigma_S0 + rowSums(u_ort^2)
      lpost[(j - burnin) %/% thin] <- sum(dinvgamma(inits$sigma^2, shape = sigma_post_a * 0.5, rate = sigma_post_b * 0.5, log = T))
    }

    # Sample A0
    if(!fix_A) {
      u_std <- (yt - B %*%xt) # change from Gaussian
      u_neg <- - u_std
      a_sample <- rep(0, K * (K - 1) /2)
      if (K > 1) {
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- a_prior[id_start:id_end]
          V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
          a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sigma[i]/ w_sqrt[i,],
                                                       xsub = matrix(u_neg[1:(i-1),] / sigma[i] / reprow(w_sqrt[i,], i-1), nrow = i-1),
                                                       a_sub = a_sub,
                                                       V_a_sub = V_a_sub)
        }
      }
      A_post <- matrix(0, nrow = K, ncol = K)
      A_post[upper.tri(A)] <- a_sample
      A <- t(A_post)
      diag(A) <- 1
    }

    if (cal_A & (j > burnin) & (j %% thin == 0)){
      u_std <- (yt - B %*%xt) # change from Gaussian
      u_neg <- - u_std
      A_star <- inits$A0[lower.tri(A)]
      if (K > 1) {
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- prior$a_prior[id_start:id_end]
          V_a_sub <- prior$V_a_prior[id_start:id_end, id_start:id_end]
          sub_A_star <- A_star[id_start:id_end]

          ysub = u_std[i,] / sigma[i]/ w_sqrt[i,]
          xsub = matrix(u_neg[1:(i-1),] / sigma[i] / reprow(w_sqrt[i,], i-1), nrow = i-1)

          n = nrow(xsub)

          V_a_sub_inv = solve(V_a_sub)

          V_a_post_inv <- V_a_sub_inv + xsub %*% t(xsub)
          V_a_post <- solve(V_a_post_inv)
          a_post <- V_a_post %*% ( V_a_sub_inv %*% a_sub + xsub %*% ysub)

          lpost[(j - burnin) %/% thin] <- lpost[(j - burnin) %/% thin] + mvnfast::dmvn(sub_A_star, mu = a_post, sigma = V_a_post, log = T)
        }
      }
    }
    # Sample w
    u <-  A %*% (yt - B %*%xt)
    w <- matrix( rinvgamma( n = K*t_max, shape = (nu+1)*0.5, rate = 0.5*( nu + (u^2)/sigma2 ) ), K, t_max )
    w_sqrt <- sqrt(w)


    # Sample nu
    if(!fix_Nu) {
      nu_temp = nu + exp(logsigma_nu)*rnorm(K)
      for (k in c(1:K)){
        if (nu_temp[k] > 4 && nu_temp[k] < 100){
          num_mh = dgamma(nu_temp[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = nu_temp[k]*0.5, rate = nu_temp[k]*0.5, log = T))
          denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))
          alpha = num_mh - denum_mh;
          temp = log(runif(1));
          if (alpha > temp){
            nu[k] = nu_temp[k]
            acount_nu[k] = acount_nu[k] + 1
          }

        }
      }
    }
    if (cal_Nu & (j > burnin) & (j %% thin == 0)){
      # Denominator
      if (fix_Nu) {

        for (k in c(1:K)){
          nu_temp = inits$nu[k] + exp(logsigma_nu[k])*rnorm(1)
          if (nu_temp > 4 && nu_temp < 100){
            num_mh = dgamma(nu_temp, shape = nu_gam_a, rate = nu_gam_b, log = T) +
              sum(dinvgamma(w[k,], shape = nu_temp*0.5, rate = nu_temp*0.5, log = T))
            denum_mh = dgamma(inits$nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
              sum(dinvgamma(w[k,], shape = inits$nu[k]*0.5, rate = inits$nu[k]*0.5, log = T))
            alpha = min(1, exp(num_mh - denum_mh))
            temp = runif(1)
            if (alpha > temp){
              acount_nu[k] = acount_nu[k] + 1
            }
            lpost[(j - burnin) %/% thin,k] <- log(alpha) # Take a negative for denominator

          }
        }


      } else {
        # Nominator
        for (k in c(1:K)){
          q_nu_nustar <- dnorm(inits$nu[k], mean = nu[k], sd = exp(logsigma_nu[k]) ) # Move nu to inits$nu
          num_mh = dgamma(inits$nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = inits$nu[k]*0.5, rate = inits$nu[k]*0.5, log = T))
          denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))
          alpha = min(1, exp(num_mh - denum_mh))
          lpost[(j - burnin) %/% thin,k] <- log( alpha * q_nu_nustar)
        }
      }
    }
    # if(j %% batchlength == 0 ){
    #   for (jj in c(1:K)) {
    #     if (acount_nu[jj] > batchlength * TARGACCEPT){
    #       logsigma_nu[jj] = logsigma_nu[jj] + adaptamount(j %/% batchlength);
    #     }
    #     if (acount_nu[jj] < batchlength * TARGACCEPT){
    #       logsigma_nu[jj] = logsigma_nu[jj] - adaptamount(j %/% batchlength);
    #     }
    #     acount_nu[jj] = 0
    #   }
    # }
    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, nu, as.numeric(w))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", round(acount_nu/j,2)," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2), " \n")
      acount_w <- rep(0,t_max)
    }
  }
  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("sigma",c(1:K), sep = ""),
                        paste("nu",c(1:K), sep = ""),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K))
  )

  return(lpost)
}

###########################################################################
#' @export
Chib.OST.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL,
                           samples = 1000, burnin = 0, thin = 1,
                           fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                           cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)

  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "OST", SV = FALSE)
  }
  # prior B
  b_prior = prior$b_prior
  V_b_prior = prior$V_b_prior
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior A
  a_prior = prior$a_prior
  V_a_prior = prior$V_a_prior
  # prior nu
  nu_gam_a = prior$nu_gam_a
  nu_gam_b = prior$nu_gam_b
  # prior gamma
  gamma_prior = prior$gamma_prior
  V_gamma_prior = prior$V_gamma_prior
  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  #samples <- inits$samples
  # burnin <- inits$burnin
  # thin <- inits$thin
  A <- inits$A0
  B <- Vec_to_Mat(inits$B0, K, p)
  b_sample <- as.numeric(B)
  a_sample <- as.numeric(A[lower.tri(A)])
  sigma <- inits$sigma
  sigma2 <- sigma^2

  V_b_prior_inv <- solve(V_b_prior)

  # Multi degrees of freedom
  nu <- inits$nu
  mu.xi <- nu/( nu - 2 )
  logsigma_nu <- inits$logsigma_nu
  acount_nu <- rep(0,K)
  acount_w <- rep(0, t_max)
  gamma <- inits$gamma
  D <- diag(gamma)

  # new precompute
  i <- K*m
  theta.prior.prec = matrix(0,i+K,i+K)
  theta.prior.prec[1:i,1:i] <- V_b_prior_inv
  theta.prior.prec[i+(1:K),i+(1:K)] <- solve( V_gamma_prior )
  theta.prior.precmean <- theta.prior.prec %*% c( b_prior, gamma_prior )

  # Init w as Gaussian
  w <- inits$w
  w_sqrt <- sqrt(w)
  w_sqrt_inv <- 1/w_sqrt
  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K*t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  lpost <- rep(0, (samples - burnin)%/% thin)
  if (cal_Nu) lpost <- matrix(0,nrow = (samples - burnin)%/% thin, ncol = K)
  for (j in c(1:samples)){

    # Sample B and gamma
    if(!fix_B) {
      wt <- rep( 1/sigma, t_max ) / as.numeric(w_sqrt)
      y.tilde <- as.vector( A %*% yt ) * wt
      # make and stack diagonal W matrices
      W.mat <- matrix(0,K*t_max,K)
      idx <- (0:(t_max-1))*K
      for ( i in 1:K ) {
        W.mat[idx + (i-1)*t_max*K + i] <- w[i,] -  mu.xi[i]
      }
      x.tilde <- cbind( kronecker( t(xt), A ), W.mat ) * wt
      theta.prec.chol <- chol( theta.prior.prec + crossprod(x.tilde) )
      theta <- backsolve( theta.prec.chol,
                          backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                     upper.tri = T, transpose = T )
                          + rnorm(K*(m+1)) )
      B <- matrix(theta[1:(m*K)],K,m)
      b_sample <- as.vector(B)
      gamma <- theta[(m*K+1):((m+1)*K)]
      D <- diag(gamma)
    }
    if (cal_B & (j > burnin) & (j %% thin == 0)){
      B_star <- c(vec(inits$B), inits$gamma)

      # Sample B and gamma
      wt <- rep( 1/sigma, t_max ) / as.numeric(w_sqrt)
      y.tilde <- as.vector( A %*% yt ) * wt
      # make and stack diagonal W matrices
      W.mat <- matrix(0,K*t_max,K)
      idx <- (0:(t_max-1))*K
      for ( i in 1:K ) {
        W.mat[idx + (i-1)*t_max*K + i] <- w[i,] -  mu.xi[i]
      }
      x.tilde <- cbind( kronecker( t(xt), A ), W.mat ) * wt

      theta.prec.chol <- chol( theta.prior.prec + crossprod(x.tilde) )
      b_star <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T ))
      #mvnfast::dmvn(B_star, mu = b_star, sigma = (solve(V_b_prior_inv + crossprod(x.tilde) )), log = T)
      lpost[(j - burnin) %/% thin] <- - length(B_star) * 0.5 * log(2*pi) + sum(log(diag(theta.prec.chol))) -
        0.5 * t(B_star - b_star) %*% (theta.prior.prec + crossprod(x.tilde)) %*% (B_star - b_star)
    }

    # Sample sigma
    if(!fix_Sigma) {
      u_ort <- (A %*% (yt - B %*%xt) - D %*% ( w - mu.xi ))/ w_sqrt    # change from Gaussian
      sigma_post_a <- sigma0_T0 + rep(t_max,K)
      sigma_post_b <- sigma0_S0 + apply(u_ort^2, 1, sum)
      for (i in c(1:K)){
        sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
      }
      sigma <- sqrt(sigma2)
    }

    if (cal_Sigma & (j > burnin) & (j %% thin == 0)){
      u_ort <- (A %*% (yt - B %*%xt) - D %*% ( w - mu.xi ))/ w_sqrt    # change from Gaussian
      sigma_post_a <- prior$sigma_T0 + rep(t_max,K)
      sigma_post_b <- prior$sigma_S0 + rowSums(u_ort^2)
      lpost[(j - burnin) %/% thin] <- sum(dinvgamma(inits$sigma^2, shape = sigma_post_a * 0.5, rate = sigma_post_b * 0.5, log = T))
    }
    # Sample A0
    if(!fix_A) {
      u_std <- (yt - B %*%xt) # change from Gaussian
      u_neg <- - u_std
      a_sample <- rep(0, K * (K - 1) /2)
      if (K > 1) {
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- a_prior[id_start:id_end]
          V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
          a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = (u_std[i,] - (w[i,]  - mu.xi[i] ) * gamma[i]) / sigma[i]/ w_sqrt[i,],
                                                       xsub = matrix(u_neg[1:(i-1),] / sigma[i] / reprow(w_sqrt[i,], i-1), nrow = i-1),
                                                       a_sub = a_sub,
                                                       V_a_sub = V_a_sub)
        }
      }
      A_post <- matrix(0, nrow = K, ncol = K)
      A_post[upper.tri(A)] <- a_sample
      A <- t(A_post)
      diag(A) <- 1
    }

    if (cal_A & (j > burnin) & (j %% thin == 0)){
      u_std <- (yt - B %*%xt) # change from Gaussian
      u_neg <- - u_std
      A_star <- inits$A0[lower.tri(A)]
      if (K > 1) {
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- prior$a_prior[id_start:id_end]
          V_a_sub <- prior$V_a_prior[id_start:id_end, id_start:id_end]
          sub_A_star <- A_star[id_start:id_end]

          ysub = (u_std[i,] - (w[i,]  - mu.xi[i] ) * gamma[i]) / sigma[i]/ w_sqrt[i,]
          xsub = matrix(u_neg[1:(i-1),] / sigma[i] / reprow(w_sqrt[i,], i-1), nrow = i-1)

          n = nrow(xsub)

          V_a_sub_inv = solve(V_a_sub)

          V_a_post_inv <- V_a_sub_inv + xsub %*% t(xsub)
          V_a_post <- solve(V_a_post_inv)
          a_post <- V_a_post %*% ( V_a_sub_inv %*% a_sub + xsub %*% ysub)

          lpost[(j - burnin) %/% thin] <- lpost[(j - burnin) %/% thin] + mvnfast::dmvn(sub_A_star, mu = a_post, sigma = V_a_post, log = T)
        }
      }
    }
    # Sample w
    u <-  A %*% (yt - B %*%xt) + gamma * mu.xi
    w <- matrix( mapply( GIGrvg::rgig, n = 1, lambda = -(nu+1)*0.5, chi = nu + (u^2)/ sigma2,
                         psi = gamma^2/sigma2 ), K, t_max )
    w_sqrt <- sqrt(w)


    # Sample nu
    if(!fix_Nu) {
      nu_temp = nu + exp(logsigma_nu)*rnorm(K)
      for (k in c(1:K)){
        if (nu_temp[k] > 4 && nu_temp[k] < 100){
          num_mh = dgamma(nu_temp[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = nu_temp[k]*0.5, rate = nu_temp[k]*0.5, log = T))
          denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))
          alpha = num_mh - denum_mh;
          temp = log(runif(1));
          if (alpha > temp){
            nu[k] = nu_temp[k]
            acount_nu[k] = acount_nu[k] + 1
          }

        }
      }
      mu.xi <- nu/( nu - 2 )
    }
    if (cal_Nu & (j > burnin) & (j %% thin == 0)){
      # Denominator
      if (fix_Nu) {

        for (k in c(1:K)){
          nu_temp = inits$nu[k] + exp(logsigma_nu[k])*rnorm(1)
          if (nu_temp > 4 && nu_temp < 100){
            num_mh = dgamma(nu_temp, shape = nu_gam_a, rate = nu_gam_b, log = T) +
              sum(dinvgamma(w[k,], shape = nu_temp*0.5, rate = nu_temp*0.5, log = T))
            denum_mh = dgamma(inits$nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
              sum(dinvgamma(w[k,], shape = inits$nu[k]*0.5, rate = inits$nu[k]*0.5, log = T))
            alpha = min(1, exp(num_mh - denum_mh))
            temp = runif(1)
            if (alpha > temp){
              acount_nu[k] = acount_nu[k] + 1
            }
            lpost[(j - burnin) %/% thin,k] <- log(alpha) # Take a negative for denominator

          }
        }


      } else {
        # Nominator
        for (k in c(1:K)){
          q_nu_nustar <- dnorm(inits$nu[k], mean = nu[k], sd = exp(logsigma_nu[k]) ) # Move nu to inits$nu
          num_mh = dgamma(inits$nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = inits$nu[k]*0.5, rate = inits$nu[k]*0.5, log = T))
          denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))
          alpha = min(1, exp(num_mh - denum_mh))
          lpost[(j - burnin) %/% thin,k] <- log( alpha * q_nu_nustar)
        }
      }
    }
    # if(j %% batchlength == 0 ){
    #   for (jj in c(1:K)) {
    #     if (acount_nu[jj] > batchlength * TARGACCEPT){
    #       logsigma_nu[jj] = logsigma_nu[jj] + adaptamount(j %/% batchlength);
    #     }
    #     if (acount_nu[jj] < batchlength * TARGACCEPT){
    #       logsigma_nu[jj] = logsigma_nu[jj] - adaptamount(j %/% batchlength);
    #     }
    #     acount_nu[jj] = 0
    #   }
    # }
    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, gamma, nu, as.numeric(w))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", round(acount_nu/j,2)," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2), " ", round(gamma,2), " \n")
      acount_w <- rep(0,t_max)
    }
  }
  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("sigma",c(1:K), sep = ""),
                        paste("gamma",c(1:K), sep = ""),
                        paste("nu",c(1:K), sep = ""),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K))
  )

  return(lpost)
}
