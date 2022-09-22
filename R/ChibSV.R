#' @export
Chib.LML.SV <- function(Chain,
                        samples = 3000, burnin = 0, thin = 1, numCores = 4){
  # samples = 3000; burnin = 0; thin = 1;
  # numCores <- 16;

  P1 <- P2 <- P3 <- P4 <- P5 <- 0

  prior <- Chain$prior
  inits <- list()
  dist <- prior$dist

  B_mat <- get_post(Chain, element = "B")
  A_mat <- get_post(Chain, element = "a")
  Gamma_mat <- get_post(Chain, element = "gamma")
  Sigma_mat <- get_post(Chain, element = "sigma")
  Nu_mat <- get_post(Chain, element = "nu")
  W_mat <- get_post(Chain, element = "w")
  H_mat <- get_post(Chain, element = "h")
  H0_mat <- get_post(Chain, element = "lh0")

  B_med <- apply(B_mat, MARGIN = 2, FUN = mean)
  A_med <- apply(A_mat, MARGIN = 2, FUN = mean)
  Gamma_med <- apply(Gamma_mat, MARGIN = 2, FUN = mean)
  Sigma_med <- apply(Sigma_mat, MARGIN = 2, FUN = mean)
  Nu_med <- apply(Nu_mat, MARGIN = 2, FUN = mean)
  W_med <- apply(W_mat, MARGIN = 2, FUN = mean)
  H_med <- apply(H_mat, MARGIN = 2, FUN = mean)
  H0_med <- apply(H0_mat, MARGIN = 2, FUN = mean)

  inits$samples <- samples # Reset
  inits$burnin <- burnin # Reset
  inits$thin <- thin # Reset

  # Commom parameters
  inits$A0 <- a0toA(A_med, K = K)
  inits$B0 <- matrix(B_med, nrow = K)
  inits$sigma_h <- diag(Sigma_med)
  inits$h <- matrix(H_med, nrow = K)
  inits$h0 <- H0_med
  #inits$sigma <- Redundant


  str(inits)

  if (prior$SV == TRUE){
    Start = Sys.time()
    if (dist == "Gaussian") {

      ChibLLP_chain <- ChibLLP(Chain, inits, ndraws = 1000, numCores = numCores)
      RhpcBLASctl::blas_set_num_threads(numCores)

      cat("Pi( B* | y, A*, Sigma*) \n ")
      P1 <- Chib.Gaussian.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                                samples = samples, burnin = burnin, thin = thin,
                                fix_B = FALSE, fix_A = TRUE, fix_Sigma = TRUE,
                                cal_B = TRUE, cal_A = FALSE, cal_Sigma = FALSE)
      cat("Pi( A* | y, Sigma*) \n ")
      P2 <- Chib.Gaussian.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                                samples = samples, burnin = burnin, thin = thin,
                                fix_B = FALSE, fix_A = FALSE, fix_Sigma = TRUE,
                                cal_B = FALSE, cal_A = TRUE, cal_Sigma = FALSE)
      cat("Pi( Sigma* | y) \n ")
      P3 <- Chib.Gaussian.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
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
      ChibLLP_chain <- ChibLLP(Chain, inits, ndraws = 1000, numCores = numCores)
      RhpcBLASctl::blas_set_num_threads(numCores)

      cat("Pi( B* | y, A*, Sigma*, Nu*) \n ")
      P1 <- Chib.Student.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                               samples = samples, burnin = burnin, thin = thin,
                               fix_B = FALSE, fix_A = TRUE, fix_Sigma = TRUE, fix_Nu = TRUE,
                               cal_B = TRUE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE)


      cat("Pi( A* | y, Sigma*, Nu*) \n ")
      P2 <- Chib.Student.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                               samples = samples, burnin = burnin, thin = thin,
                               fix_B = FALSE, fix_A = FALSE, fix_Sigma = TRUE, fix_Nu = TRUE,
                               cal_B = FALSE, cal_A = TRUE, cal_Sigma = FALSE, cal_Nu = FALSE)
      cat(" Pi( Sigma* | y, Nu*) \n ")
      P3 <- Chib.Student.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                               samples = samples, burnin = burnin, thin = thin,
                               fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                               cal_B = FALSE, cal_A = FALSE, cal_Sigma = TRUE, cal_Nu = FALSE)

      cat(" Pi( Nu* | y) --> Nominator \n ")
      P4 <- Chib.Student.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                               samples = samples, burnin = burnin, thin = thin,
                               fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                               cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      cat("Pi( Nu* | y) --> Denominator \n ")
      P5 <- Chib.Student.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                               samples = samples, burnin = burnin, thin = thin,
                               fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                               cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      CML_Chain <- ChibLLP_chain$LL - logmeanexp(P1) - logmeanexp(P2) - logmeanexp(P3) - logmeanexp(P4) + logmeanexp(P5) # P5 Denominator

    }
    if (dist == "Skew.Student") {
      inits$w <- W_med
      ChibLLP_chain <- ChibLLP(Chain, inits, ndraws = 1000, numCores = numCores)
      RhpcBLASctl::blas_set_num_threads(numCores)

      cat("Pi( B* | y, A*, Sigma*, Nu*) \n ")
      P1 <- Chib.Skew.Student.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                                    samples = samples, burnin = burnin, thin = thin,
                                    fix_B = FALSE, fix_A = TRUE, fix_Sigma = TRUE, fix_Nu = TRUE,
                                    cal_B = TRUE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE)


      cat("Pi( A* | y, Sigma*, Nu*) \n ")
      P2 <- Chib.Skew.Student.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                                    samples = samples, burnin = burnin, thin = thin,
                                    fix_B = FALSE, fix_A = FALSE, fix_Sigma = TRUE, fix_Nu = TRUE,
                                    cal_B = FALSE, cal_A = TRUE, cal_Sigma = FALSE, cal_Nu = FALSE)
      cat(" Pi( Sigma* | y, Nu*) \n ")
      P3 <- Chib.Skew.Student.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                                    samples = samples, burnin = burnin, thin = thin,
                                    fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                                    cal_B = FALSE, cal_A = FALSE, cal_Sigma = TRUE, cal_Nu = FALSE)

      cat(" Pi( Nu* | y) --> Nominator \n ")
      P4 <- Chib.Skew.Student.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                                    samples = samples, burnin = burnin, thin = thin,
                                    fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                                    cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      cat("Pi( Nu* | y) --> Denominator \n ")
      P5 <- Chib.Skew.Student.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                                    samples = samples, burnin = burnin, thin = thin,
                                    fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                                    cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      CML_Chain <- ChibLLP_chain$LL - logmeanexp(P1) - logmeanexp(P2) - logmeanexp(P3) - logmeanexp(P4) + logmeanexp(P5) # P5 Denominator


    }
    if (dist == "MT") {
      inits$w <- matrix(W_med, nrow = K)
      ChibLLP_chain <- ChibLLP(Chain, inits, ndraws = 1000, numCores = numCores)
      RhpcBLASctl::blas_set_num_threads(numCores)


      cat("Pi( B* | y, A*, Sigma*, Nu*) \n ")
      P1 <- Chib.MT.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = TRUE, fix_Sigma = TRUE, fix_Nu = TRUE,
                          cal_B = TRUE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE)


      cat("Pi( A* | y, Sigma*, Nu*) \n ")
      P2 <- Chib.MT.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = TRUE, fix_Nu = TRUE,
                          cal_B = FALSE, cal_A = TRUE, cal_Sigma = FALSE, cal_Nu = FALSE)
      cat("Pi( Sigma* | y, Nu*) \n ")
      P3 <- Chib.MT.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                          cal_B = FALSE, cal_A = FALSE, cal_Sigma = TRUE, cal_Nu = FALSE)

      cat("Pi( Nu* | y) --> Nominator \n ")#
      P4 <- Chib.MT.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                          cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      cat("Pi( Nu* | y) --> Denominator \n ")
      P5 <- Chib.MT.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                          cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      CML_Chain <- ChibLLP_chain$LL - logmeanexp(P1) - logmeanexp(P2) - logmeanexp(P3) - sum(logmeanexp(P4)) + sum(logmeanexp(P5) ) # P5 Denominator

    }
    if (dist == "MST") {
      inits$w <- matrix(W_med, nrow = K)
      ChibLLP_chain <- ChibLLP(Chain, inits, ndraws = 1000, numCores = numCores)
      RhpcBLASctl::blas_set_num_threads(numCores)


      cat("Pi( B* | y, A*, Sigma*, Nu*) \n ")
      P1 <- Chib.MST.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                           samples = samples, burnin = burnin, thin = thin,
                           fix_B = FALSE, fix_A = TRUE, fix_Sigma = TRUE, fix_Nu = TRUE,
                           cal_B = TRUE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE)


      cat("Pi( A* | y, Sigma*, Nu*) \n ")
      P2 <- Chib.MST.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                           samples = samples, burnin = burnin, thin = thin,
                           fix_B = FALSE, fix_A = FALSE, fix_Sigma = TRUE, fix_Nu = TRUE,
                           cal_B = FALSE, cal_A = TRUE, cal_Sigma = FALSE, cal_Nu = FALSE)
      cat("Pi( Sigma* | y, Nu*) \n ")
      P3 <- Chib.MST.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                           samples = samples, burnin = burnin, thin = thin,
                           fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                           cal_B = FALSE, cal_A = FALSE, cal_Sigma = TRUE, cal_Nu = FALSE)

      cat("Pi( Nu* | y) --> Nominator \n ")#
      P4 <- Chib.MST.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                           samples = samples, burnin = burnin, thin = thin,
                           fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                           cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      cat("Pi( Nu* | y) --> Denominator \n ")
      P5 <- Chib.MST.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                           samples = samples, burnin = burnin, thin = thin,
                           fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                           cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      CML_Chain <- ChibLLP_chain$LL - logmeanexp(P1) - logmeanexp(P2) - logmeanexp(P3) - sum(logmeanexp(P4)) + sum(logmeanexp(P5) ) # P5 Denominator

    }
    if (dist == "OT") {
      inits$w <- matrix(W_med, nrow = K)
      ChibLLP_chain <- ChibLLP(Chain, inits, ndraws = 1000, numCores = numCores)
      RhpcBLASctl::blas_set_num_threads(numCores)


      cat("Pi( B* | y, A*, Sigma*, Nu*) \n ")
      P1 <- Chib.OT.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = TRUE, fix_Sigma = TRUE, fix_Nu = TRUE,
                          cal_B = TRUE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE)


      cat("Pi( A* | y, Sigma*, Nu*) \n ")
      P2 <- Chib.OT.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = TRUE, fix_Nu = TRUE,
                          cal_B = FALSE, cal_A = TRUE, cal_Sigma = FALSE, cal_Nu = FALSE)
      cat("Pi( Sigma* | y, Nu*) \n ")
      P3 <- Chib.OT.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                          cal_B = FALSE, cal_A = FALSE, cal_Sigma = TRUE, cal_Nu = FALSE)

      cat("Pi( Nu* | y) --> Nominator \n ")#
      P4 <- Chib.OT.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                          cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      cat("Pi( Nu* | y) --> Denominator \n ")
      P5 <- Chib.OT.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                          samples = samples, burnin = burnin, thin = thin,
                          fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                          cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      CML_Chain <- ChibLLP_chain$LL - logmeanexp(P1) - logmeanexp(P2) - logmeanexp(P3) - sum(logmeanexp(P4)) + sum(logmeanexp(P5) ) # P5 Denominator

    }
    if (dist == "OST") {
      inits$w <- matrix(W_med, nrow = K)
      ChibLLP_chain <- ChibLLP(Chain, inits, ndraws = 1000, numCores = numCores)
      RhpcBLASctl::blas_set_num_threads(numCores)


      cat("Pi( B* | y, A*, Sigma*, Nu*) \n ")
      P1 <- Chib.OST.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                           samples = samples, burnin = burnin, thin = thin,
                           fix_B = FALSE, fix_A = TRUE, fix_Sigma = TRUE, fix_Nu = TRUE,
                           cal_B = TRUE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = FALSE)


      cat("Pi( A* | y, Sigma*, Nu*) \n ")
      P2 <- Chib.OST.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                           samples = samples, burnin = burnin, thin = thin,
                           fix_B = FALSE, fix_A = FALSE, fix_Sigma = TRUE, fix_Nu = TRUE,
                           cal_B = FALSE, cal_A = TRUE, cal_Sigma = FALSE, cal_Nu = FALSE)
      cat("Pi( Sigma* | y, Nu*) \n ")
      P3 <- Chib.OST.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                           samples = samples, burnin = burnin, thin = thin,
                           fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = TRUE,
                           cal_B = FALSE, cal_A = FALSE, cal_Sigma = TRUE, cal_Nu = FALSE)

      cat("Pi( Nu* | y) --> Nominator \n ")#
      P4 <- Chib.OST.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
                           samples = samples, burnin = burnin, thin = thin,
                           fix_B = FALSE, fix_A = FALSE, fix_Sigma = FALSE, fix_Nu = FALSE,
                           cal_B = FALSE, cal_A = FALSE, cal_Sigma = FALSE, cal_Nu = TRUE)

      cat("Pi( Nu* | y) --> Denominator \n ")
      P5 <- Chib.OST.SV(y, K = K, p = p, y0 = y0, prior = prior, inits = inits,
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
                esttime = elapsedTime)
    class(out) <- c("ChibMLM")

    } else {
      warning("prior$SV is FALSE")
    }
    return(out)
}
#' @export
Chib.Gaussian.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL,
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
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Gaussian", SV = TRUE)
  }
  # prior B
  b_prior = prior$b_prior
  V_b_prior = prior$V_b_prior
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior h0_mean
  h0_mean <- 2 * log(prior$sigma)
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
  h <- inits$h
  sigma_h <- inits$sigma_h

  V_b_prior_inv <- solve(V_b_prior)
  # new precompute
  theta.prior.precmean <- V_b_prior_inv %*% b_prior

  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K*t_max,
                     ncol = (samples - burnin)%/% thin)
  lpost <- rep(0, (samples - burnin)%/% thin)
  for (j in c(1:samples)){

    # Sample B
    if(!fix_B) {
      wt <- as.vector(exp(-h/2))
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

      wt <- as.vector(exp(-h/2))
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
    # Sample vol
    ytilde <- A%*% (yt - B %*% xt)
    aux <- sample_h_mod(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    sqrtvol <- aux$sigt
    h0 <- as.numeric(aux$h0)

    if(!fix_Sigma) {
      if (K>1) {
        sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
      } else {
        sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
      }

      sigma_h <- diag(mapply( GIGrvg::rgig, n = 1, lambda = - (t_max - 1)*0.5, chi = sse_2,
                              psi = 1/prior$sigma_S0 ), nrow = K)

    }

    if (cal_Sigma & (j > burnin) & (j %% thin == 0)){

      if (K>1) {
        sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
      } else {
        sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
      }


      lpost[(j - burnin) %/% thin] <- sum(mapply(GIGrvg::dgig, x = diag(inits$sigma_h), lambda = - (t_max - 1)*0.5,
                                                 chi = sse_2, psi = 1/prior$sigma_S0, log = T))
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
          a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sqrtvol[i,],
                                                       xsub = matrix(u_neg[1:(i-1),] / sqrtvol[i,], nrow = i-1),
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

          ysub = u_std[i,] / sqrtvol[i,]
          xsub = matrix(u_neg[1:(i-1),] / sqrtvol[i,], nrow = i-1)

          n = nrow(xsub)

          V_a_sub_inv = solve(V_a_sub)

          V_a_post_inv <- V_a_sub_inv + xsub %*% t(xsub)
          V_a_post <- solve(V_a_post_inv)
          a_post <- V_a_post %*% ( V_a_sub_inv %*% a_sub + xsub %*% ysub)

          lpost[(j - burnin) %/% thin] <- lpost[(j - burnin) %/% thin] + mvnfast::dmvn(sub_A_star, mu = a_post, sigma = V_a_post, log = T)
        }
      }
    }
    if ((j > inits$burnin) & (j %% inits$thin == 0)){
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, diag(sigma_h), h0, as.numeric(h))
    }

    if (j %% 1000 == 0) { cat(" Iteration ", j, " \n")}
  }
  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("sigma_h",c(1:K), sep = ""),
                        paste("lh0",c(1:K), sep = ""),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K))
  )
  return(lpost)
}

###########################################################################
#' @export
Chib.Student.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL,
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
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Student", SV = TRUE)
  }
  # prior B
  b_prior = prior$b_prior
  V_b_prior = prior$V_b_prior
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior h0_mean
  h0_mean <- 2 * log(prior$sigma)
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
  B <- inits$B0
  b_sample <- as.numeric(B)
  a_sample <- as.numeric(A[lower.tri(A)])
  h <- inits$h
  sigma_h <- inits$sigma_h

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
  mcmc <-  matrix(NA, nrow = m*K + 0.5*K*(K-1) + 1 + K + K + K*t_max + t_max,
                  ncol = (samples - burnin)%/% thin)
  lpost <- rep(0, (samples - burnin)%/% thin)
  for (j in c(1:samples)){


    # Sample B
    if(!fix_B) {
      wt <- as.vector(1/(exp(h/2)*w_sqrt))
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

      wt <- as.vector(1/(exp(h/2)*w_sqrt))
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

    # Sample vol
    ytilde <- A%*% ((yt - B %*% xt)/w_sqrt)
    aux <- sample_h_mod(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    sqrtvol <- aux$sigt
    h0 <- as.numeric(aux$h0)

    if(!fix_Sigma) {
      if (K>1) {
        sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
      } else {
        sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
      }

      sigma_h <- diag(mapply( GIGrvg::rgig, n = 1, lambda = - (t_max - 1)*0.5, chi = sse_2,
                              psi = 1/prior$sigma_S0 ), nrow = K)

    }

    if (cal_Sigma & (j > burnin) & (j %% thin == 0)){
      if (K>1) {
        sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
      } else {
        sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
      }


      lpost[(j - burnin) %/% thin] <- sum(mapply(GIGrvg::dgig, x = diag(inits$sigma_h), lambda = - (t_max - 1)*0.5,
                                                 chi = sse_2, psi = 1/prior$sigma_S0, log = T))
    }


    # Sample A0
    if(!fix_A) {
      u_std <- (yt - B %*%xt)/ w_sqrt # change from Gaussian
      u_neg <- - u_std
      a_sample <- rep(0, K * (K - 1) /2)
      if (K > 1) {
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- a_prior[id_start:id_end]
          V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
          a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sqrtvol[i,],
                                                       xsub = matrix(u_neg[1:(i-1),] / sqrtvol[i,], nrow = i-1),
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
      u_std <- (yt - B %*%xt)/ w_sqrt # change from Gaussian
      u_neg <- - u_std
      A_star <- inits$A0[lower.tri(A)]
      if (K > 1) {
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- prior$a_prior[id_start:id_end]
          V_a_sub <- prior$V_a_prior[id_start:id_end, id_start:id_end]
          sub_A_star <- A_star[id_start:id_end]

          ysub = u_std[i,] / sqrtvol[i,]
          xsub = matrix(u_neg[1:(i-1),] / sqrtvol[i,], nrow = i-1)

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
    w_sample <- rinvgamma( n = t_max, shape = (nu+K)*0.5, rate = 0.5*( nu + colSums((u^2)/exp(h)) ) )
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
        # Move from nu_star to nu, using truncated proposal
        while (TRUE){
          nu_temp = inits$nu + exp(logsigma_nu)*rnorm(1)
          if (nu_temp > 4 && nu_temp < 100){
            break
          }
        }
          num_mh = dgamma(nu_temp, shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w_sample, shape = nu_temp*0.5, rate = nu_temp*0.5, log = T)) -
            log( pnorm(100, inits$nu, exp(logsigma_nu)) - pnorm(4, inits$nu, exp(logsigma_nu)) ) # truncation proposal

          denum_mh = dgamma(inits$nu, shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w_sample, shape = inits$nu*0.5, rate = inits$nu*0.5, log = T)) -
          log( pnorm(100, nu_temp,  exp(logsigma_nu)) - pnorm(4, nu_temp, exp(logsigma_nu)) )  # truncation proposal

          alpha = min(1, exp(num_mh - denum_mh))
          temp = runif(1)
          if (alpha > temp){
            acount_nu = acount_nu + 1
          }
          lpost[(j - burnin) %/% thin] <- log(alpha) # Take a negative for denominator


      } else {
        # Nominator
        q_nu_nustar <- dnorm(inits$nu, mean = nu, sd = exp(logsigma_nu) ) /
          ( pnorm(100, nu,  exp(logsigma_nu)) - pnorm(4, nu, exp(logsigma_nu)) ) # Move nu to inits$nu

        num_mh = dgamma(inits$nu, shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w_sample, shape = inits$nu*0.5, rate = inits$nu*0.5, log = T))
        denum_mh = dgamma(nu, shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w_sample, shape = nu*0.5, rate = nu*0.5, log = T))
        alpha = min(1, exp(num_mh - denum_mh))
        lpost[(j - burnin) %/% thin] <- log( alpha * q_nu_nustar)
      }
    }

    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, nu, diag(sigma_h), h0, as.numeric(h), as.numeric(w_sample))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", round(acount_nu/j,2)," ", round(nu,2)," \n")
      acount_w <- rep(0,t_max)
    }
  }
  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("nu"),
                        paste("sigma_h",c(1:K), sep = ""),
                        paste("lh0",c(1:K), sep = ""),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        paste("w",c(1:t_max), sep = ""))

  return(lpost)
}

#############################################################################################
#' @export
Chib.Skew.Student.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL,
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
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Skew.Student", SV = TRUE)
  }
  # prior B
  b_prior = prior$b_prior
  V_b_prior = prior$V_b_prior
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior h0_mean
  h0_mean <- 2 * log(prior$sigma)
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
  B <- inits$B0
  b_sample <- as.numeric(B)
  a_sample <- as.numeric(A[lower.tri(A)])
  h <- inits$h
  sigma_h <- inits$sigma_h

  V_b_prior_inv <- solve(V_b_prior)

  nu <- inits$nu
  mu.xi <- nu/( nu - 2 )
  logsigma_nu <- inits$logsigma_nu
  acount_nu <- 0
  acount_w <- rep(0, t_max)
  gamma <- inits$gamma
  D <- diag(gamma, nrow = K)

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
  mcmc <-  matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + 1 + K + K + K*t_max + t_max,
                  ncol = (samples - burnin)%/% thin)
  lpost <- rep(0, (samples - burnin)%/% thin)
  for (j in c(1:samples)){

    # Sample B and gamma
    if(!fix_B) {
      wt <- as.vector(1/(exp(h/2)*w_sqrt))
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

      wt <- as.vector(1/(exp(h/2)*w_sqrt))
      y.tilde <- as.vector( A %*% yt ) * wt
      x.tilde <- kronecker( cbind( t(xt), w_sample - mu.xi ), A ) * wt
      theta.prec.chol <- chol( theta.prior.prec + crossprod(x.tilde) )
      b_star <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T ))
      #mvnfast::dmvn(B_star, mu = b_star, sigma = (solve(theta.prior.prec + crossprod(x.tilde) )), log = T)
      lpost[(j - burnin) %/% thin] <- - length(B_star) * 0.5 * log(2*pi) + sum(log(diag(theta.prec.chol))) -
        0.5 * t(B_star - b_star) %*% (theta.prior.prec + crossprod(x.tilde)) %*% (B_star - b_star)
    }
    # Sample vol
    ytilde <- A%*% ((yt - B %*% xt  - D%*% (w-mu.xi) )/w_sqrt)
    aux <- sample_h_mod(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    sqrtvol <- aux$sigt
    h0 <- as.numeric(aux$h0)

    if(!fix_Sigma) {
      if (K>1) {
        sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
      } else {
        sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
      }

      sigma_h <- diag(mapply( GIGrvg::rgig, n = 1, lambda = - (t_max - 1)*0.5, chi = sse_2,
                              psi = 1/prior$sigma_S0 ), nrow = K)

    }

    if (cal_Sigma & (j > burnin) & (j %% thin == 0)){
      if (K>1) {
        sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
      } else {
        sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
      }


      lpost[(j - burnin) %/% thin] <- sum(mapply(GIGrvg::dgig, x = diag(inits$sigma_h), lambda = - (t_max - 1)*0.5,
                                                 chi = sse_2, psi = 1/prior$sigma_S0, log = T))
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
          a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sqrtvol[i,],
                                                       xsub = matrix(u_neg[1:(i-1),] / sqrtvol[i,], nrow = i-1),
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

          ysub = u_std[i,] / sqrtvol[i,]
          xsub = matrix(u_neg[1:(i-1),] / sqrtvol[i,], nrow = i-1)

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
    q2 <- colSums( ( ( A %*% ( yt - B%*%xt + mu.xi * gamma) ) / sqrtvol )^2 )
    p2 <- colSums( ( c( A %*% gamma ) / sqrtvol )^2 )
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
        # Move from nu_star to nu, using truncated proposal
        while (TRUE){
          nu_temp = inits$nu + exp(logsigma_nu)*rnorm(1)
          if (nu_temp > 4 && nu_temp < 100){
            break
          }
        }

          num_mh = dgamma(nu_temp, shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w_sample, shape = nu_temp*0.5, rate = nu_temp*0.5, log = T)) -
            log( pnorm(100, inits$nu, exp(logsigma_nu)) - pnorm(4, inits$nu, exp(logsigma_nu)) ) # truncation proposal

          denum_mh = dgamma(inits$nu, shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w_sample, shape = inits$nu*0.5, rate = inits$nu*0.5, log = T)) -
            log( pnorm(100, nu_temp,  exp(logsigma_nu)) - pnorm(4, nu_temp, exp(logsigma_nu)) )  # truncation proposal
          alpha = min(1, exp(num_mh - denum_mh))
          temp = runif(1)
          if (alpha > temp){
            acount_nu = acount_nu + 1
          }
          lpost[(j - burnin) %/% thin] <- log(alpha) # Take a negative for denominator


      } else {
        # Nominator
        q_nu_nustar <- dnorm(inits$nu, mean = nu, sd = exp(logsigma_nu) ) /
          ( pnorm(100, nu,  exp(logsigma_nu)) - pnorm(4, nu, exp(logsigma_nu)) ) # Move nu to inits$nu
        num_mh = dgamma(inits$nu, shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w_sample, shape = inits$nu*0.5, rate = inits$nu*0.5, log = T))
        denum_mh = dgamma(nu, shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w_sample, shape = nu*0.5, rate = nu*0.5, log = T))
        alpha = min(1, exp(num_mh - denum_mh))
        lpost[(j - burnin) %/% thin] <- log( alpha * q_nu_nustar)
      }
    }

    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, gamma, nu, diag(sigma_h), h0, as.numeric(h), as.numeric(w_sample))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", round(acount_nu/j,2)," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2) , " ", gamma ,  " \n")
      acount_w <- rep(0,t_max)
    }
  }
  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("gamma",c(1:K), sep = ""),
                        paste("nu"),
                        paste("sigma_h",c(1:K), sep = ""),
                        paste("lh0",c(1:K), sep = ""),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        paste("w",c(1:t_max), sep = ""))

  return(lpost)
}

###########################################################################
#' @export
Chib.MT.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL,
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
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "MT", SV = TRUE)
  }
  # prior B
  b_prior = prior$b_prior
  V_b_prior = prior$V_b_prior
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior h0_mean
  h0_mean <- 2 * log(prior$sigma)
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
  B <- inits$B0
  b_sample <- as.numeric(B)
  a_sample <- as.numeric(A[lower.tri(A)])
  h <- inits$h
  sigma_h <- inits$sigma_h

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
  lpost <- rep(0, (samples - burnin)%/% thin)
  repeat_time <- 1
  if (cal_Nu) {
    lpost <- matrix(0,nrow = (samples - burnin)%/% thin, ncol = K)
    repeat_time <- K
  }

  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K*t_max + K*t_max,
                 ncol = (samples - burnin)%/% thin)

  for (nu_id in c(1:repeat_time)){
    # Multi degrees of freedom
    nu <- inits$nu
    logsigma_nu <- inits$logsigma_nu
    acount_nu <- rep(0,K)
    acount_w <- rep(0, t_max)
    # Init w as Gaussian
    w <- inits$w
    w_sqrt <- sqrt(w)

    for (j in c(1:samples)){

      # Sample B
      if(!fix_B) {
        wt <- as.vector(exp(-h/2))
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

        wt <- as.vector(exp(-h/2))
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

      # Sample vol
      ytilde <- A%*% ((yt - B %*%xt)/ w_sqrt)
      aux <- sample_h_mod(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
      h <- aux$Sigtdraw
      sqrtvol <- aux$sigt
      h0 <- as.numeric(aux$h0)

      if(!fix_Sigma) {
        if (K>1) {
          sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
        } else {
          sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
        }

        sigma_h <- diag(mapply( GIGrvg::rgig, n = 1, lambda = - (t_max - 1)*0.5, chi = sse_2,
                                psi = 1/prior$sigma_S0 ), nrow = K)

      }

      if (cal_Sigma & (j > burnin) & (j %% thin == 0)){
        if (K>1) {
          sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
        } else {
          sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
        }


        lpost[(j - burnin) %/% thin] <- sum(mapply(GIGrvg::dgig, x = diag(inits$sigma_h), lambda = - (t_max - 1)*0.5,
                                                   chi = sse_2, psi = 1/prior$sigma_S0, log = T))
      }


      # Sample A0
      if(!fix_A) {
        u_std <- (yt - B %*%xt) / w_sqrt # change from Gaussian
        u_neg <- - u_std
        a_sample <- rep(0, K * (K - 1) /2)
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- a_prior[id_start:id_end]
          V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
          a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sqrtvol[i,],
                                                       xsub = matrix(u_neg[1:(i-1),] / sqrtvol[i,], nrow = i-1),
                                                       a_sub = a_sub,
                                                       V_a_sub = V_a_sub)
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

            ysub = u_std[i,] / sqrtvol[i,]
            xsub = matrix(u_neg[1:(i-1),] / sqrtvol[i,], nrow = i-1)

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
      u_proposal <- A %*% (yt - B %*%xt)
      a_target <- (nu*0.5 + 1*0.5) * 0.75
      b_target <- (nu*0.5 + 0.5 * u_proposal^2 / sqrtvol^2)*0.75  # adjust by 0.75
      w_temp <- matrix( rinvgamma( n = K*t_max, shape = a_target, rate = b_target ), K, t_max ) # rinvgamma recycles arguments
      w_temp_sqrt <- sqrt(w_temp)

      A.w.u.s.prop <- ( A %*% ( u/w_temp_sqrt ) )/sqrtvol
      num_mh <- colSums(dinvgamma(w_temp, shape = nu*0.5, rate = nu*0.5, log = T)) - #prior
                0.5*( colSums(log(w_temp)) + colSums(A.w.u.s.prop^2) ) - # posterior
                colSums( dinvgamma(w_temp, shape = a_target, rate = b_target, log = T)) # proposal
      A.w.u.s.curr <- ( A %*% ( u/w_sqrt ) )/sqrtvol
      denum_mh <- colSums(dinvgamma(w, shape = nu*0.5, rate = nu*0.5, log = T)) - #prior
                  0.5*( colSums(log(w)) + colSums(A.w.u.s.curr^2) ) - # posterior
                  colSums( dinvgamma(w, shape = a_target, rate = b_target, log = T)) # proposal
      acc <- (num_mh-denum_mh) > log(runif(t_max))
      w[,acc] <- w_temp[,acc]
      w_sqrt[,acc] <- w_temp_sqrt[,acc]
      acount_w[acc] <- acount_w[acc] + 1


      # Sample nu nominator
      if(!fix_Nu) {
        # if not Chib calculation of nu, do it as normal
        if (!cal_Nu){
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
        } else {

          # change nu up to nu_id, fix all the rest
          nu_temp = nu + exp(logsigma_nu)*rnorm(K)
          for (k in c(1:nu_id)){
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

      }
      # Sample nu denominator
      if (cal_Nu & fix_Nu){
        # change nu up to nu_id, fix all the rest
        if (nu_id > 1){
          nu_temp = nu + exp(logsigma_nu)*rnorm(K)
          for (k in c(1:(nu_id-1) )){
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


      }


      # Now calculate the Chib estimation of p(nu[nu_id] | rest)
      if (cal_Nu & (j > burnin) & (j %% thin == 0)){
        # Denominator
        if (fix_Nu) {

          # nu up nu_id only
          k <- nu_id
          # Move from nu_star to nu, using truncated proposal
          while (TRUE){
            nu_temp = inits$nu[k] + exp(logsigma_nu[k])*rnorm(1)
            if (nu_temp > 4 && nu_temp < 100){
              break
            }
          }
            num_mh = dgamma(nu_temp, shape = nu_gam_a, rate = nu_gam_b, log = T) +
              sum(dinvgamma(w[k,], shape = nu_temp*0.5, rate = nu_temp*0.5, log = T)) -
              log( pnorm(100, inits$nu[k], exp(logsigma_nu[k])) - pnorm(4, inits$nu[k], exp(logsigma_nu[k])) ) # truncation proposal

            denum_mh = dgamma(inits$nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
              sum(dinvgamma(w[k,], shape = inits$nu[k]*0.5, rate = inits$nu[k]*0.5, log = T)) -
              log( pnorm(100, nu_temp,  exp(logsigma_nu[k])) - pnorm(4, nu_temp, exp(logsigma_nu[k])) )  # truncation proposal

            alpha = min(1, exp(num_mh - denum_mh))
            temp = runif(1)
            if (alpha > temp){
              acount_nu[k] = acount_nu[k] + 1
            }
            lpost[(j - burnin) %/% thin,k] <- log(alpha) # Take a negative for denominator

        } else {
          # Nominator
          k <- nu_id
          q_nu_nustar <- dnorm(inits$nu[k], mean = nu[k], sd = exp(logsigma_nu[k]) ) /
            ( pnorm(100, nu[k],  exp(logsigma_nu[k])) - pnorm(4, nu[k], exp(logsigma_nu[k])) ) # Move nu to inits$nu
          num_mh = dgamma(inits$nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = inits$nu[k]*0.5, rate = inits$nu[k]*0.5, log = T))
          denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))
          alpha = min(1, exp(num_mh - denum_mh))
          lpost[(j - burnin) %/% thin,k] <- log( alpha * q_nu_nustar)
        }

      }


      if ((j > inits$burnin) & (j %% inits$thin == 0))
        mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, nu, diag(sigma_h), h0, as.numeric(h), as.numeric(w))

      if (j %% 1000 == 0) {
        cat(" Iteration ", j, " ", round(acount_nu/j,2)," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2), " \n")
        acount_w <- rep(0,t_max)
      }
      }
  }
  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("nu",c(1:K), sep = ""),
                        paste("sigma_h",c(1:K), sep = ""),
                        paste("lh0",c(1:K), sep = ""),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)))

  return(lpost)
}

#############################################################################################
#' @export
Chib.MST.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL,
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
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "MST", SV = TRUE)
  }
  # prior B
  b_prior = prior$b_prior
  V_b_prior = prior$V_b_prior
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior h0_mean
  h0_mean <- 2 * log(prior$sigma)
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
  B <- inits$B0
  b_sample <- as.numeric(B)
  a_sample <- as.numeric(A[lower.tri(A)])
  h <- inits$h
  sigma_h <- inits$sigma_h

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

  # new precompute
  i <- K*m
  theta.prior.prec = matrix(0,i+K,i+K)
  theta.prior.prec[1:i,1:i] <- V_b_prior_inv
  theta.prior.prec[i+(1:K),i+(1:K)] <- solve( V_gamma_prior )
  theta.prior.precmean <- theta.prior.prec %*% c( b_prior, gamma_prior )

  # Output
  lpost <- rep(0, (samples - burnin)%/% thin)
  repeat_time <- 1
  if (cal_Nu) {
    lpost <- matrix(0,nrow = (samples - burnin)%/% thin, ncol = K)
    repeat_time <- K
  }
  mcmc <-  matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K + K*t_max + K*t_max,
                  ncol = (samples - burnin)%/% thin)
  if (cal_Nu) lpost <- matrix(0,nrow = (samples - burnin)%/% thin, ncol = K)
  for (nu_id in c(1:repeat_time)){
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

    for (j in c(1:samples)){

      # Sample B and gamma
      if(!fix_B) {
        wt <- as.vector(exp(-h/2))
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

        wt <- as.vector(exp(-h/2))
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

      # Sample vol
      ytilde <- A%*% ((yt - B %*%xt - D %*% ( w - mu.xi ))/ w_sqrt)
      aux <- sample_h_mod(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
      h <- aux$Sigtdraw
      sqrtvol <- aux$sigt
      h0 <- as.numeric(aux$h0)

      if(!fix_Sigma) {
        if (K>1) {
          sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
        } else {
          sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
        }

        sigma_h <- diag(mapply( GIGrvg::rgig, n = 1, lambda = - (t_max - 1)*0.5, chi = sse_2,
                                psi = 1/prior$sigma_S0 ), nrow = K)

      }

      if (cal_Sigma & (j > burnin) & (j %% thin == 0)){
        if (K>1) {
          sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
        } else {
          sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
        }


        lpost[(j - burnin) %/% thin] <- sum(mapply(GIGrvg::dgig, x = diag(inits$sigma_h), lambda = - (t_max - 1)*0.5,
                                                   chi = sse_2, psi = 1/prior$sigma_S0, log = T))
      }

      # Sample A0
      if(!fix_A) {
        u_std <- (yt - B %*% xt - D%*% ( w - mu.xi )) / w_sqrt # change from Gaussian
        u_neg <- - u_std
        a_sample <- rep(0, K * (K - 1) /2)
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- a_prior[id_start:id_end]
          V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
          a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sqrtvol[i,],
                                                       xsub = matrix(u_neg[1:(i-1),] / sqrtvol[i,], nrow = i-1),
                                                       a_sub = a_sub,
                                                       V_a_sub = V_a_sub)
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

            ysub = u_std[i,] / sqrtvol[i,]
            xsub = matrix(u_neg[1:(i-1),] / sqrtvol[i,], nrow = i-1)

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
      u <- (yt - B %*%xt + mu.xi * gamma )
      u_proposal <- A %*% (yt - B %*%xt) + mu.xi * gamma
      a_target <- (nu*0.5 + 1*0.5) * 0.75
      b_target <- (nu*0.5 + 0.5 * u_proposal^2 / sqrtvol^2)*0.75  # adjust by 0.75
      w_temp <- matrix( rinvgamma( n = K*t_max, shape = a_target, rate = b_target ), K, t_max ) # rinvgamma recycles arguments
      log_w_temp <- log(w_temp)
      w_temp_sqrt <- sqrt(w_temp)

      A.w.u.s.prop <- ( A %*% ( u/w_temp_sqrt - w_temp_sqrt*gamma ) )/sqrtvol
      num_mh <- colSums(dinvgamma(w_temp, shape = nu*0.5, rate = nu*0.5, log = T)) - #prior
        0.5*( colSums(log(w_temp)) + colSums(A.w.u.s.prop^2) ) - # posterior
        colSums( dinvgamma(w_temp, shape = a_target, rate = b_target, log = T)) # proposal
      A.w.u.s.curr <- ( A %*% ( u/w_sqrt - w_sqrt*gamma ) )/sqrtvol
      denum_mh <- colSums(dinvgamma(w, shape = nu*0.5, rate = nu*0.5, log = T)) - #prior
        0.5*( colSums(log(w)) + colSums(A.w.u.s.curr^2) ) - # posterior
        colSums( dinvgamma(w, shape = a_target, rate = b_target, log = T)) # proposal
      acc <- (num_mh-denum_mh) > log(runif(t_max))
      w[,acc] <- w_temp[,acc]
      w_sqrt[,acc] <- w_temp_sqrt[,acc]
      acount_w[acc] <- acount_w[acc] + 1

      # Sample nu nominator
      if(!fix_Nu) {
        # if not Chib calculation of nu, do it as normal
        if (!cal_Nu){
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

        } else {

          # change nu up to nu_id, fix all the rest
          nu_temp = nu + exp(logsigma_nu)*rnorm(K)
          for (k in c(1:nu_id)){
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

      }
      # Sample nu denominator
      if (cal_Nu & fix_Nu){
        # change nu up to nu_id, fix all the rest
        if (nu_id > 1){
          nu_temp = nu + exp(logsigma_nu)*rnorm(K)
          for (k in c(1:(nu_id-1) )){
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


      }


      # Now calculate the Chib estimation of p(nu[nu_id] | rest)
      if (cal_Nu & (j > burnin) & (j %% thin == 0)){
        # Denominator
        if (fix_Nu) {

          # nu up nu_id only
          k <- nu_id
          # Move from nu_star to nu, using truncated proposal
          while (TRUE){
            nu_temp = inits$nu[k] + exp(logsigma_nu[k])*rnorm(1)
            if (nu_temp > 4 && nu_temp < 100){
              break
            }
          }
            num_mh = dgamma(nu_temp, shape = nu_gam_a, rate = nu_gam_b, log = T) +
              sum(dinvgamma(w[k,], shape = nu_temp*0.5, rate = nu_temp*0.5, log = T)) -
              log( pnorm(100, inits$nu[k], exp(logsigma_nu[k])) - pnorm(4, inits$nu[k], exp(logsigma_nu[k])) ) # truncation proposal
            denum_mh = dgamma(inits$nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
              sum(dinvgamma(w[k,], shape = inits$nu[k]*0.5, rate = inits$nu[k]*0.5, log = T))-
              log( pnorm(100, nu_temp,  exp(logsigma_nu[k])) - pnorm(4, nu_temp, exp(logsigma_nu[k])) )  # truncation proposal
            alpha = min(1, exp(num_mh - denum_mh))
            temp = runif(1)
            if (alpha > temp){
              acount_nu[k] = acount_nu[k] + 1
            }
            lpost[(j - burnin) %/% thin,k] <- log(alpha) # Take a negative for denominator

        } else {
          # Nominator
          k <- nu_id
          q_nu_nustar <- dnorm(inits$nu[k], mean = nu[k], sd = exp(logsigma_nu[k]) ) /
            ( pnorm(100, nu[k],  exp(logsigma_nu[k])) - pnorm(4, nu[k], exp(logsigma_nu[k])) ) # Move nu to inits$nu
          num_mh = dgamma(inits$nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = inits$nu[k]*0.5, rate = inits$nu[k]*0.5, log = T))
          denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))
          alpha = min(1, exp(num_mh - denum_mh))
          lpost[(j - burnin) %/% thin,k] <- log( alpha * q_nu_nustar)
        }

      }

      if ((j > inits$burnin) & (j %% inits$thin == 0))
        mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, gamma, nu, diag(sigma_h), h0, as.numeric(h), as.numeric(w))
      if (j %% 1000 == 0) {
        cat(" Iteration ", j, " ", round(acount_nu/j,2)," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2), " gamma ", round(gamma,2) ," \n")
        acount_w <- rep(0,t_max)
      }
      }
  }
  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("gamma",c(1:K), sep = ""),
                        paste("nu",c(1:K), sep = ""),
                        paste("sigma_h",c(1:K), sep = ""),
                        paste("lh0",c(1:K), sep = ""),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K))
  )

  return(lpost)
}

###########################################################################
#' @export
Chib.OT.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL,
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
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "OT", SV = TRUE)
  }
  # prior B
  b_prior = prior$b_prior
  V_b_prior = prior$V_b_prior
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior h0_mean
  h0_mean <- 2 * log(prior$sigma)
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
  B <- inits$B0
  b_sample <- as.numeric(B)
  a_sample <- as.numeric(A[lower.tri(A)])
  h <- inits$h
  sigma_h <- inits$sigma_h

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
  lpost <- rep(0, (samples - burnin)%/% thin)
  repeat_time <- 1
  if (cal_Nu) {
    lpost <- matrix(0,nrow = (samples - burnin)%/% thin, ncol = K)
    repeat_time <- K
  }

  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K*t_max + K*t_max,
                 ncol = (samples - burnin)%/% thin)

  for (nu_id in c(1:repeat_time)){
    # Multi degrees of freedom
    nu <- inits$nu
    logsigma_nu <- inits$logsigma_nu
    acount_nu <- rep(0,K)
    acount_w <- rep(0, t_max)
    # Init w as Gaussian
    w <- inits$w
    w_sqrt <- sqrt(w)
    w_sqrt_inv <- 1/w_sqrt

    for (j in c(1:samples)){

      # Sample B
      if(!fix_B) {
        wt <- as.vector(1/(exp(h/2)*w_sqrt))
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

        wt <- as.vector(1/(exp(h/2)*w_sqrt))
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

      # Sample vol
      ytilde <- (A %*% (yt - B %*% xt)) / w_sqrt
      aux <- sample_h_mod(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
      h <- aux$Sigtdraw
      sqrtvol <- aux$sigt
      h0 <- as.numeric(aux$h0)

      if(!fix_Sigma) {
        if (K>1) {
          sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
        } else {
          sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
        }

        sigma_h <- diag(mapply( GIGrvg::rgig, n = 1, lambda = - (t_max - 1)*0.5, chi = sse_2,
                                psi = 1/prior$sigma_S0 ), nrow = K)

      }

      if (cal_Sigma & (j > burnin) & (j %% thin == 0)){
        if (K>1) {
          sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
        } else {
          sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
        }


        lpost[(j - burnin) %/% thin] <- sum(mapply(GIGrvg::dgig, x = diag(inits$sigma_h), lambda = - (t_max - 1)*0.5,
                                                   chi = sse_2, psi = 1/prior$sigma_S0, log = T))
      }


      # Sample A0
      if(!fix_A) {
        u_std <- (yt - B %*%xt) # change from Gaussian
        u_neg <- - u_std
        a_sample <- rep(0, K * (K - 1) /2)
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- a_prior[id_start:id_end]
          V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
          a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = u_std[i,] / sqrtvol[i,] / w_sqrt[i,],
                                                       xsub = matrix(u_neg[1:(i-1),] / reprow(sqrtvol[i,],i-1) / reprow(w_sqrt[i,], i-1), nrow = i-1),
                                                       a_sub = a_sub,
                                                       V_a_sub = V_a_sub)
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

            ysub = u_std[i,] / sqrtvol[i,] / w_sqrt[i,]
            xsub = matrix(u_neg[1:(i-1),] / reprow(sqrtvol[i,],i-1) / reprow(w_sqrt[i,], i-1), nrow = i-1)

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
      w <- matrix( rinvgamma( n = K*t_max, shape = nu*0.5 + 1*0.5, rate = nu*0.5 + 0.5 * (u^2) / exp(h) ),
                   K, t_max )
      w_sqrt <- sqrt(w)


      # Sample nu nominator
      if(!fix_Nu) {
        # if not Chib calculation of nu, do it as normal
        if (!cal_Nu){
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
        } else {

          # change nu up to nu_id, fix all the rest
          nu_temp = nu + exp(logsigma_nu)*rnorm(K)
          for (k in c(1:nu_id)){
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

      }
      # Sample nu denominator
      if (cal_Nu & fix_Nu){
        # change nu up to nu_id, fix all the rest
        if (nu_id > 1){
          nu_temp = nu + exp(logsigma_nu)*rnorm(K)
          for (k in c(1:(nu_id-1) )){
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


      }


      # Now calculate the Chib estimation of p(nu[nu_id] | rest)
      if (cal_Nu & (j > burnin) & (j %% thin == 0)){
        # Denominator
        if (fix_Nu) {

          # nu up nu_id only
          k <- nu_id
          # Move from nu_star to nu, using truncated proposal
          while (TRUE){
            nu_temp = inits$nu[k] + exp(logsigma_nu[k])*rnorm(1)
            if (nu_temp > 4 && nu_temp < 100){
              break
            }
          }
            num_mh = dgamma(nu_temp, shape = nu_gam_a, rate = nu_gam_b, log = T) +
              sum(dinvgamma(w[k,], shape = nu_temp*0.5, rate = nu_temp*0.5, log = T)) -
              log( pnorm(100, inits$nu[k], exp(logsigma_nu[k])) - pnorm(4, inits$nu[k], exp(logsigma_nu[k])) ) # truncation proposal
            denum_mh = dgamma(inits$nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
              sum(dinvgamma(w[k,], shape = inits$nu[k]*0.5, rate = inits$nu[k]*0.5, log = T)) -
              log( pnorm(100, nu_temp,  exp(logsigma_nu[k])) - pnorm(4, nu_temp, exp(logsigma_nu[k])) )  # truncation proposal
            alpha = min(1, exp(num_mh - denum_mh))
            temp = runif(1)
            if (alpha > temp){
              acount_nu[k] = acount_nu[k] + 1
            }
            lpost[(j - burnin) %/% thin,k] <- log(alpha) # Take a negative for denominator

        } else {
          # Nominator
          k <- nu_id
          q_nu_nustar <- dnorm(inits$nu[k], mean = nu[k], sd = exp(logsigma_nu[k]) ) /
            ( pnorm(100, nu[k],  exp(logsigma_nu[k])) - pnorm(4, nu[k], exp(logsigma_nu[k])) ) # Move nu to inits$nu
          num_mh = dgamma(inits$nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = inits$nu[k]*0.5, rate = inits$nu[k]*0.5, log = T))
          denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))
          alpha = min(1, exp(num_mh - denum_mh))
          lpost[(j - burnin) %/% thin,k] <- log( alpha * q_nu_nustar)
        }

      }

      if ((j > inits$burnin) & (j %% inits$thin == 0))
        mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, nu, diag(sigma_h), h0, as.numeric(h), as.numeric(w))

      if (j %% 1000 == 0) {
        cat(" Iteration ", j, " ", round(acount_nu/j,2)," ", round(nu,2)," \n")
        acount_w <- rep(0,t_max)
      }
      }
  }
  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("nu",c(1:K), sep = ""),
                        paste("sigma_h",c(1:K), sep = ""),
                        paste("lh0",c(1:K), sep = ""),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)))

  return(lpost)
}

#############################################################################################
#' @export
Chib.OST.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL,
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
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "OST", SV = TRUE)
  }
  # prior B
  b_prior = prior$b_prior
  V_b_prior = prior$V_b_prior
  # prior sigma
  sigma0_T0 <- prior$sigma_T0
  sigma0_S0 <- prior$sigma_S0
  # prior h0_mean
  h0_mean <- 2 * log(prior$sigma)
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
  B <- inits$B0
  b_sample <- as.numeric(B)
  a_sample <- as.numeric(A[lower.tri(A)])
  h <- inits$h
  sigma_h <- inits$sigma_h

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

  # new precompute
  i <- K*m
  theta.prior.prec = matrix(0,i+K,i+K)
  theta.prior.prec[1:i,1:i] <- V_b_prior_inv
  theta.prior.prec[i+(1:K),i+(1:K)] <- solve( V_gamma_prior )
  theta.prior.precmean <- theta.prior.prec %*% c( b_prior, gamma_prior )

  # Output
  lpost <- rep(0, (samples - burnin)%/% thin)
  repeat_time <- 1
  if (cal_Nu) {
    lpost <- matrix(0,nrow = (samples - burnin)%/% thin, ncol = K)
    repeat_time <- K
  }

  mcmc <-  matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K + K*t_max + K*t_max,
                  ncol = (samples - burnin)%/% thin)
  if (cal_Nu) lpost <- matrix(0,nrow = (samples - burnin)%/% thin, ncol = K)
  for (nu_id in c(1:repeat_time)){
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
    w_sqrt_inv <- 1/w_sqrt

    for (j in c(1:samples)){

      # Sample B and gamma
      if(!fix_B) {
        wt <- as.vector( 1/sqrt(w*exp(h)) )
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
        wt <- as.vector( 1/sqrt(w*exp(h)) )
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
      # Sample vol
      ytilde <- (A %*% (yt - B %*% xt) - D %*% ( w - mu.xi )) / w_sqrt
      aux <- sample_h_mod(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
      h <- aux$Sigtdraw
      sqrtvol <- aux$sigt
      h0 <- as.numeric(aux$h0)

      if(!fix_Sigma) {
        if (K>1) {
          sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
        } else {
          sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
        }

        sigma_h <- diag(mapply( GIGrvg::rgig, n = 1, lambda = - (t_max - 1)*0.5, chi = sse_2,
                                psi = 1/prior$sigma_S0 ), nrow = K)

      }

      if (cal_Sigma & (j > burnin) & (j %% thin == 0)){
        if (K>1) {
          sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
        } else {
          sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
        }


        lpost[(j - burnin) %/% thin] <- sum(mapply(GIGrvg::dgig, x = diag(inits$sigma_h), lambda = - (t_max - 1)*0.5,
                                                   chi = sse_2, psi = 1/prior$sigma_S0, log = T))
      }

      # Sample A0
      if(!fix_A) {
        u_std <- (yt - B %*% xt) # change from Gaussian
        u_neg <- - u_std
        a_sample <- rep(0, K * (K - 1) /2)
        for (i in c(2:K)){
          id_end <- i*(i-1)/2
          id_start <- id_end - i + 2
          a_sub <- a_prior[id_start:id_end]
          V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
          a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = (u_std[i,] - (w[i,]  - mu.xi[i] ) * gamma[i]) / sqrtvol[i,] / w_sqrt[i,],
                                                       xsub = matrix(u_neg[1:(i-1),] / reprow(sqrtvol[i,],i-1) / reprow(w_sqrt[i,], i-1), nrow = i-1),
                                                       a_sub = a_sub,
                                                       V_a_sub = V_a_sub)
        }
        A_post <- matrix(0, nrow = K, ncol = K)
        A_post[upper.tri(A)] <- a_sample
        A <- t(A_post)
        diag(A) <- 1
      }

      if (cal_A & (j > burnin) & (j %% thin == 0)){
        u_std <- (yt - B %*% xt) # change from Gaussian
        u_neg <- - u_std
        A_star <- inits$A0[lower.tri(A)]
        if (K > 1) {
          for (i in c(2:K)){
            id_end <- i*(i-1)/2
            id_start <- id_end - i + 2
            a_sub <- prior$a_prior[id_start:id_end]
            V_a_sub <- prior$V_a_prior[id_start:id_end, id_start:id_end]
            sub_A_star <- A_star[id_start:id_end]

            ysub = (u_std[i,] - (w[i,]  - mu.xi[i] ) * gamma[i]) / sqrtvol[i,] / w_sqrt[i,]
            xsub = matrix(u_neg[1:(i-1),] / reprow(sqrtvol[i,],i-1) / reprow(w_sqrt[i,], i-1), nrow = i-1)

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
      w <- matrix( mapply( GIGrvg::rgig, n = 1, lambda = -(nu+1)*0.5, chi = nu + (u^2)/exp(h),
                           psi = gamma^2/exp(h) ), K, t_max )

      w_sqrt <- sqrt(w)

      # Sample nu nominator
      if(!fix_Nu) {
        # if not Chib calculation of nu, do it as normal
        if (!cal_Nu){
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

        } else {

          # change nu up to nu_id, fix all the rest
          nu_temp = nu + exp(logsigma_nu)*rnorm(K)
          for (k in c(1:nu_id)){
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

      }
      # Sample nu denominator
      if (cal_Nu & fix_Nu){
        # change nu up to nu_id, fix all the rest
        if (nu_id > 1){
          nu_temp = nu + exp(logsigma_nu)*rnorm(K)
          for (k in c(1:(nu_id-1) )){
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


      }


      # Now calculate the Chib estimation of p(nu[nu_id] | rest)
      if (cal_Nu & (j > burnin) & (j %% thin == 0)){
        # Denominator
        if (fix_Nu) {

          # nu up nu_id only
          k <- nu_id
          # Move from nu_star to nu, using truncated proposal
          while (TRUE){
            nu_temp = inits$nu[k] + exp(logsigma_nu[k])*rnorm(1)
            if (nu_temp > 4 && nu_temp < 100){
              break
            }
          }
            num_mh = dgamma(nu_temp, shape = nu_gam_a, rate = nu_gam_b, log = T) +
              sum(dinvgamma(w[k,], shape = nu_temp*0.5, rate = nu_temp*0.5, log = T)) -
              log( pnorm(100, inits$nu[k], exp(logsigma_nu[k])) - pnorm(4, inits$nu[k], exp(logsigma_nu[k])) ) # truncation proposal
            denum_mh = dgamma(inits$nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
              sum(dinvgamma(w[k,], shape = inits$nu[k]*0.5, rate = inits$nu[k]*0.5, log = T)) -
              log( pnorm(100, nu_temp,  exp(logsigma_nu[k])) - pnorm(4, nu_temp, exp(logsigma_nu[k])) )  # truncation proposal
            alpha = min(1, exp(num_mh - denum_mh))
            temp = runif(1)
            if (alpha > temp){
              acount_nu[k] = acount_nu[k] + 1
            }
            lpost[(j - burnin) %/% thin,k] <- log(alpha) # Take a negative for denominator

        } else {
          # Nominator
          k <- nu_id
          q_nu_nustar <- dnorm(inits$nu[k], mean = nu[k], sd = exp(logsigma_nu[k]) ) /
            ( pnorm(100, nu[k],  exp(logsigma_nu[k])) - pnorm(4, nu[k], exp(logsigma_nu[k])) ) # Move nu to inits$nu
          num_mh = dgamma(inits$nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = inits$nu[k]*0.5, rate = inits$nu[k]*0.5, log = T))
          denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
            sum(dinvgamma(w[k,], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))
          alpha = min(1, exp(num_mh - denum_mh))
          lpost[(j - burnin) %/% thin,k] <- log( alpha * q_nu_nustar)
        }

      }
      if ((j > inits$burnin) & (j %% inits$thin == 0))
        mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, gamma, nu, diag(sigma_h), h0, as.numeric(h), as.numeric(w))
      if (j %% 1000 == 0) {
        cat(" Iteration- ", j, " ", round(acount_nu/j,2)," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2) , " ", round(gamma,2) ,  " \n")
        acount_w <- rep(0,t_max)
      }
      }
  }
  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("gamma",c(1:K), sep = ""),
                        paste("nu",c(1:K), sep = ""),
                        paste("sigma_h",c(1:K), sep = ""),
                        paste("lh0",c(1:K), sep = ""),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)))

  return(lpost)
}

