#' Bayesian inference of VAR model with SV
#'
#' Bayesian inference of VAR model with SV
#' \deqn{y_t = B x_t + SQRT(w_t) A^(-1) Sigma eps_t}
#' @param y The input data as matrix T x K where T is the number of observations and K is the number of variables
#' @param K The number of variables in BVAR model.
#' @param p The number of lags in BVAR model.
#' @param dist The variable specifies the BVAR error distribution. It should be one of
#' c("Gaussian","Student","Skew.Student","Hyper.Student", "multiStudent","Skew.multiStudent","Hyper.multiStudent").
#' @param y0 The number of observations
#' @param prior The prior specification of BVAR.
#' @param inits The initial values of BVAR.
#' @return A coda object of the posterior samples.
#' @export
#' @examples
#' \dontrun{
#' datagen <- sim.VAR.SV(dist="Gaussian")
#' y <- datagen$y
#' prior <- get_prior(y, p = 2, dist="Gaussian", SV = T)
#' inits <- get_init(prior)
#' Chain1 <- BVAR.SV(y, K = 5, p = 2, dist = "Gaussian", y0 = NULL, prior = prior, inits = inits)
#' plot(Chain1)
#' }
BVAR.SV <- function(y, K, p, dist, y0 = NULL, prior = NULL, inits = NULL){

  if (dist == "Gaussian") Chain <- BVAR.Gaussian.SV(y, K, p, y0, prior, inits)
  if (dist == "Student") Chain <- BVAR.Student.SV(y, K, p, y0, prior, inits)
  if (dist == "Skew.Student") Chain <- BVAR.Skew.Student.SV(y, K, p, y0, prior, inits)
  if (dist == "Hyper.Student") Chain <- BVAR.Hyper.Student.SV(y, K, p, y0, prior, inits)
  if (dist == "multiStudent") Chain <- BVAR.multiStudent.SV(y, K, p, y0, prior, inits)
  if (dist == "Hyper.multiStudent") Chain <- BVAR.Hyper.multiStudent.SV(y, K, p, y0, prior, inits)
  return(Chain)
}
#' @export
BVAR.Gaussian.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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
  # prior A
  a_prior = prior$a_prior
  V_a_prior = prior$V_a_prior


  # Initial values
  if (is.null(inits)){
    inits <- get_init(prior)
  }
  samples <- inits$samples
  A <- inits$A0
  B <- inits$B0
  h <- inits$h
  sigma_h <- inits$sigma_h

  V_b_prior_inv <- solve(V_b_prior)

  # svdraw <- list()
  # paravol <- matrix(0, ncol = 3, nrow = K)
  # for (i in c(1:K)){
  #   svdraw[[i]] <- list(para = c(mu = 0, phi = 0.95, sigma = 0.2),
  #                       latent = h[i,])
  # }

  # Output
  mcmc <- NULL
  for (j in c(1:samples)){
    # Sample B, use reduce sum here
    b_post = rep(0, m*K)
    V_b_post_inv = solve(V_b_prior)
    for (i in c(1:t_max)){
      V_b_post_inv <- V_b_post_inv + kronecker(xt[,i] %*% t(xt[,i]), t(A)%*% diag(1/exp(h[,i])) %*% A)
      b_post <- b_post + kronecker(xt[,i], (t(A)%*% diag(1/exp(h[,i])) %*% A) %*% yt[,i])
    }

    V_b_post <- solve(V_b_post_inv)
    b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    B <- Vec_to_Mat(b_sample, K,p)

    # Sample vol
    ytilde <- A%*% (yt - B %*% xt)
    aux <- sample_h_ele( ytilde, sigma_h,  h, K, t_max)
    h <- aux$Sigtdraw
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h

    # for (i in c(1:K)){
    #   svdraw[[i]] <- svsample2(ytilde[i,], startpara = para(svdraw[[i]]),
    #                            startlatent = latent(svdraw[[i]]), priormu = priormu,
    #                            priorphi = priorphi, priorsigma = priorsigma)
    #   paravol[i,] <- para(svdraw[[i]])
    #   h[i,] <- as.vector(latent(svdraw[[i]]))
    # }
    # sqrtvol <- exp(h/2)

    # Sample A0
    u_std <- (yt - B %*%xt)
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
    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc <- cbind(mcmc, c(b_sample, a_sample, as.vector(sigma_h), as.vector(h)))
    if (j %% 100 == 0) { cat(" Iteration ", j, " \n")}
  }
  row.names(mcmc) <- c( paste("b",c(1:(m*K)), sep = ""),
                        paste("a",c(1:(K * (K - 1) /2)), sep = ""),
                        sprintf("sigmaH_%d_%d",rep.col(c(1:K),K), rep.row(c(1:K),K)),
                        sprintf("h_%d_%d", rep.col(c(1:K),t_max), rep.row(c(1:t_max),K))
  )
  return(as.mcmc(t(mcmc)))
}

###########################################################################
#' @export
BVAR.Student.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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
  samples <- inits$samples
  A <- inits$A0
  B <- inits$B0
  h <- inits$h
  sigma_h <- inits$sigma_h

  V_b_prior_inv <- solve(V_b_prior)

  nu <- inits$nu
  logsigma_nu <- 0
  acount_nu <- 0
  acount_w <- rep(0, t_max)
  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- rep.row(w_sample, K)
  w_sqrt <- sqrt(w)

  # svdraw <- list()
  # paravol <- matrix(0, ncol = 3, nrow = K)
  # for (i in c(1:K)){
  #   svdraw[[i]] <- list(para = c(mu = 0, phi = 0.95, sigma = 0.2),
  #                       latent = h[i,])
  # }

  # Output
  mcmc <- NULL
  for (j in c(1:samples)){
    # Sample B
    b_post = rep(0, m*K)
    V_b_post_inv = solve(V_b_prior)
    for (i in c(1:t_max)){
      V_b_post_inv <- V_b_post_inv + kronecker(xt[,i] %*% t(xt[,i]), t(A)%*% diag(1/exp(h[,i])) %*% A /w_sample[i])
      b_post <- b_post + kronecker(xt[,i], (t(A)%*% diag(1/exp(h[,i])) %*% A /w_sample[i]) %*% yt[,i])
    }

    V_b_post <- solve(V_b_post_inv)
    b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    B <- Vec_to_Mat(b_sample, K,p)

    # Sample vol
    ytilde <- A%*% (yt - B %*% xt)/w_sqrt
    aux <- sample_h_ele( ytilde, sigma_h,  h, K, t_max)
    h <- aux$Sigtdraw
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h

    # for (i in c(1:K)){
    #   svdraw[[i]] <- svsample2(ytilde[i,], startpara = para(svdraw[[i]]),
    #                            startlatent = latent(svdraw[[i]]), priormu = priormu,
    #                            priorphi = priorphi, priorsigma = priorsigma)
    #   paravol[i,] <- para(svdraw[[i]])
    #   h[i,] <- as.vector(latent(svdraw[[i]]))
    # }
    # sqrtvol <- exp(h/2)

    # Sample A0
    u_std <- (yt - B %*%xt)/ w_sqrt # change from Gaussian
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

    # Sample w
    u <- (yt - B %*%xt)
    for (i in c(1:t_max)){
      w_sample[i] <- rinvgamma(1, shape = nu*0.5 + K*0.5, rate = nu*0.5 + 0.5 * t(u[,i]) %*% (t(A)%*% diag(1/exp(h[,i])) %*% A) %*% u[,i])
    }
    w <- rep.row(w_sample, K)
    w_sqrt <- sqrt(w)


    # Sample nu
    nu_temp = nu + exp(logsigma_nu)*rnorm(1)
    if (nu_temp > 2 && nu_temp < 100){
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

    if(j %% batchlength == 0 ){
      if (acount_nu > batchlength * TARGACCEPT){
        logsigma_nu = logsigma_nu + adaptamount(j %/% batchlength);
      }
      if (acount_nu < batchlength * TARGACCEPT){
        logsigma_nu = logsigma_nu - adaptamount(j %/% batchlength);
      }
      acount_nu = 0
    }

    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc <- cbind(mcmc, c(b_sample, a_sample, nu, as.vector(sigma_h), as.vector(h), as.vector(w_sample)))
    if (j %% 100 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", nu," \n")
      acount_w <- rep(0,t_max)
    }
  }
  row.names(mcmc) <- c( paste("b",c(1:(m*K)), sep = ""),
                        paste("a",c(1:(K * (K - 1) /2)), sep = ""),
                        paste("nu"),
                        sprintf("sigma_h_%d_%d",rep.col(c(1:K),K), rep.row(c(1:K),K)),
                        sprintf("h_%d_%d", rep.col(c(1:K),t_max), rep.row(c(1:t_max),K)),
                        paste("w",c(1:t_max), sep = ""))

  return(as.mcmc(t(mcmc)))
}

#############################################################################################
# library(Matrix)
# library(magic)
# library(tmvnsim)
#' @export
BVAR.Skew.Student.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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
  samples <- inits$samples
  A <- inits$A0
  B <- inits$B0
  h <- inits$h
  sigma_h <- inits$sigma_h

  V_b_prior_inv <- solve(V_b_prior)

  nu <- inits$nu
  logsigma_nu <- 0
  acount_nu <- 0
  acount_w <- rep(0, t_max)
  gamma <- inits$gamma
  D <- diag(gamma)

  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- rep.row(w_sample, K)
  w_sqrt <- sqrt(w)
  # Init z as Truncated Gaussian
  z <- matrix(abs(rnorm(K * t_max)), ncol = t_max, nrow = K)

  # svdraw <- list()
  # paravol <- matrix(0, ncol = 3, nrow = K)
  # for (i in c(1:K)){
  #   svdraw[[i]] <- list(para = c(mu = 0, phi = 0.95, sigma = 0.2),
  #                       latent = h[i,])
  # }

  # Output
  mcmc <- NULL
  bg_prior <- c(gamma_prior, b_prior)
  #V_bg_prior <-  bdiag(V_gamma_prior, V_b_prior)  # Sparse matrix
  V_bg_prior <-  adiag(V_gamma_prior, V_b_prior)
  V_bg_prior_inv <- solve(V_bg_prior)

  for (j in c(1:samples)){
    # Sample B and gamma
    bg_post = rep(0,K + m*K) # first K elements are gamma
    V_bg_post_inv = V_bg_prior_inv
    for (i in c(1:t_max)){
      zx <- cbind(diag(z[,i]), kronecker(t(xt[,i]), diag(K)))
      V_bg_post_inv <- V_bg_post_inv +  t(zx) %*% ( (t(A)%*% diag(1/exp(h[,i])) %*% A) /w_sample[i]) %*% zx
      bg_post <- bg_post + t(zx) %*% ( (t(A)%*% diag(1/exp(h[,i])) %*% A) /w_sample[i]) %*% yt[,i]
    }
    V_bg_post <- solve(V_bg_post_inv)
    V_bg_post <- (V_bg_post + t(V_bg_post))/2
    bg_post <- V_bg_post %*% ( V_bg_prior_inv %*% bg_prior + bg_post)
    bg_sample <- bg_post + t(chol(V_bg_post)) %*% rnorm(K + m*K)
    gamma <- bg_sample[1:K]
    D <- diag(gamma)
    b_sample <- bg_sample[(K+1):(K+m*K)]
    B <- Vec_to_Mat(b_sample, K,p)

    # Sample Z
    lb <- rep(0, K)
    for (i in c(1:t_max)){
      V_Zt_inv <- diag(rep(1,K)) + D %*% ( (t(A)%*% diag(1/exp(h[,i])) %*% A) /w_sample[i]) %*% D
      V_Zt <- solve(V_Zt_inv)
      mu_Zt <- as.vector(V_Zt %*% D %*% ( (t(A)%*% diag(1/exp(h[,i])) %*% A) /w_sample[i]) %*% (yt[,i] - B %*% xt[,i]) )
      z[,i] <- TruncatedNormal::rtmvnorm(n = 1, mu = mu_Zt, sigma = V_Zt, lb = lb) # library(TruncatedNormal)
      # z[,i] <- tmvtnorm::rtmvnorm(1, mean = mu_Zt, sigma = V_Zt, lower=rep(0, length = K)) # library(tmvtnorm)
      # z[,i] <- tmvnsim(1, K, means = mu_Zt, sigma = V_Zt, lower=rep(0, length = K))$samp # library(tmvnsim)
    }

    # Sample vol
    ytilde <- A%*% (yt - B %*% xt - D%*% z)/w_sqrt    # change from Gaussian
    aux <- sample_h_ele( ytilde, sigma_h,  h, K, t_max)
    h <- aux$Sigtdraw
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h

    # for (i in c(1:K)){
    #   svdraw[[i]] <- svsample2(ytilde[i,], startpara = para(svdraw[[i]]),
    #                            startlatent = latent(svdraw[[i]]), priormu = priormu,
    #                            priorphi = priorphi, priorsigma = priorsigma)
    #   paravol[i,] <- para(svdraw[[i]])
    #   h[i,] <- as.vector(latent(svdraw[[i]]))
    # }
    # sqrtvol <- exp(h/2)

    # Sample A0
    u_std <- (yt - B %*%xt - D%*% z)/ w_sqrt # change from Gaussian
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

    # Sample w
    u <- (yt - B %*%xt - D%*% z)
    for (i in c(1:t_max)){
      w_sample[i] <- rinvgamma(1, shape = nu*0.5 + K*0.5, rate = nu*0.5 + 0.5 * t(u[,i]) %*% (t(A)%*% diag(1/exp(h[,i])) %*% A) %*% u[,i])
    }
    w <- rep.row(w_sample, K)
    w_sqrt <- sqrt(w)

    # Sample nu
    nu_temp = nu + exp(logsigma_nu)*rnorm(1)
    if (nu_temp > 2 && nu_temp < 100){
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

    if(j %% batchlength == 0 ){
      if (acount_nu > batchlength * TARGACCEPT){
        logsigma_nu = logsigma_nu + adaptamount(j %/% batchlength);
      }
      if (acount_nu < batchlength * TARGACCEPT){
        logsigma_nu = logsigma_nu - adaptamount(j %/% batchlength);
      }
      acount_nu = 0
    }
    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc <- cbind(mcmc, c(b_sample, a_sample, gamma, nu, as.vector(sigma_h), as.vector(h), as.vector(w_sample)))
    if (j %% 100 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", " ", nu ," \n")
      acount_w <- rep(0,t_max)
    }
  }
  row.names(mcmc) <- c( paste("b",c(1:(m*K)), sep = ""),
                        paste("a",c(1:(K * (K - 1) /2)), sep = ""),
                        paste("gamma",c(1:K), sep = ""),
                        paste("nu"),
                        sprintf("sigma_h_%d_%d",rep.col(c(1:K),K), rep.row(c(1:K),K)),
                        sprintf("h_%d_%d", rep.col(c(1:K),t_max), rep.row(c(1:t_max),K)),
                        paste("w",c(1:t_max), sep = ""))
  return(as.mcmc(t(mcmc)))
}

#############################################################################################
#' @export
BVAR.Hyper.Student.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)

  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Hyper.Student", SV = TRUE)
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
  samples <- inits$samples
  A <- inits$A0
  B <- inits$B0
  h <- inits$h
  sigma_h <- inits$sigma_h

  V_b_prior_inv <- solve(V_b_prior)

  nu <- inits$nu
  logsigma_nu <- 0
  acount_nu <- 0
  acount_w <- rep(0, t_max)
  gamma <- inits$gamma
  D <- diag(gamma)

  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- rep.row(w_sample, K)
  w_sqrt <- sqrt(w)

  # svdraw <- list()
  # paravol <- matrix(0, ncol = 3, nrow = K)
  # for (i in c(1:K)){
  #   svdraw[[i]] <- list(para = c(mu = 0, phi = 0.95, sigma = 0.2),
  #                       latent = h[i,])
  # }

  # Output
  mcmc <- NULL
  for (j in c(1:samples)){
    # Sample B
    b_post = rep(0, m*K)
    V_b_post_inv = solve(V_b_prior)
    for (i in c(1:t_max)){
      V_b_post_inv <- V_b_post_inv + kronecker(xt[,i] %*% t(xt[,i]), t(A)%*% diag(1/exp(h[,i])) %*% A /w_sample[i])
      b_post <- b_post + kronecker(xt[,i], (t(A)%*% diag(1/exp(h[,i])) %*% A /w_sample[i]) %*% (yt[,i] - gamma * w_sample[i] ) )
    }

    V_b_post <- solve(V_b_post_inv)
    b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    B <- Vec_to_Mat(b_sample, K,p)

    # Sample gamma
    gamma_post = rep(0, K)
    V_gamma_post_inv = solve(V_gamma_prior)
    for (i in c(1:t_max)){
      V_gamma_post_inv <- V_gamma_post_inv +  w_sample[i] * t(A)%*% diag(1/exp(h[,i])) %*% A # ww' sqrt(w_inv) Signma_inv sqrt(w_inv)
      gamma_post <- gamma_post + (t(A)%*% diag(1/exp(h[,i])) %*% A) %*% (yt[,i] - B %*% xt[,i] )
    }

    V_gamma_post <- solve(V_gamma_post_inv)
    gamma_post <- V_gamma_post %*% ( solve(V_gamma_prior) %*% gamma_prior + gamma_post)
    gamma <- as.vector(gamma_post + t(chol(V_gamma_post)) %*% rnorm(K))
    D <- diag(gamma)

    # Sample vol
    ytilde <- A%*% (yt - B %*% xt  - D%*% w)/w_sqrt
    aux <- sample_h_ele( ytilde, sigma_h,  h, K, t_max)
    h <- aux$Sigtdraw
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h

    # for (i in c(1:K)){
    #   svdraw[[i]] <- svsample2(ytilde[i,], startpara = para(svdraw[[i]]),
    #                            startlatent = latent(svdraw[[i]]), priormu = priormu,
    #                            priorphi = priorphi, priorsigma = priorsigma)
    #   paravol[i,] <- para(svdraw[[i]])
    #   h[i,] <- as.vector(latent(svdraw[[i]]))
    # }
    # sqrtvol <- exp(h/2)

    # Sample A0
    u_std <- (yt - B %*% xt - D%*% w)/ w_sqrt # change from Gaussian
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

    # Sample w
    u <- (yt - B %*%xt)
    for (i in c(1:t_max)){
      ran_stu = rt(1, df = 4)
      Sigma2_inv_t = (t(A)%*% diag(1/exp(h[,i])) %*% A)
      a_eq <- -0.5* (t(gamma) %*% Sigma2_inv_t %*% gamma)
      b_eq <- -((0.5*K + nu/2+1))
      c_eq =  0.5 * t(u[,i]) %*% Sigma2_inv_t %*% u[,i] + nu/2
      w_mode <- (- b_eq - sqrt(b_eq^2-4*a_eq*c_eq))/(2*a_eq)
      mode_target <- log( w_mode )
      cov_target <- 1/( c_eq / w_mode - a_eq * w_mode) * 1.2 # scale by 1.2

      log_w_temp = mode_target + sqrt(cov_target)*ran_stu
      w_temp = exp(log_w_temp)

      num_mh =  dinvgamma(w_temp, shape = nu*0.5, rate = nu*0.5, log = T) +  # prior
        ( - 0.5 * K * log_w_temp - 0.5 * t((u[,i] - gamma*w_temp)) %*% Sigma2_inv_t %*% (u[,i] - gamma*w_temp) / w_temp) - # posterior
        dt( (log_w_temp - mode_target)/sqrt(cov_target), df = 4, log = T) + log_w_temp # proposal
      denum_mh = dinvgamma(w_sample[i], shape = nu*0.5, rate = nu*0.5, log = T) +  # prior
        ( - 0.5 * K * log(w_sample[i]) - 0.5 * t((u[,i] - gamma*w_sample[i])) %*% Sigma2_inv_t %*% (u[,i] - gamma*w_sample[i]) / w_sample[i]) - # posterior
        dt( (log(w_sample[i]) - mode_target)/sqrt(cov_target), df = 4, log = T) + log(w_sample[i]) # proposal
      alpha = num_mh - denum_mh;
      temp = log(runif(1))
      if (alpha > temp){
        w_sample[i] <- w_temp
        w[,i] <- w_temp
        w_sqrt[,i] <- sqrt(w_temp)
        acount_w[i] <- acount_w[i] + 1
      }
    }

    # Sample nu
    nu_temp = nu + exp(logsigma_nu)*rnorm(1)
    if (nu_temp > 2 && nu_temp < 100){
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

    if(j %% batchlength == 0 ){
      if (acount_nu > batchlength * TARGACCEPT){
        logsigma_nu = logsigma_nu + adaptamount(j %/% batchlength);
      }
      if (acount_nu < batchlength * TARGACCEPT){
        logsigma_nu = logsigma_nu - adaptamount(j %/% batchlength);
      }
      acount_nu = 0
    }
    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc <- cbind(mcmc, c(b_sample, a_sample, gamma, nu, as.vector(sigma_h), as.vector(h), as.vector(w_sample)))
    if (j %% 100 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", nu ,  " \n")
      acount_w <- rep(0,t_max)
    }
  }
  row.names(mcmc) <- c( paste("b",c(1:(m*K)), sep = ""),
                        paste("a",c(1:(K * (K - 1) /2)), sep = ""),
                        paste("gamma",c(1:K), sep = ""),
                        paste("nu"),
                        sprintf("sigma_h_%d_%d",rep.col(c(1:K),K), rep.row(c(1:K),K)),
                        sprintf("h_%d_%d", rep.col(c(1:K),t_max), rep.row(c(1:t_max),K)),
                        paste("w",c(1:t_max), sep = ""))

  return(as.mcmc(t(mcmc)))
}

###########################################################################
#' @export
BVAR.multiStudent.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)

  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "multiStudent", SV = TRUE)
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
  samples <- inits$samples
  A <- inits$A0
  B <- inits$B0
  h <- inits$h
  sigma_h <- inits$sigma_h

  V_b_prior_inv <- solve(V_b_prior)

  # Multi degrees of freedom
  nu <- inits$nu
  logsigma_nu <- rep(0,K)
  acount_nu <- rep(0,K)
  acount_w <- rep(0, t_max)
  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- rep.row(w_sample, K)
  w_sqrt <- sqrt(w)
  w_sqrt_inv <- 1/w_sqrt

  # svdraw <- list()
  # paravol <- matrix(0, ncol = 3, nrow = K)
  # for (i in c(1:K)){
  #   svdraw[[i]] <- list(para = c(mu = 0, phi = 0.95, sigma = 0.2),
  #                       latent = h[i,])
  # }

  # Output
  mcmc <- NULL
  for (j in c(1:samples)){
    # Sample B
    b_post = rep(0, m*K)
    V_b_post_inv = solve(V_b_prior)
    for (i in c(1:t_max)){
      V_b_post_inv <- V_b_post_inv + kronecker(xt[,i] %*% t(xt[,i]), diag(w_sqrt_inv[,i]) %*% t(A)%*% diag(1/exp(h[,i])) %*% A %*% diag(w_sqrt_inv[,i]))
      b_post <- b_post + kronecker(xt[,i], (diag(w_sqrt_inv[,i]) %*% t(A)%*% diag(1/exp(h[,i])) %*% A %*% diag(w_sqrt_inv[,i])) %*% yt[,i])
    }

    V_b_post <- solve(V_b_post_inv)
    b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    B <- Vec_to_Mat(b_sample, K,p)

    # Sample vol
    ytilde <- A%*% ((yt - B %*%xt)/ w_sqrt)
    aux <- sample_h_ele( ytilde, sigma_h,  h, K, t_max)
    h <- aux$Sigtdraw
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h

    # for (i in c(1:K)){
    #   svdraw[[i]] <- svsample2(ytilde[i,], startpara = para(svdraw[[i]]),
    #                            startlatent = latent(svdraw[[i]]), priormu = priormu,
    #                            priorphi = priorphi, priorsigma = priorsigma)
    #   paravol[i,] <- para(svdraw[[i]])
    #   h[i,] <- as.vector(latent(svdraw[[i]]))
    # }
    # sqrtvol <- exp(h/2)

    # Sample A0
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

    # Sample w
    u <- (yt - B %*%xt)
    u_proposal <- A %*% (yt - B %*%xt)
    for (i in c(1:t_max)){
      a_target <- (nu*0.5 + 1*0.5) * 0.75 # adjust by 0.75
      b_target <- (nu*0.5 + 0.5 * u_proposal[,i]^2 / sqrtvol[,i]^2) * 0.75 # adjust by 0.75
      w_temp = w[,i]
      for (k in c(1:K)){
        w_temp[k] <- rinvgamma(1, shape = a_target[k], rate = b_target[k])
      }
      log_w_temp <- log(w_temp)
      w_temp_sqrt_inv = 1 / sqrt(w_temp)

      #total_log_dens <- (-0.5) * sum(log(w[,i])) - 0.5 * t(u[,i]) %*% diag(w_sqrt_inv[,i]) %*% Sigma2_inv %*% diag(w_sqrt_inv[,i]) %*% u[,i]  + sum(dinvgamma(w[,i], shape = nu*0.5, rate = nu*0.5, log = T))
      Sigma2_inv_t = t(A)%*% diag(1/exp(h[,i])) %*% A
      num_mh =  sum(dinvgamma(w_temp, shape = nu*0.5, rate = nu*0.5, log = T)) +  # prior
        ( - 0.5 * sum(log_w_temp) - 0.5 * t(u[,i]) %*% diag(w_temp_sqrt_inv) %*% Sigma2_inv_t %*% diag(w_temp_sqrt_inv) %*% u[,i]) - # posterior
        sum( dinvgamma(w_temp, shape = a_target, rate = b_target, log = T)) # proposal
      denum_mh = sum(dinvgamma(w[,i], shape = nu*0.5, rate = nu*0.5, log = T)) + # prior
        (-0.5) * sum(log(w[,i])) - 0.5 * t(u[,i]) %*% diag(w_sqrt_inv[,i]) %*% Sigma2_inv_t %*% diag(w_sqrt_inv[,i]) %*% u[,i]  - # posterior
        sum( dinvgamma(w[,i], shape = a_target, rate = b_target, log = T))
      alpha = num_mh - denum_mh;
      temp = log(runif(1))
      if (alpha > temp){
        w[,i] <- w_temp
        w_sqrt[,i] <- sqrt(w_temp)
        w_sqrt_inv[,i] <- 1/sqrt(w_temp)
        acount_w[i] <- acount_w[i] + 1
      }
    }



    # Sample nu
    nu_temp = nu + exp(logsigma_nu)*rnorm(K)
    for (k in c(1:K)){
      if (nu_temp[k] > 2 && nu_temp[k] < 100){
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

    if(j %% batchlength == 0 ){
      for (jj in c(1:K)) {
        if (acount_nu[jj] > batchlength * TARGACCEPT){
          logsigma_nu[jj] = logsigma_nu[jj] + adaptamount(j %/% batchlength);
        }
        if (acount_nu[jj] < batchlength * TARGACCEPT){
          logsigma_nu[jj] = logsigma_nu[jj] - adaptamount(j %/% batchlength);
        }
        acount_nu[jj] = 0
      }
    }

    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc <- cbind(mcmc, c(b_sample, a_sample, nu, as.vector(sigma_h), as.vector(h), as.vector(w)))

    if (j %% 100 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", nu ,  " \n")
      acount_w <- rep(0,t_max)
    }
  }
  row.names(mcmc) <- c( paste("b",c(1:(m*K)), sep = ""),
                        paste("a",c(1:(K * (K - 1) /2)), sep = ""),
                        paste("nu",c(1:K), sep = ""),
                        sprintf("sigma_h_%d_%d",rep.col(c(1:K),K), rep.row(c(1:K),K)),
                        sprintf("h_%d_%d", rep.col(c(1:K),t_max), rep.row(c(1:t_max),K)),
                        sprintf("w_%d_%d", rep.col(c(1:K),t_max), rep.row(c(1:t_max),K)))

  return(as.mcmc(t(mcmc)))
}

#############################################################################################
#' @export
BVAR.Hyper.multiStudent.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)

  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "Hyper.multiStudent", SV = TRUE)
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
  samples <- inits$samples
  A <- inits$A0
  B <- inits$B0
  h <- inits$h
  sigma_h <- inits$sigma_h

  V_b_prior_inv <- solve(V_b_prior)

  # Multi degrees of freedom
  nu <- inits$nu
  logsigma_nu <- rep(0,K)
  acount_nu <- rep(0,K)
  acount_w <- rep(0, t_max)
  gamma <- inits$gamma
  D <- diag(gamma)

  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- rep.row(w_sample, K)
  w_sqrt <- sqrt(w)
  w_sqrt_inv <- 1/w_sqrt

  # svdraw <- list()
  # paravol <- matrix(0, ncol = 3, nrow = K)
  # for (i in c(1:K)){
  #   svdraw[[i]] <- list(para = c(mu = 0, phi = 0.95, sigma = 0.2),
  #                       latent = h[i,])
  # }

  # Output
  mcmc <- NULL
  for (j in c(1:samples)){
    # Sample B
    b_post = rep(0, m*K)
    V_b_post_inv = solve(V_b_prior)
    for (i in c(1:t_max)){
      V_b_post_inv <- V_b_post_inv + kronecker(xt[,i] %*% t(xt[,i]), diag(w_sqrt_inv[,i]) %*% (t(A)%*% diag(1/exp(h[,i])) %*% A) %*% diag(w_sqrt_inv[,i]) )
      b_post <- b_post + kronecker(xt[,i], (diag(w_sqrt_inv[,i]) %*% (t(A)%*% diag(1/exp(h[,i])) %*% A) %*% diag(w_sqrt_inv[,i])) %*% (yt[,i] - gamma * w[,i]))
    }

    V_b_post <- solve(V_b_post_inv)
    b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    B <- Vec_to_Mat(b_sample, K,p)

    gamma_post = rep(0, K)
    V_gamma_post_inv = solve(V_gamma_prior)
    for (i in c(1:t_max)){
      V_gamma_post_inv <- V_gamma_post_inv +  diag(w_sqrt[,i]) %*% (t(A)%*% diag(1/exp(h[,i])) %*% A) %*% diag(w_sqrt[,i]) # ww' sqrt(w_inv) Signma_inv sqrt(w_inv)
      gamma_post <- gamma_post + (diag(w_sqrt[,i]) %*% (t(A)%*% diag(1/exp(h[,i])) %*% A) %*% diag(w_sqrt_inv[,i])) %*% (yt[,i] - B %*% xt[,i] )
    }

    V_gamma_post <- solve(V_gamma_post_inv)
    gamma_post <- V_gamma_post %*% ( solve(V_gamma_prior) %*% gamma_prior + gamma_post)
    gamma <- as.vector(gamma_post + t(chol(V_gamma_post)) %*% rnorm(K))
    D <- diag(gamma)

    # Sample vol
    ytilde <- A%*% ((yt - B %*%xt - D %*% w)/ w_sqrt)
    aux <- sample_h_ele( ytilde, sigma_h,  h, K, t_max)
    h <- aux$Sigtdraw
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h

    # for (i in c(1:K)){
    #   svdraw[[i]] <- svsample2(ytilde[i,], startpara = para(svdraw[[i]]),
    #                            startlatent = latent(svdraw[[i]]), priormu = priormu,
    #                            priorphi = priorphi, priorsigma = priorsigma)
    #   paravol[i,] <- para(svdraw[[i]])
    #   h[i,] <- as.vector(latent(svdraw[[i]]))
    # }
    # sqrtvol <- exp(h/2)

    # Sample A0
    u_std <- (yt - B %*% xt - D%*% w) / w_sqrt # change from Gaussian
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

    # Sample w
    u <- (yt - B %*%xt)
    u_proposal <- A %*% (yt - B %*%xt)
    for (i in c(1:t_max)){
      a_target <- (nu*0.5 + 1*0.5) * 0.75
      b_target <- (nu*0.5 + 0.5 * u_proposal[,i]^2 / sqrtvol[,i]^2)*0.75  # adjust by 0.75
      w_temp = w[,i]
      for (k in c(1:K)){
        w_temp[k] <- rinvgamma(1, shape = a_target[k], rate = b_target[k])
      }
      log_w_temp <- log(w_temp)
      w_temp_sqrt_inv = 1 / sqrt(w_temp)

      #total_log_dens <- (-0.5) * sum(log(w[,i])) - 0.5 * t(u[,i]) %*% diag(w_sqrt_inv[,i]) %*% Sigma2_inv %*% diag(w_sqrt_inv[,i]) %*% u[,i]  + sum(dinvgamma(w[,i], shape = nu*0.5, rate = nu*0.5, log = T))
      Sigma2_inv_t = t(A)%*% diag(1/exp(h[,i])) %*% A
      num_mh =  sum(dinvgamma(w_temp, shape = nu*0.5, rate = nu*0.5, log = T)) +  # prior
        ( - 0.5 * sum(log_w_temp) - 0.5 * t(u[,i] - gamma*w_temp) %*% diag(w_temp_sqrt_inv) %*% Sigma2_inv_t %*% diag(w_temp_sqrt_inv) %*% (u[,i] - gamma*w_temp)) - # posterior
        sum( dinvgamma(w_temp, shape = a_target, rate = b_target, log = T)) # proposal
      denum_mh = sum(dinvgamma(w[,i], shape = nu*0.5, rate = nu*0.5, log = T)) + # prior
        (-0.5) * sum(log(w[,i])) - 0.5 * t(u[,i] - gamma*w[,i]) %*% diag(w_sqrt_inv[,i]) %*% Sigma2_inv_t %*% diag(w_sqrt_inv[,i]) %*% (u[,i] - gamma*w[,i])  - # posterior
        sum( dinvgamma(w[,i], shape = a_target, rate = b_target, log = T))
      alpha = num_mh - denum_mh;
      temp = log(runif(1))
      if (alpha > temp){
        w[,i] <- w_temp
        w_sqrt[,i] <- sqrt(w_temp)
        w_sqrt_inv[,i] <- 1/sqrt(w_temp)
        acount_w[i] <- acount_w[i] + 1
      }
    }

    # Sample nu
    nu_temp = nu + exp(logsigma_nu)*rnorm(K)
    for (k in c(1:K)){
      if (nu_temp[k] > 2 && nu_temp[k] < 100){
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

    if(j %% batchlength == 0 ){
      for (jj in c(1:K)) {
        if (acount_nu[jj] > batchlength * TARGACCEPT){
          logsigma_nu[jj] = logsigma_nu[jj] + adaptamount(j %/% batchlength);
        }
        if (acount_nu[jj] < batchlength * TARGACCEPT){
          logsigma_nu[jj] = logsigma_nu[jj] - adaptamount(j %/% batchlength);
        }
        acount_nu[jj] = 0
      }
    }
    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc <- cbind(mcmc, c(b_sample, a_sample, gamma, nu, sigma_h,  as.vector(h), as.vector(w)))
    if (j %% 100 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", nu ,  " \n")
      acount_w <- rep(0,t_max)
    }
  }
  row.names(mcmc) <- c( paste("b",c(1:(m*K)), sep = ""),
                        paste("a",c(1:(K * (K - 1) /2)), sep = ""),
                        paste("gamma",c(1:K), sep = ""),
                        paste("nu",c(1:K), sep = ""),
                        sprintf("sigma_h_%d_%d",rep.col(c(1:K),K), rep.row(c(1:K),K)),
                        sprintf("h_%d_%d", rep.col(c(1:K),t_max), rep.row(c(1:t_max),K)),
                        sprintf("w_%d_%d", rep.col(c(1:K),t_max), rep.row(c(1:t_max),K))
  )

  return(as.mcmc(t(mcmc)))
}
