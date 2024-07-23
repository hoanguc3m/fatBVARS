#' Bayesian inference of VAR model with RW-SV
#'
#' Bayesian inference of VAR model with RW-SV
#' \deqn{y_t = B x_t + SQRT(w_t) A^(-1) H_t eps_t}
#' @param y The input data as matrix T x K where T is the number of observations and K is the number of variables
#' @param K The number of variables in BVAR model.
#' @param p The number of lags in BVAR model.
#' @param dist The variable specifies the BVAR error distribution. It should be one of
#' c("Gaussian","Student","Skew.Student","Skew.Student", "MT","Skew.MT","MST").
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
  if (!(dist %in% c("Gaussian","Student","Skew.Student",
                    "MT","MST",
                    "OT","OST",
                    "dynSkew.Student", "dynMST", "dynOST") ))
    stop("dist is not implemented.")
  if (prior$SV == TRUE){
    Start = Sys.time()
    if (dist == "Gaussian") Chain <- BVAR.Gaussian.SV(y, K, p, y0, prior, inits)
    if (dist == "Student") Chain <- BVAR.Student.SV(y, K, p, y0, prior, inits)
    if (dist == "Skew.Student") Chain <- BVAR.Skew.Student.SV(y, K, p, y0, prior, inits)
    if (dist == "MT") Chain <- BVAR.MT.SV(y, K, p, y0, prior, inits)
    if (dist == "MST") Chain <- BVAR.MST.SV(y, K, p, y0, prior, inits)
    if (dist == "OT") Chain <- BVAR.OT.SV(y, K, p, y0, prior, inits)
    if (dist == "OST") Chain <- BVAR.OST.SV(y, K, p, y0, prior, inits)

    if (dist == "dynSkew.Student") Chain <- BVAR.dynSkew.Student.SV(y, K, p, y0, prior, inits)
    if (dist == "dynMST") Chain <- BVAR.dynMST.SV(y, K, p, y0, prior, inits)
    if (dist == "dynOST") Chain <- BVAR.dynOST.SV(y, K, p, y0, prior, inits)
    elapsedTime = Sys.time() - Start
    print(elapsedTime)
    out <- list(mcmc = Chain,
                y = y,
                y0 = y0,
                K = K,
                p = p,
                dist = dist,
                prior = prior,
                inits = inits,
                esttime = elapsedTime)
    class(out) <- c("fatBVARSV")
    return(out)
  } else {
    warning("prior$SV is TRUE")
  }
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
  # prior h0_mean
  h0_mean <- 2 * log(prior$sigma)
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
  # new precompute
  theta.prior.precmean <- V_b_prior_inv %*% b_prior

  # svdraw <- list()
  # paravol <- matrix(0, ncol = 3, nrow = K)
  # for (i in c(1:K)){
  #   svdraw[[i]] <- list(para = c(mu = 0, phi = 0.95, sigma = 0.2),
  #                       latent = h[i,])
  # }

  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K*t_max,
                     ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # Sample B, use reduce sum here
    # b_post = rep(0, m*K)
    # V_b_post_inv = V_b_prior_inv
    # for (i in c(1:t_max)){
    #   V_b_post_inv <- V_b_post_inv + kronecker(xt[,i] %*% t(xt[,i]), t(A)%*% diag(1/exp(h[,i])) %*% A)
    #   b_post <- b_post + kronecker(xt[,i], (t(A)%*% diag(1/exp(h[,i])) %*% A) %*% yt[,i])
    # }

    # V_b_post_inv <- V_b_prior_inv +
    #   Reduce( f = "+",
    #           x = lapply(1:t_max, function(i) kronecker(xt[,i] %*% t(xt[,i]), t(A)%*% diag(1/exp(h[,i]), nrow = K) %*% A)),
    #           accumulate = FALSE)
    # b_post <- Reduce(f = "+",
    #                  x = lapply(1:t_max, function(i) kronecker(xt[,i], (t(A)%*% diag(1/exp(h[,i]), nrow = K) %*% A) %*% yt[,i])),
    #                  accumulate = FALSE)
    #
    # V_b_post <- solve(V_b_post_inv)
    # b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    # b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    # B <- Vec_to_Mat(b_sample, K,p)

    # Sample B
    wt <- as.vector(exp(-h/2))
    y.tilde <- as.vector( A %*% yt ) * wt
    x.tilde <- kronecker( t(xt), A ) * wt
    theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
    b_sample <- backsolve( theta.prec.chol,
                        backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                   upper.tri = T, transpose = T )
                        + rnorm(K*m) )
    B <- matrix(b_sample,K,m)

    # Sample vol
    ytilde <- A%*% (yt - B %*% xt)
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h


    # for (i in c(1:K)){
    #   svdraw[[i]] <- svsample2(ytilde[i,], startpara = para(svdraw[[i]]),
    #                            startlatent = latent(svdraw[[i]]), priormu = priormu,
    #                            priorphi = priorphi, priorsigma = priorsigma)
    #   paravol[i,] <- para(svdraw[[i]])
    #   h[i,] <- as.numeric(latent(svdraw[[i]]))
    # }
    # sqrtvol <- exp(h/2)

    # Sample A0
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
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)

  # new precompute
  theta.prior.precmean <- V_b_prior_inv %*% b_prior

  # svdraw <- list()
  # paravol <- matrix(0, ncol = 3, nrow = K)
  # for (i in c(1:K)){
  #   svdraw[[i]] <- list(para = c(mu = 0, phi = 0.95, sigma = 0.2),
  #                       latent = h[i,])
  # }

  # Output
  mcmc <-  matrix(NA, nrow = m*K + 0.5*K*(K-1) + 1 + K + K + K*t_max + t_max,
                  ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # # Sample B
    # b_post = rep(0, m*K)
    # V_b_post_inv = V_b_prior_inv
    # for (i in c(1:t_max)){
    #   V_b_post_inv <- V_b_post_inv + kronecker(xt[,i] %*% t(xt[,i]), t(A)%*% diag(1/exp(h[,i]), nrow = K) %*% A /w_sample[i])
    #   b_post <- b_post + kronecker(xt[,i], (t(A)%*% diag(1/exp(h[,i]), nrow = K) %*% A /w_sample[i]) %*% yt[,i])
    # }
    #
    # V_b_post <- solve(V_b_post_inv)
    # b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    # b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    # B <- Vec_to_Mat(b_sample, K,p)

    # Sample B
    wt <- as.vector(1/(exp(h/2)*w_sqrt))
    y.tilde <- as.vector( A %*% yt ) * wt
    x.tilde <- kronecker( t(xt), A ) * wt
    theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
    b_sample <- backsolve( theta.prec.chol,
                        backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                   upper.tri = T, transpose = T )
                        + rnorm(K*m) )
    B <- matrix(b_sample,K,m)

    # Sample vol
    ytilde <- A%*% ((yt - B %*% xt)/w_sqrt)
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h

    # for (i in c(1:K)){
    #   svdraw[[i]] <- svsample2(ytilde[i,], startpara = para(svdraw[[i]]),
    #                            startlatent = latent(svdraw[[i]]), priormu = priormu,
    #                            priorphi = priorphi, priorsigma = priorsigma)
    #   paravol[i,] <- para(svdraw[[i]])
    #   h[i,] <- as.numeric(latent(svdraw[[i]]))
    #   h0[i,] <- as.numeric(latent0(svdraw[[i]]))
    # }
    # sqrtvol <- exp(h/2)

    # Sample A0
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

    # # Sample w
    # u <- (yt - B %*%xt)
    # for (i in c(1:t_max)){
    #   w_sample[i] <- rinvgamma(1, shape = nu*0.5 + K*0.5, rate = nu*0.5 + 0.5 * t(u[,i]) %*% (t(A)%*% diag(1/exp(h[,i]), nrow = K) %*% A) %*% u[,i])
    # }
    # w <- reprow(w_sample, K)
    # w_sqrt <- sqrt(w)

    # Sample w
    u <-  A %*% (yt - B %*%xt)
    w_sample <- rinvgamma( n = t_max, shape = (nu+K)*0.5, rate = 0.5*( nu + colSums((u^2)/exp(h)) ) )
    w <- reprow(w_sample, K)
    w_sqrt <- sqrt(w)

    # Sample nu
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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, nu, diag(sigma_h), h0, as.numeric(h), as.numeric(w_sample))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", round(nu,2)," \n")
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

  return(as.mcmc(t(mcmc)))
}

#############################################################################################
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
  samples <- inits$samples
  A <- inits$A0
  B <- inits$B0
  h <- inits$h
  sigma_h <- inits$sigma_h

  V_b_prior_inv <- solve(V_b_prior)

  nu <- inits$nu
  mu.xi <- nu/( nu - 2 )
  logsigma_nu <- 0
  acount_nu <- 0
  acount_w <- rep(0, t_max)
  gamma <- inits$gamma
  D <- diag(gamma, nrow = K)

  # Init w as Gaussian
  w_sample <- rep(1, t_max)
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
                  ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # # Sample B
    # b_post = rep(0, m*K)
    # V_b_post_inv = V_b_prior_inv
    # for (i in c(1:t_max)){
    #   V_b_post_inv <- V_b_post_inv + kronecker(xt[,i] %*% t(xt[,i]), t(A)%*% diag(1/exp(h[,i]), nrow = K) %*% A /w_sample[i])
    #   b_post <- b_post + kronecker(xt[,i], (t(A)%*% diag(1/exp(h[,i]), nrow = K) %*% A /w_sample[i]) %*% (yt[,i] - gamma * w_sample[i] ) )
    # }
    #
    # V_b_post <- solve(V_b_post_inv)
    # b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    # b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    # B <- Vec_to_Mat(b_sample, K,p)
    #
    # # Sample gamma
    # gamma_post = rep(0, K)
    # V_gamma_post_inv = solve(V_gamma_prior)
    # for (i in c(1:t_max)){
    #   V_gamma_post_inv <- V_gamma_post_inv +  w_sample[i] * t(A)%*% diag(1/exp(h[,i]), nrow = K) %*% A # ww' sqrt(w_inv) Signma_inv sqrt(w_inv)
    #   gamma_post <- gamma_post + (t(A)%*% diag(1/exp(h[,i]), nrow = K) %*% A) %*% (yt[,i] - B %*% xt[,i] )
    # }
    #
    # V_gamma_post <- solve(V_gamma_post_inv)
    # gamma_post <- V_gamma_post %*% ( solve(V_gamma_prior) %*% gamma_prior + gamma_post)
    # gamma <- as.numeric(gamma_post + t(chol(V_gamma_post)) %*% rnorm(K))
    # D <- diag(gamma, nrow = K)

    # Sample B and gamma
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
    D <- diag(gamma, K)

    # Sample vol
    ytilde <- A%*% ((yt - B %*% xt  - D%*% (w-mu.xi) )/w_sqrt)
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h

    # Sample A0
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


    # Sample w
    q2 <- colSums( ( ( A %*% ( yt - B%*%xt + mu.xi * gamma ) ) / sqrtvol )^2 )
    p2 <- colSums( ( c( A %*% gamma ) / sqrtvol )^2 )
    w_sample <- mapply( GIGrvg::rgig, n = 1, lambda = -(nu+K)*0.5, chi = nu + q2,
                         psi = p2 )
    w <- reprow(w_sample, K)
    w_sqrt <- sqrt(w)

    # Sample nu
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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, gamma, nu, diag(sigma_h), h0, as.numeric(h), as.numeric(w_sample))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2) , " ", gamma ,  " \n")
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

  return(as.mcmc(t(mcmc)))
}

###########################################################################
#' @export
BVAR.MT.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)
#  w_sqrt_inv <- 1/w_sqrt

  # svdraw <- list()
  # paravol <- matrix(0, ncol = 3, nrow = K)
  # for (i in c(1:K)){
  #   svdraw[[i]] <- list(para = c(mu = 0, phi = 0.95, sigma = 0.2),
  #                       latent = h[i,])
  # }

  # new precompute
  theta.prior.precmean <- V_b_prior_inv %*% b_prior

  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K*t_max + K*t_max,
                ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # # Sample B
    # b_post = rep(0, m*K)
    # V_b_post_inv = solve(V_b_prior)
    # for (i in c(1:t_max)){
    #   V_b_post_inv <- V_b_post_inv + kronecker(xt[,i] %*% t(xt[,i]), diag(w_sqrt_inv[,i]) %*% t(A)%*% diag(1/exp(h[,i])) %*% A %*% diag(w_sqrt_inv[,i]))
    #   b_post <- b_post + kronecker(xt[,i], (diag(w_sqrt_inv[,i]) %*% t(A)%*% diag(1/exp(h[,i])) %*% A %*% diag(w_sqrt_inv[,i])) %*% yt[,i])
    # }
    #
    # V_b_post <- solve(V_b_post_inv)
    # b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    # b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    # B <- Vec_to_Mat(b_sample, K,p)

    # Sample B
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

    # Sample vol
    ytilde <- A%*% ((yt - B %*%xt)/ w_sqrt)
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h

    # for (i in c(1:K)){
    #   svdraw[[i]] <- svsample2(ytilde[i,], startpara = para(svdraw[[i]]),
    #                            startlatent = latent(svdraw[[i]]), priormu = priormu,
    #                            priorphi = priorphi, priorsigma = priorsigma)
    #   paravol[i,] <- para(svdraw[[i]])
    #   h[i,] <- as.numeric(latent(svdraw[[i]]))
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
    a_target <- (nu*0.5 + 1*0.5) * 0.75
    b_target <- (nu*0.5 + 0.5 * u_proposal^2 / sqrtvol^2)*0.75  # adjust by 0.75
    w_temp <- matrix( rinvgamma( n = K*t_max, shape = a_target, rate = b_target ), K, t_max ) # rinvgamma recycles arguments
    w_temp_sqrt <- sqrt(w_temp)
    # num_mh <- sapply(1:t_max, function (i) { Sigma2_inv_t = t(A)%*% diag(1/exp(h[,i])) %*% A; u.w <- w_temp_sqrt_inv[,i] * u[,i]
    #                                          sum(dinvgamma(w_temp[,i], shape = nu*0.5, rate = nu*0.5, log = T)) + #prior
    #                                           ( - 0.5 * sum(log_w_temp[,i]) - 0.5 * t(u.w) %*% Sigma2_inv_t %*% u.w ) - # posterior
    #                                           sum( dinvgamma(w_temp[,i], shape = a_target, rate = b_target[,i], log = T)) # proposal
    #                                        }
    #                   )
    # denum_mh <- sapply(1:t_max, function (i) { Sigma2_inv_t = t(A)%*% diag(1/exp(h[,i])) %*% A; u.w <- u[,i]/sqrt(w[,i])
    #                                            sum(dinvgamma(w[,i], shape = nu*0.5, rate = nu*0.5, log = T)) + #prior
    #                                            ( - 0.5 * sum(log(w[,i])) - 0.5 * t(u.w) %*% Sigma2_inv_t %*% u.w ) - # posterior
    #                                            sum( dinvgamma(w[,i], shape = a_target, rate = b_target[,i], log = T)) # proposal
    #                                          }
    #                   )
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

    # # Sample w
    # u <- (yt - B %*%xt)
    # u_proposal <- A %*% (yt - B %*%xt)
    # for (i in c(1:t_max)){
    #   a_target <- (nu*0.5 + 1*0.5) * 0.75 # adjust by 0.75
    #   b_target <- (nu*0.5 + 0.5 * u_proposal[,i]^2 / sqrtvol[,i]^2) * 0.75 # adjust by 0.75
    #   w_temp = w[,i]
    #   for (k in c(1:K)){
    #     w_temp[k] <- rinvgamma(1, shape = a_target[k], rate = b_target[k])
    #   }
    #   log_w_temp <- log(w_temp)
    #   w_temp_sqrt_inv = 1 / sqrt(w_temp)
    #
    #   #total_log_dens <- (-0.5) * sum(log(w[,i])) - 0.5 * t(u[,i]) %*% diag(w_sqrt_inv[,i]) %*% Sigma2_inv %*% diag(w_sqrt_inv[,i]) %*% u[,i]  + sum(dinvgamma(w[,i], shape = nu*0.5, rate = nu*0.5, log = T))
    #   Sigma2_inv_t = t(A)%*% diag(1/exp(h[,i])) %*% A
    #   num_mh =  sum(dinvgamma(w_temp, shape = nu*0.5, rate = nu*0.5, log = T)) +  # prior
    #     ( - 0.5 * sum(log_w_temp) - 0.5 * t(u[,i]) %*% diag(w_temp_sqrt_inv) %*% Sigma2_inv_t %*% diag(w_temp_sqrt_inv) %*% u[,i]) - # posterior
    #     sum( dinvgamma(w_temp, shape = a_target, rate = b_target, log = T)) # proposal
    #   denum_mh = sum(dinvgamma(w[,i], shape = nu*0.5, rate = nu*0.5, log = T)) + # prior
    #     (-0.5) * sum(log(w[,i])) - 0.5 * t(u[,i]) %*% diag(w_sqrt_inv[,i]) %*% Sigma2_inv_t %*% diag(w_sqrt_inv[,i]) %*% u[,i]  - # posterior
    #     sum( dinvgamma(w[,i], shape = a_target, rate = b_target, log = T))
    #   alpha = num_mh - denum_mh;
    #   temp = log(runif(1))
    #   if (alpha > temp){
    #     w[,i] <- w_temp
    #     w_sqrt[,i] <- sqrt(w_temp)
    #     w_sqrt_inv[,i] <- 1/sqrt(w_temp)
    #     acount_w[i] <- acount_w[i] + 1
    #   }
    # }



    # Sample nu
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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, nu, diag(sigma_h), h0, as.numeric(h), as.numeric(w))

    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2), " \n")
      acount_w <- rep(0,t_max)
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

  return(as.mcmc(t(mcmc)))
}

#############################################################################################
#' @export
BVAR.MST.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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
  samples <- inits$samples
  A <- inits$A0
  B <- inits$B0
  h <- inits$h
  sigma_h <- inits$sigma_h

  V_b_prior_inv <- solve(V_b_prior)

  # Multi degrees of freedom
  nu <- inits$nu
  mu.xi <- nu/( nu - 2 )
  logsigma_nu <- rep(0,K)
  acount_nu <- rep(0,K)
  acount_w <- rep(0, t_max)
  gamma <- inits$gamma
  D <- diag(gamma)

  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)

  # new precompute
  i <- K*m
  theta.prior.prec = matrix(0,i+K,i+K)
  theta.prior.prec[1:i,1:i] <- V_b_prior_inv
  theta.prior.prec[i+(1:K),i+(1:K)] <- solve( V_gamma_prior )
  theta.prior.precmean <- theta.prior.prec %*% c( b_prior, gamma_prior )

  # Output
  mcmc <-  matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K + K*t_max + K*t_max,
                  ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # # Sample B
    # b_post = rep(0, m*K)
    # V_b_post_inv = V_b_prior_inv
    # for (i in c(1:t_max)){
    #   V_b_post_inv <- V_b_post_inv + kronecker(xt[,i] %*% t(xt[,i]), diag(w_sqrt_inv[,i]) %*% (t(A)%*% diag(1/exp(h[,i])) %*% A) %*% diag(w_sqrt_inv[,i]) )
    #   b_post <- b_post + kronecker(xt[,i], (diag(w_sqrt_inv[,i]) %*% (t(A)%*% diag(1/exp(h[,i])) %*% A) %*% diag(w_sqrt_inv[,i])) %*% (yt[,i] - gamma * w[,i]))
    # }
    #
    # V_b_post <- solve(V_b_post_inv)
    # b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    # b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    # B <- Vec_to_Mat(b_sample, K,p)
    #
    # gamma_post = rep(0, K)
    # V_gamma_post_inv = solve(V_gamma_prior)
    # for (i in c(1:t_max)){
    #   V_gamma_post_inv <- V_gamma_post_inv +  diag(w_sqrt[,i]) %*% (t(A)%*% diag(1/exp(h[,i])) %*% A) %*% diag(w_sqrt[,i]) # ww' sqrt(w_inv) Signma_inv sqrt(w_inv)
    #   gamma_post <- gamma_post + (diag(w_sqrt[,i]) %*% (t(A)%*% diag(1/exp(h[,i])) %*% A) %*% diag(w_sqrt_inv[,i])) %*% (yt[,i] - B %*% xt[,i] )
    # }
    #
    # V_gamma_post <- solve(V_gamma_post_inv)
    # gamma_post <- V_gamma_post %*% ( solve(V_gamma_prior) %*% gamma_prior + gamma_post)
    # gamma <- as.numeric(gamma_post + t(chol(V_gamma_post)) %*% rnorm(K))
    # D <- diag(gamma)


    # Sample B and gamma
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

    # Sample vol
    ytilde <- A%*% ((yt - B %*%xt - D %*% ( w - mu.xi ))/ w_sqrt)
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h

    # Sample A0
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

    # Sample w
    u <- (yt - B %*%xt + mu.xi * gamma  )
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

    # Sample nu
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
    mu.xi <- nu / ( nu - 2 )

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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, gamma, nu, diag(sigma_h), h0, as.numeric(h), as.numeric(w))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2), " gamma ", round(gamma,2) ," \n")
      acount_w <- rep(0,t_max)
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

  return(as.mcmc(t(mcmc)))
}

###########################################################################
#' @export
BVAR.OT.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)

  # svdraw <- list()
  # paravol <- matrix(0, ncol = 3, nrow = K)
  # for (i in c(1:K)){
  #   svdraw[[i]] <- list(para = c(mu = 0, phi = 0.95, sigma = 0.2),
  #                       latent = h[i,])
  # }

  # new precompute
  theta.prior.precmean <- V_b_prior_inv %*% b_prior

  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K*t_max + K*t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # # Sample B
    # b_post = rep(0, m*K)
    # V_b_post_inv = V_b_prior_inv
    # for (i in c(1:t_max)){
    #   V_b_post_inv <- V_b_post_inv + kronecker(xt[,i] %*% t(xt[,i]), t(A)%*% diag(1/exp(h[,i]) / w[,i] ) %*% A )
    #   b_post <- b_post + kronecker(xt[,i], (t(A)%*% diag(1/exp(h[,i])/ w[,i] ) %*% A) %*% yt[,i])
    # }
    #
    # V_b_post <- solve(V_b_post_inv)
    # b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    # b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    # B <- Vec_to_Mat(b_sample, K,p)

    # Sample B
    wt <- as.vector(1/(exp(h/2)*w_sqrt))
    y.tilde <- as.vector( A %*% yt ) * wt
    x.tilde <- kronecker( t(xt), A ) * wt
    theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
    b_sample <- backsolve( theta.prec.chol,
                        backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                   upper.tri = T, transpose = T )
                            + rnorm(K*m) )
    B <- matrix(b_sample,K,m)

    # Sample vol
    ytilde <- (A %*% (yt - B %*% xt)) / w_sqrt
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h

    # for (i in c(1:K)){
    #   svdraw[[i]] <- svsample2(ytilde[i,], startpara = para(svdraw[[i]]),
    #                            startlatent = latent(svdraw[[i]]), priormu = priormu,
    #                            priorphi = priorphi, priorsigma = priorsigma)
    #   paravol[i,] <- para(svdraw[[i]])
    #   h[i,] <- as.numeric(latent(svdraw[[i]]))
    # }
    # sqrtvol <- exp(h/2)

    # Sample A0
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

    # # Sample w
    # u <-  A %*% (yt - B %*%xt)
    # for (i in c(1:t_max)){
    #   w[,i] <- mapply(rinvgamma, n = 1, shape = nu*0.5 + 1*0.5, rate = nu*0.5 + 0.5 * (u[,i]^2) / exp(h[,i]))
    # }
    # w_sqrt <- sqrt(w)

    # Sample w
    u <-  A %*% (yt - B %*%xt)
    w <- matrix( rinvgamma( n = K*t_max, shape = nu*0.5 + 1*0.5, rate = nu*0.5 + 0.5 * (u^2) / exp(h) ),
                 K, t_max )
    w_sqrt <- sqrt(w)


    # Sample nu
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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, nu, diag(sigma_h), h0, as.numeric(h), as.numeric(w))

    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", round(nu,2)," \n")
      acount_w <- rep(0,t_max)
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

  return(as.mcmc(t(mcmc)))
}

#############################################################################################
#' @export
BVAR.OST.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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
  samples <- inits$samples
  A <- inits$A0
  B <- inits$B0
  h <- inits$h
  sigma_h <- inits$sigma_h

  V_b_prior_inv <- solve(V_b_prior)

  # Multi degrees of freedom
  nu <- inits$nu
  mu.xi <- nu/( nu - 2 )
  logsigma_nu <- rep(0,K)
  acount_nu <- rep(0,K)
  acount_w <- rep(0, t_max)
  gamma <- inits$gamma
  D <- diag(gamma)

  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)

  # new precompute
  i <- K*m
  theta.prior.prec = matrix(0,i+K,i+K)
  theta.prior.prec[1:i,1:i] <- V_b_prior_inv
  theta.prior.prec[i+(1:K),i+(1:K)] <- solve( V_gamma_prior )
  theta.prior.precmean <- theta.prior.prec %*% c( b_prior, gamma_prior )

  # Output
  mcmc <-  matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K + K*t_max + K*t_max,
                  ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # # Sample B
    # b_post = rep(0, m*K)
    # V_b_post_inv = V_b_prior_inv
    # inv_A <- solve(A)
    # for (i in c(1:t_max)){
    #   V_b_post_inv <- V_b_post_inv + kronecker(xt[,i] %*% t(xt[,i]), t(A)%*% diag(1/exp(h[,i]) / w[,i] ) %*% A )
    #   b_post <- b_post + kronecker(xt[,i], (t(A)%*% diag(1/exp(h[,i])/ w[,i] ) %*% A) %*% (yt[,i] - inv_A %*% (gamma * w[,i])) )
    # }
    #
    # V_b_post <- solve(V_b_post_inv)
    # b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    # b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    # B <- Vec_to_Mat(b_sample, K,p)
    #
    # # Sample gamma
    # gamma_post = rep(0, K)
    # V_gamma_post_inv = solve(V_gamma_prior)
    # u <- A %*% (yt - B %*%xt)
    #
    # V_gamma_post_inv <- solve(V_gamma_prior) + diag(apply(w / exp(h), MARGIN = 1, FUN = sum))
    # gamma_post <- apply(u / exp(h), MARGIN = 1, FUN = sum)
    #
    # V_gamma_post <- solve(V_gamma_post_inv)
    # gamma_post <- V_gamma_post %*% ( solve(V_gamma_prior) %*% gamma_prior + gamma_post)
    # gamma <- as.numeric(gamma_post + t(chol(V_gamma_post)) %*% rnorm(K))
    # D <- diag(gamma)

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
    theta <- backsolve( theta.prec.chol,
                        backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                   upper.tri = T, transpose = T )
                        + rnorm(K*(m+1)) )
    B <- matrix(theta[1:(m*K)],K,m)
    b_sample <- as.vector(B)
    gamma <- theta[(m*K+1):((m+1)*K)]
    D <- diag(gamma)

    # Sample vol
    ytilde <- (A %*% (yt - B %*% xt) - D %*% ( w - mu.xi ) ) / w_sqrt
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h

    # Sample A0
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

    # Sample w
    u <-  A %*% (yt - B %*%xt) + gamma * mu.xi
    w <- matrix( mapply( GIGrvg::rgig, n = 1, lambda = -(nu+1)*0.5, chi = nu + (u^2)/exp(h),
                         psi = gamma^2/exp(h) ), K, t_max )
    w_sqrt <- sqrt(w)

    # Sample nu
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
    mu.xi <- nu / ( nu - 2 )

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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, gamma, nu, diag(sigma_h), h0, as.numeric(h), as.numeric(w))
    if (j %% 1000 == 0) {
      cat(" Iteration- ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2) , " ", round(gamma,2) ,  " \n")
      acount_w <- rep(0,t_max)
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

  return(as.mcmc(t(mcmc)))
}


#############################################################################################
#' @export
BVAR.dynMST.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)

  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "dynMST", SV = TRUE)
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
  Gamma_T <- inits$gamma
  if (length(Gamma_T) <= K){
    Gamma_T = repcol(Gamma_T, t_max)
  }
  gamma0 = Gamma_T[,1]
  sigma_G = diag(1e-4, nrow = K)
  #D <- diag(gamma, nrow = K)


  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)
  w_sqrt_inv <- 1/w_sqrt

  # svdraw <- list()
  # paravol <- matrix(0, ncol = 3, nrow = K)
  # for (i in c(1:K)){
  #   svdraw[[i]] <- list(para = c(mu = 0, phi = 0.95, sigma = 0.2),
  #                       latent = h[i,])
  # }

  # Output
  mcmc <-  matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K + K + K * t_max + K*t_max + K*t_max,
                  ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # Sample B
    b_post = rep(0, m*K)
    V_b_post_inv = V_b_prior_inv
    for (i in c(1:t_max)){
      V_b_post_inv <- V_b_post_inv + kronecker(xt[,i] %*% t(xt[,i]), diag(w_sqrt_inv[,i]) %*% (t(A)%*% diag(1/exp(h[,i])) %*% A) %*% diag(w_sqrt_inv[,i]) )
      b_post <- b_post + kronecker(xt[,i], (diag(w_sqrt_inv[,i]) %*% (t(A)%*% diag(1/exp(h[,i])) %*% A) %*% diag(w_sqrt_inv[,i])) %*% (yt[,i] - Gamma_T[,i] * w[,i]))
    }

    V_b_post <- solve(V_b_post_inv)
    b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    B <- Vec_to_Mat(b_sample, K,p)

    # Sample gamma : carterkohn(y,Z,Ht)
    # y: n*T
    # Z: (nT)*(n_max*K)
    # Ht: (nT)*n
    # ysim <-  Z %*% B_t + t(chol(Hsim)) %*% rnorm(M)
    Sigma_T <- matrix(0, nrow = K * t_max, ncol = K)
    Z_T <- matrix(0, nrow = K * t_max, ncol = K)
    A_inv <- solve(A)
    for (i in c(1:t_max)){
      Sigma_T[ ((i-1)*K + 1):(i*K) , ] <- diag(w_sqrt[,i], nrow = K) %*% (A_inv %*% diag(exp(h[,i]), nrow = K) %*% t(A_inv)) %*% diag(w_sqrt[,i], nrow = K)
      Z_T[ ((i-1)*K + 1):(i*K) , ] <- diag(w[,i], nrow = K)
    }
    Gamma_filter <- fatBVARS::carterkohn(yt - B %*% xt, Z_T,
                                         Sigma_T, sigma_G, K, K, t_max,
                                         gamma_prior,V_gamma_prior)

    Gamma_T <- Gamma_filter$bdraws
    gamma0 <- Gamma_filter$b0draws

    # Draw sigma_G, the covariance of B(t)
    if (K>1) {
      sse_2 <- apply( (Gamma_T[,1:t_max] - cbind(gamma0,Gamma_T[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
    } else {
      sse_2 <- sum( (Gamma_T[,1:t_max] - c(gamma0,Gamma_T[,1:(t_max-1)]) )^2)
    }
    sigma_G_a <- rep(t_max,K) # prior of sigma_G
    sigma_G_b <- sse_2 # prior of sigma_G

    for (i in c(1:K)){
      sigma_new <- rinvgamma(1, shape = sigma_G_a[i] * 0.5, rate = sigma_G_b[i] * 0.5)
      alpha = (sigma_G[i,i] - sigma_new) / 2 + 0.5 * (log(sigma_new) - log(sigma_G[i,i])) # B_sigma = 1
      temp = log(runif(1))
      if (alpha > temp){
        sigma_G[i,i] <- sigma_new
      }
    }


    # Sample vol
    ytilde <- A%*% ((yt - B %*%xt - Gamma_T * w)/ w_sqrt)
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h

    # for (i in c(1:K)){
    #   svdraw[[i]] <- svsample2(ytilde[i,], startpara = para(svdraw[[i]]),
    #                            startlatent = latent(svdraw[[i]]), priormu = priormu,
    #                            priorphi = priorphi, priorsigma = priorsigma)
    #   paravol[i,] <- para(svdraw[[i]])
    #   h[i,] <- as.numeric(latent(svdraw[[i]]))
    # }
    # sqrtvol <- exp(h/2)

    # Sample A0
    u_std <- (yt - B %*% xt - Gamma_T * w) / w_sqrt # change from Gaussian
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
        ( - 0.5 * sum(log_w_temp) - 0.5 * t(u[,i] - Gamma_T[,i]*w_temp) %*% diag(w_temp_sqrt_inv) %*% Sigma2_inv_t %*% diag(w_temp_sqrt_inv) %*% (u[,i] - Gamma_T[,i]*w_temp)) - # posterior
        sum( dinvgamma(w_temp, shape = a_target, rate = b_target, log = T)) # proposal
      denum_mh = sum(dinvgamma(w[,i], shape = nu*0.5, rate = nu*0.5, log = T)) + # prior
        (-0.5) * sum(log(w[,i])) - 0.5 * t(u[,i] - Gamma_T[,i]*w[,i]) %*% diag(w_sqrt_inv[,i]) %*% Sigma2_inv_t %*% diag(w_sqrt_inv[,i]) %*% (u[,i] - Gamma_T[,i]*w[,i])  - # posterior
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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, nu, diag(sigma_h), h0, diag(sigma_G), gamma0, as.numeric(Gamma_T), as.numeric(h), as.numeric(w))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2),
          round(apply(Gamma_T,MARGIN = 1,FUN = mean),2) ,  " \n")
      acount_w <- rep(0,t_max)
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
                        paste("sigma_G",c(1:K), sep = ""),
                        paste("G0",c(1:K), sep = ""),
                        sprintf("gamma_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K))
  )

  return(as.mcmc(t(mcmc)))
}

#############################################################################################
#' @export
BVAR.dynOST.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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
  Gamma_T <- inits$gamma
  if (length(Gamma_T) <= K){
    Gamma_T = repcol(Gamma_T, t_max)
  }
  gamma0 = Gamma_T[,1]
  sigma_G = diag(1e-4, nrow = K)
  #D <- diag(gamma, nrow = K)

  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)

  # svdraw <- list()
  # paravol <- matrix(0, ncol = 3, nrow = K)
  # for (i in c(1:K)){
  #   svdraw[[i]] <- list(para = c(mu = 0, phi = 0.95, sigma = 0.2),
  #                       latent = h[i,])
  # }

  # Output
  mcmc <-  matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K + K + K * t_max + K*t_max + K*t_max,
                  ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # Sample B
    b_post = rep(0, m*K)
    V_b_post_inv = V_b_prior_inv
    inv_A <- solve(A)
    for (i in c(1:t_max)){
      V_b_post_inv <- V_b_post_inv + kronecker(xt[,i] %*% t(xt[,i]), t(A)%*% diag(1/exp(h[,i]) / w[,i] ) %*% A )
      b_post <- b_post + kronecker(xt[,i], (t(A)%*% diag(1/exp(h[,i])/ w[,i] ) %*% A)  ) %*% (yt[,i] - inv_A %*% (Gamma_T[,i] * w[,i]))
    }

    V_b_post <- solve(V_b_post_inv)
    b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    B <- Vec_to_Mat(b_sample, K,p)

    # Sample gamma : carterkohn(y,Z,Ht)
    # y: n*T
    # Z: (nT)*(n_max*K)
    # Ht: (nT)*n
    # ysim <-  Z %*% B_t + t(chol(Hsim)) %*% rnorm(M)
    u <- A %*% (yt - B %*%xt)
    Sigma_T <- matrix(0, nrow = K * t_max, ncol = K)
    Z_T <- matrix(0, nrow = K * t_max, ncol = K)
    A_inv <- solve(A)
    for (i in c(1:t_max)){
      Sigma_T[ ((i-1)*K + 1):(i*K) , ] <- diag(w[,i] * exp(h[,i]), nrow = K)
      Z_T[ ((i-1)*K + 1):(i*K) , ] <- diag(w[,i], nrow = K)
    }
    Gamma_filter <- fatBVARS::carterkohn(u, Z_T,
                                         Sigma_T, sigma_G, K, K, t_max,
                                         gamma_prior,V_gamma_prior)

    Gamma_T <- Gamma_filter$bdraws
    gamma0 <- Gamma_filter$b0draws

    # Draw sigma_G, the covariance of B(t)
    if (K>1) {
      sse_2 <- apply( (Gamma_T[,1:t_max] - cbind(gamma0,Gamma_T[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
    } else {
      sse_2 <- sum( (Gamma_T[,1:t_max] - c(gamma0,Gamma_T[,1:(t_max-1)]) )^2)
    }
    sigma_G_a <- rep(t_max,K) # prior of sigma_G
    sigma_G_b <- sse_2 # prior of sigma_G

    for (i in c(1:K)){
      sigma_new <- rinvgamma(1, shape = sigma_G_a[i] * 0.5, rate = sigma_G_b[i] * 0.5)
      alpha = (sigma_G[i,i] - sigma_new) / 2 + 0.5 * (log(sigma_new) - log(sigma_G[i,i])) # B_sigma = 1
      temp = log(runif(1))
      if (alpha > temp){
        sigma_G[i,i] <- sigma_new
      }
    }

    # Sample vol
    ytilde <- (A %*% (yt - B %*% xt) - Gamma_T * w) / w_sqrt
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h

    # Sample A0
    u_std <- (yt - B %*% xt) # change from Gaussian
    u_neg <- - u_std
    a_sample <- rep(0, K * (K - 1) /2)
    for (i in c(2:K)){
      id_end <- i*(i-1)/2
      id_start <- id_end - i + 2
      a_sub <- a_prior[id_start:id_end]
      V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
      a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = (u_std[i,] - w[i,] * Gamma_T[i,]) / sqrtvol[i,] / w_sqrt[i,],
                                                   xsub = matrix(u_neg[1:(i-1),] / reprow(sqrtvol[i,],i-1) / reprow(w_sqrt[i,], i-1), nrow = i-1),
                                                   a_sub = a_sub,
                                                   V_a_sub = V_a_sub)
    }
    A_post <- matrix(0, nrow = K, ncol = K)
    A_post[upper.tri(A)] <- a_sample
    A <- t(A_post)
    diag(A) <- 1

    # Sample w
    u <- A %*% (yt - B %*%xt)
    for (i in c(1:t_max)){
      ran_stu = rt(K, df = 4)
      Sigma2_inv_t = 1/exp(h[,i])
      a_eq <- -0.5* (Gamma_T[,i]^2 * Sigma2_inv_t)
      b_eq <- -((0.5*1 + nu/2+1))
      c_eq =  0.5 * (u[,i]^2 * Sigma2_inv_t) + nu/2
      w_mode <- (- b_eq - sqrt(b_eq^2-4*a_eq*c_eq))/(2*a_eq)
      w_mode <- ifelse(abs(w_mode)<1e-5, - c_eq/b_eq, w_mode)
      mode_target <- log( w_mode )
      cov_target <- 1/( c_eq / w_mode - a_eq * w_mode) * 1.2 # scale by 1.2

      log_w_temp = mode_target + sqrt(cov_target)*ran_stu
      w_temp = exp(log_w_temp)

      num_mh =  mapply( dinvgamma, x = w_temp, shape = nu*0.5, rate = nu*0.5, log = T) +  # prior
        ( - 0.5 * 1 * log_w_temp - 0.5 * (u[,i] - Gamma_T[,i]*w_temp)^2 * Sigma2_inv_t / w_temp) - # posterior
        dt( (log_w_temp - mode_target)/sqrt(cov_target), df = 4, log = T) + log_w_temp # proposal
      denum_mh = mapply( dinvgamma, x = w[,i], shape = nu*0.5, rate = nu*0.5, log = T) +  # prior
        ( - 0.5 * 1 * log(w[,i]) - 0.5 * (u[,i] - Gamma_T[,i]*w[,i])^2 * Sigma2_inv_t / w[,i]) - # posterior
        dt( (log(w[,i]) - mode_target)/sqrt(cov_target), df = 4, log = T) + log(w[,i]) # proposal
      alpha = num_mh - denum_mh
      temp = log(runif(K))
      acre = alpha > temp
      w[,i] <- ifelse(acre, w_temp, w[,i])
      w_sqrt[,i] <- sqrt(w[,i])
      acount_w[i] <- acount_w[i] + mean(acre)
    }

    # Sample nu
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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, nu, diag(sigma_h), h0, diag(sigma_G), gamma0, as.numeric(Gamma_T), as.numeric(h), as.numeric(w))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2) , " ",
          round(apply(Gamma_T,MARGIN = 1,FUN = mean),2) ,  " \n")
      acount_w <- rep(0,t_max)
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
                        paste("sigma_G",c(1:K), sep = ""),
                        paste("G0",c(1:K), sep = ""),
                        sprintf("gamma_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)))

  return(as.mcmc(t(mcmc)))
}

