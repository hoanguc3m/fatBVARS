#############################################################################################
#' @export
BVAR.dynSkew.Student.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)

  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "dynSkew.Student", SV = FALSE)
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
  B <- Vec_to_Mat(inits$B0, K, p)
  sigma <- inits$sigma
  # Sigma <- solve(inits$A0) %*% diag(inits$sigma, nrow = K)
  # Sigma2 <- Sigma %*% t(Sigma)
  # Sigma2_inv <- solve(Sigma2)
  V_b_prior_inv <- solve(V_b_prior)

  nu <- inits$nu
  mu.xi <- nu/( nu - 2 )
  logsigma_nu <- 0
  acount_nu <- 0
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

  # new precompute
  i <- K*m
  theta.prior.prec <- V_b_prior_inv
  theta.prior.precmean <- theta.prior.prec %*% b_prior

  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + 1 + K + K * t_max + t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # Sample B
    #yt_G <- (yt - Gamma_T * w)  / w_sqrt
    # Sample B and gamma
    wt <- as.vector(1/(sigma*w_sqrt))
    y.tilde <- as.vector( A %*% (yt - Gamma_T * (w- mu.xi) ) ) * wt
    x.tilde <- kronecker( t(xt), A ) * wt
    theta.prec.chol <- chol( theta.prior.prec + crossprod(x.tilde) )
    theta <- backsolve( theta.prec.chol,
                        backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                   upper.tri = T, transpose = T )
                        + rnorm(K*m) )
    B <- matrix(theta,K,m)
    b_sample <- as.vector(B)
    # gamma <- theta[(m*K+1):((m+1)*K)]
    # D <- diag(gamma, K)

    # Sample gamma : carterkohn(y,Z,Ht)
    # y: n*T
    # Z: (nT)*(n_max*K)
    # Ht: (nT)*n
    # ysim <-  Z %*% B_t + t(chol(Hsim)) %*% rnorm(M)
    Sigma_T <- matrix(0, nrow = K * t_max, ncol = K)
    A_inv <- solve(A)
    for (i in c(1:t_max)){
      Sigma_T[ ((i-1)*K + 1):(i*K) , ] <- A_inv %*% diag(sigma^2, nrow = K) %*% t(A_inv) * w_sample[i]
    }
    Gamma_filter <- fatBVARS::carterkohn(yt - B %*% xt, (w_sample - mu.xi) %x% diag(K),
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


    # Sample sigma
    sigma2 <- rep(0,K)
    u_ort <- A %*% ((yt - B %*% xt - Gamma_T * (w- mu.xi) )/ w_sqrt)    # change from Gaussian
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
    #    sigma_post_b <- sigma0_S0 + apply(u_ort^2, 1, sum)
    sigma_post_b <- sigma0_S0 + rowSums(u_ort^2)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)

    # Sample A0
    u_std <- (yt - B %*% xt - Gamma_T * (w- mu.xi))/ w_sqrt # change from Gaussian
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

    Sigma <- solve(A) %*% diag(sigma, nrow = length(sigma))
    Sigma2 <- Sigma %*% t(Sigma)
    Sigma2_inv <- solve(Sigma2)
    Sigma2_inv <- (Sigma2_inv + t(Sigma2_inv))*0.5

    # Sample w
    q2 <- colSums( ( ( A %*% ( yt - B%*%xt + mu.xi * Gamma_T) ) / sigma )^2 )
    p2 <- sum( ( c( A %*% Gamma_T ) / sigma )^2 )
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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, nu, diag(sigma_G), gamma0, as.numeric(Gamma_T), as.numeric(w_sample))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2), " ", round(apply(Gamma_T,MARGIN = 1,FUN = mean),2), " \n")
      acount_w <- rep(0,t_max)
    }
  }
  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("sigma",c(1:K), sep = ""),
                        paste("nu"),
                        paste("sigma_G",c(1:K), sep = ""),
                        paste("G0",c(1:K), sep = ""),
                        sprintf("gamma_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        paste("w",c(1:t_max), sep = ""))

  return(as.mcmc(t(mcmc)))
}



#############################################################################################
#' @export
BVAR.dynSkew.Student.SV <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)

  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "dynSkew.Student", SV = TRUE)
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

  # new precompute
  i <- K*m
  theta.prior.prec <- V_b_prior_inv
  theta.prior.precmean <- theta.prior.prec %*% b_prior

  # Output
  mcmc <-  matrix(NA, nrow = m*K + 0.5*K*(K-1) + 1 + K + K + K + K + K * t_max + K*t_max + t_max,
                  ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # Sample B
    wt <- as.vector(1/(exp(h/2)*w_sqrt))
    y.tilde <- as.vector( A %*% (yt - Gamma_T * (w- mu.xi) ) ) * wt
    x.tilde <- kronecker( t(xt), A ) * wt
    theta.prec.chol <- chol( theta.prior.prec + crossprod(x.tilde) )
    theta <- backsolve( theta.prec.chol,
                        backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                   upper.tri = T, transpose = T )
                        + rnorm(K*m) )
    B <- matrix(theta,K,m)
    b_sample <- as.vector(B)
    # gamma <- theta[(m*K+1):((m+1)*K)]
    # D <- diag(gamma, K)

    # Sample gamma : carterkohn(y,Z,Ht)
    # y: n*T
    # Z: (nT)*(n_max*K)
    # Ht: (nT)*n
    # ysim <-  Z %*% B_t + t(chol(Hsim)) %*% rnorm(M)
    Sigma_T <- matrix(0, nrow = K * t_max, ncol = K)
    A_inv <- solve(A)
    for (i in c(1:t_max)){
      Sigma_T[ ((i-1)*K + 1):(i*K) , ] <- A_inv %*% diag(exp(h[,i]), nrow = K) %*% t(A_inv) * w_sample[i]
    }
    Gamma_filter <- fatBVARS::carterkohn(yt - B %*% xt, (w_sample - mu.xi) %x% diag(K),
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
    ytilde <- A%*% ((yt - B %*% xt  - Gamma_T * (w-mu.xi) )/w_sqrt)
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h

    # Sample A0
    u_std <- (yt - B %*% xt - Gamma_T * (w- mu.xi))/ w_sqrt # change from Gaussian
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
    q2 <- colSums( ( ( A %*% ( yt - B%*%xt + mu.xi * Gamma_T) ) / sqrtvol )^2 )
    p2 <- colSums( ( c( A %*% Gamma_T ) / sqrtvol )^2 )
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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, nu, diag(sigma_h), h0, diag(sigma_G), gamma0, as.numeric(Gamma_T), as.numeric(h), as.numeric(w_sample))
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
                        paste("nu"),
                        paste("sigma_h",c(1:K), sep = ""),
                        paste("lh0",c(1:K), sep = ""),
                        paste("sigma_G",c(1:K), sep = ""),
                        paste("G0",c(1:K), sep = ""),
                        sprintf("gamma_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        sprintf("h_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)),
                        paste("w",c(1:t_max), sep = ""))

  return(as.mcmc(t(mcmc)))
}

