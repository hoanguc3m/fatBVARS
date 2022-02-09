
# Center the W gamma - mean(W) gamma
# Does not seem to improve the mixing
BVAR.OST.SV.center <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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
  logsigma_nu <- rep(0,K)
  acount_nu <- rep(0,K)
  acount_w <- rep(0, t_max)
  gamma <- inits$gamma
  D <- diag(gamma)
  mu.xi <- nu/( nu - 2 )

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

  nu.prop.std <- rep(3,K)
  nu_temp <- nu

  # Output
  mcmc <-  matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K + K*t_max + K*t_max,
                  ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){


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
    # ff <- function(nu,k) {
    #   x <- dgamma(nu, shape = nu_gam_a, rate = nu_gam_b, log = T) +
    #               sum(dinvgamma(w[k,], shape = nu*0.5, rate = nu*0.5, log = T))
    #   return(x)
    # }
    tmpx <- nu
    for ( k in 1:K ) {
      tmpx[k] <- (nu_gam_a - 1 + 0.5 * t_max) / ( nu_gam_b + 0.5 * (sum(1/w[k,]) + sum(log(w[k,])) - t_max) ) # mode

      if ( runif(1) < 0.1 ) {
        nu.prop.std[k] <- tmpx[k] / sqrt(nu_gam_a - 1 + 0.5 * t_max) * 4 # big move, sd at mode
      } else {
        nu.prop.std[k] <- tmpx[k] / sqrt(nu_gam_a - 1 + 0.5 * t_max) # sd at mode
      }
      nu_temp[k] <- tmpx[k] + rnorm(1)*nu.prop.std[k]

    }
    # }
    # nu_temp = nu + exp(logsigma_nu)*rnorm(K)
    for (k in c(1:K)){
      if (nu_temp[k] > 8 && nu_temp[k] < 100){
        #      if (nu_temp[k] >= 5 && nu_temp[k] < 100){
        num_mh = dgamma(nu_temp[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w[k,], shape = nu_temp[k]*0.5, rate = nu_temp[k]*0.5, log = T))
        denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w[k,], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))

        propn_mh <- dnorm(nu[k], mean = tmpx[k], sd = nu.prop.std[k], log = T)
        propd_mh <- dnorm(nu_temp[k], mean = tmpx[k], sd = nu.prop.std[k], log = T)

        alpha = num_mh - denum_mh + propn_mh - propd_mh;
        temp = log(runif(1));
        if (alpha > temp){
          nu[k] = nu_temp[k]
          acount_nu[k] = acount_nu[k] + 1
        }

      }
    }
    mu.xi <- nu / ( nu - 2 )

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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, gamma, nu, diag(sigma_h), h0, as.numeric(h), as.numeric(w))
    if (j %% 1000 == 0) {
      #      cat(" Iteration- ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2) , " ", round(gamma,2) ,  " \n")
      cat(" Iteration- ", j, " ", round(acount_nu/j,2)," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2) , " ", round(gamma,2) ,  " \n")
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

nu.cond.post <- function(nu,lambda,alpha,n) {
  logpost <- log(nu)*( alpha - 1 ) + log(nu/2)*n*nu/2  - nu*lambda - n*lgamma(nu/2)
  return(logpost)
}


atau_post <- function(atau=a_tau,lambda2=lambda2_tau,thetas=theta,k=length(theta),rat=1){
  logpost <- sum(dgamma(thetas,atau,(atau*lambda2/2),log=TRUE))+dexp(atau,rate=rat,log=TRUE)
  return(logpost)
}

# BVAR.Gaussian.SV using Huber and Feldkircher prior: Adaptive Shrinkage in Bayesian Vector Autoregressive Models
BVAR.Gaussian.SV.HB <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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

  # Hierachical priors
  psi_b <- rep(1, K*m) # local shrinkage
  psi_a <- rep(10, 0.5*K*(K-1)) # local shrinkage

  varsilon_b <- 0.1 # Initial
  varsilon_a <- 0.005 # Initial

  lambda_b <- 0.1 # global shrinkage
  lambda_a <- 0.3 # global shrinkage


  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K*t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)

  accept <- 0
  accept2 <- 0
  scale1 <- .25
  scale2 <- .45
  c_tau=.1; d_tau=.1; # Hyper parameter

  for (j in c(1:samples)){

    # new precompute change here
    V_b_prior_inv <- diag( 1/ psi_b) # New hierarchical prior
    theta.prior.precmean <- V_b_prior_inv %*% b_prior

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


    # Sample A0
    u_std <- (yt - B %*%xt)
    u_neg <- - u_std
    a_sample <- rep(0, K * (K - 1) /2)
    V_a_prior <- diag( psi_a) # New hierarchical prior
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


    # shrink.type=="global"
    #Step III: Sample the prior scaling factors from GIG
    psi_b <- mapply(GIGrvg::rgig, n=1, lambda=varsilon_b-0.5, psi = varsilon_b * lambda_b, chi = (b_sample-b_prior)^2)
    psi_b[psi_b < 1e-9] <- 1e-9

    #Step IV: Sample varsilon_b through a simple RWMH step
    varsilon_b_prop <- exp(rnorm(1,0,scale1))*varsilon_b

    # post.diff <- log(varsilon_b_prop)-log(varsilon_b) +  # Jacob
    #       (K * m) * ( - lgamma(varsilon_b_prop) + lgamma(varsilon_b)) +
    #       (K * m) * varsilon_b_prop * log(varsilon_b_prop * lambda_b * 0.5) -
    #       (K * m) * varsilon_b * log(varsilon_b * lambda_b * 0.5) +
    #       (varsilon_b_prop - varsilon_b) * sum(log(psi_b)) -
    #       (varsilon_b_prop - varsilon_b) * (lambda_b * 0.5 * sum(psi_b)) +
    #       (varsilon_b - varsilon_b_prop)
    #
    # sum(dgamma(psi_b, shape = varsilon_b_prop, rate = (varsilon_b_prop*lambda_b/2),log=TRUE)) -
    #   sum(dgamma(psi_b, shape = varsilon_b, rate = (varsilon_b*lambda_b/2),log=TRUE)) +
    #   dexp(varsilon_b_prop,rate=1,log=TRUE) - dexp(varsilon_b,rate=1,log=TRUE) + log(varsilon_b_prop) - log(varsilon_b)

    post_a_prop <- atau_post(atau=varsilon_b_prop,thetas = psi_b,lambda2 = lambda_b)
    post_a_old <- atau_post(atau=varsilon_b,thetas = psi_b,lambda2 = lambda_b)
    post.diff <- post_a_prop-post_a_old+log(varsilon_b_prop)-log(varsilon_b)

    post.diff <- ifelse(is.nan(post.diff),-Inf,post.diff)
    if (post.diff > log(runif(1,0,1))){
      varsilon_b <- varsilon_b_prop
      accept <- accept+1
    }

    #Step V: Sample hyperparameter lambda from G(a,b)
    lambda_b <- rgamma(1, shape = c_tau + varsilon_b * K*m,
                       rate = d_tau + varsilon_b/2*sum(psi_b))
    #############
    #
    #     #Step III: Sample the prior scaling factors from GIG
    #     psi_a <- mapply(GIGrvg::rgig, n=1, lambda=varsilon_a-0.5, psi = varsilon_a * lambda_a, chi = a_sample^2)
    #     psi_a[psi_a < 1e-9] <- 1e-9
    #     # psi_a <- mapply(GIGrvg::rgig, n=1, lambda=0-0.5, psi = 0, chi = a_sample^2)
    #
    #     #Step IV: Sample varsilon_a through a simple RWMH step
    #     varsilon_a_prop <- exp(rnorm(1,0,scale2))*varsilon_a
    #
    #     # post.diff <- log(varsilon_a_prop)-log(varsilon_a) +  # Jacob
    #     #   (K * (K-1) * 0.5) * ( - lgamma(varsilon_a_prop) + lgamma(varsilon_a)) +
    #     #   (K * (K-1) * 0.5) * varsilon_a_prop * log(varsilon_a_prop * lambda_a * 0.5) -
    #     #   (K * (K-1) * 0.5) * varsilon_a * log(varsilon_a * lambda_a * 0.5) +
    #     #   (varsilon_a_prop - varsilon_a) * sum(log(psi_a)) -
    #     #   (varsilon_a_prop - varsilon_a) * (lambda_a * 0.5 * sum(psi_a) + 1)
    #     # No idea dlnorm
    #     post_a_prop <- atau_post(atau=varsilon_a_prop,thetas = psi_a,lambda2 = lambda_a)  # + dlnorm(varsilon_a,varsilon_a_prop,scale2,log=TRUE)
    #     post_a_old <- atau_post(atau=varsilon_a,thetas = psi_a,lambda2 = lambda_a)  # + dlnorm(varsilon_a_prop,varsilon_a,scale2,log=TRUE)
    #     post.diff <- post_a_prop-post_a_old+log(varsilon_a_prop)-log(varsilon_a)
    #
    #     post.diff <- ifelse(is.nan(post.diff),-Inf,post.diff)
    #     if (post.diff > log(runif(1,0,1))){
    #       varsilon_a <- varsilon_a_prop
    #       accept2 <- accept2 + 1
    #     }
    #
    #     #Step V: Sample hyperparameter lambda from G(a,b)
    #     lambda_a <- rgamma(1, shape = c_tau + varsilon_a * (K * (K-1) * 0.5),
    #                        rate = d_tau + varsilon_a/2*sum(psi_a))

    batchlength = 100;
    if(j < inits$burnin){

      if ((accept/j)>0.3) scale1 <- 1.01*scale1
      if ((accept/j)<0.15) scale1 <- 0.99*scale1
      if ((accept2/j)>0.3)  scale2 <- 1.01*scale2
      if ((accept2/j)<0.15) scale2 <- 0.99*scale2
    }

    if(j %% 1000 == 0){
      cat(" Iteration ", j, " accept =  ", accept, " scale1 =  ", scale1, " accept2 =  ", accept2, " scale2 =  ", scale2, " \n")
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
  return(as.mcmc(t(mcmc)))
}


BVAR.MST.SV.HB <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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
  gamma <- inits$gamma
  D <- diag(gamma)

  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)
  #w_sqrt_inv <- 1/w_sqrt

  # Hierachical priors
  psi_b <- rep(1, K*m) # local shrinkage
  psi_a <- rep(10, 0.5*K*(K-1)) # local shrinkage

  varsilon_b <- 0.1 # Initial
  varsilon_a <- 0.005 # Initial

  lambda_b <- 0.1 # global shrinkage
  lambda_a <- 0.3 # global shrinkage

  # new precompute
  i <- K*m
  theta.prior.prec = matrix(0,i+K,i+K)
  theta.prior.prec[1:i,1:i] <- V_b_prior_inv
  theta.prior.prec[i+(1:K),i+(1:K)] <- solve( V_gamma_prior )
  theta.prior.precmean <- theta.prior.prec %*% c( b_prior, gamma_prior )

  # Output
  mcmc <-  matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K + K*t_max + K*t_max,
                  ncol = (samples - inits$burnin)%/% inits$thin)

  accept <- 0
  accept2 <- 0
  scale1 <- .25
  scale2 <- .45
  c_tau=.01; d_tau=.01; # Hyper parameter

  for (j in c(1:samples)){

    # new precompute change here
    i <- K*m
    theta.prior.prec = matrix(0,i+K,i+K)
    theta.prior.prec[1:i,1:i] <- diag( 1/ psi_b) # New hierarchical prior
    theta.prior.prec[i+(1:K),i+(1:K)] <- solve( V_gamma_prior )
    theta.prior.precmean <- theta.prior.prec %*% c( b_prior, gamma_prior )


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
      W.mat[idx + (i-1)*t_max*K + i] <- sqrt(w[i,]) # W^(1/2)
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
    ytilde <- A%*% ((yt - B %*%xt - D %*% w)/ w_sqrt)
    aux <- sample_h_ele(ytilde = ytilde, sigma_h = sigma_h, h0_mean = h0_mean, h = h, K = K, t_max = t_max, prior = prior)
    h <- aux$Sigtdraw
    h0 <- as.numeric(aux$h0)
    sqrtvol <- aux$sigt
    sigma_h <- aux$sigma_h


    # Sample A0
    u_std <- (yt - B %*% xt - D%*% w) / w_sqrt # change from Gaussian
    u_neg <- - u_std
    a_sample <- rep(0, K * (K - 1) /2)
    V_a_prior <- diag( psi_a) # New hierarchical prior
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



    # shrink.type=="global"
    #Step III: Sample the prior scaling factors from GIG
    psi_b <- mapply(GIGrvg::rgig, n=1, lambda=varsilon_b-0.5, psi = varsilon_b * lambda_b, chi = (b_sample-b_prior)^2)
    psi_b[psi_b < 1e-9] <- 1e-9

    #Step IV: Sample varsilon_b through a simple RWMH step
    varsilon_b_prop <- exp(rnorm(1,0,scale1))*varsilon_b

    # post.diff <- log(varsilon_b_prop)-log(varsilon_b) +  # Jacob
    #       (K * m) * ( - lgamma(varsilon_b_prop) + lgamma(varsilon_b)) +
    #       (K * m) * varsilon_b_prop * log(varsilon_b_prop * lambda_b * 0.5) -
    #       (K * m) * varsilon_b * log(varsilon_b * lambda_b * 0.5) +
    #       (varsilon_b_prop - varsilon_b) * sum(log(psi_b)) -
    #       (varsilon_b_prop - varsilon_b) * (lambda_b * 0.5 * sum(psi_b)) +
    #       (varsilon_b - varsilon_b_prop)
    #
    # sum(dgamma(psi_b, shape = varsilon_b_prop, rate = (varsilon_b_prop*lambda_b/2),log=TRUE)) -
    #   sum(dgamma(psi_b, shape = varsilon_b, rate = (varsilon_b*lambda_b/2),log=TRUE)) +
    #   dexp(varsilon_b_prop,rate=1,log=TRUE) - dexp(varsilon_b,rate=1,log=TRUE) + log(varsilon_b_prop) - log(varsilon_b)

    post_a_prop <- atau_post(atau=varsilon_b_prop,thetas = psi_b,lambda2 = lambda_b)
    post_a_old <- atau_post(atau=varsilon_b,thetas = psi_b,lambda2 = lambda_b)
    post.diff <- post_a_prop-post_a_old+log(varsilon_b_prop)-log(varsilon_b)

    post.diff <- ifelse(is.nan(post.diff),-Inf,post.diff)
    if (post.diff > log(runif(1,0,1))){
      varsilon_b <- varsilon_b_prop
      accept <- accept+1
    }

    #Step V: Sample hyperparameter lambda from G(a,b)
    lambda_b <- rgamma(1, shape = c_tau + varsilon_b * K*m,
                       rate = d_tau + varsilon_b/2*sum(psi_b))
    #############

    # #Step III: Sample the prior scaling factors from GIG
    # psi_a <- mapply(GIGrvg::rgig, n=1, lambda=varsilon_a-0.5, psi = varsilon_a * lambda_a, chi = a_sample^2)
    # psi_a[psi_a < 1e-9] <- 1e-9
    #
    #
    # #Step IV: Sample varsilon_a through a simple RWMH step
    # varsilon_a_prop <- exp(rnorm(1,0,scale2))*varsilon_a
    #
    # # post.diff <- log(varsilon_a_prop)-log(varsilon_a) +  # Jacob
    # #   (K * (K-1) * 0.5) * ( - lgamma(varsilon_a_prop) + lgamma(varsilon_a)) +
    # #   (K * (K-1) * 0.5) * varsilon_a_prop * log(varsilon_a_prop * lambda_a * 0.5) -
    # #   (K * (K-1) * 0.5) * varsilon_a * log(varsilon_a * lambda_a * 0.5) +
    # #   (varsilon_a_prop - varsilon_a) * sum(log(psi_a)) -
    # #   (varsilon_a_prop - varsilon_a) * (lambda_a * 0.5 * sum(psi_a) + 1)
    # # No idea dlnorm
    # post_a_prop <- atau_post(atau=varsilon_a_prop,thetas = psi_a,lambda2 = lambda_a)  + dlnorm(varsilon_a,varsilon_a_prop,scale2,log=TRUE)
    # post_a_old <- atau_post(atau=varsilon_a,thetas = psi_a,lambda2 = lambda_a)  + dlnorm(varsilon_a_prop,varsilon_a,scale2,log=TRUE)
    # post.diff <- post_a_prop-post_a_old+log(varsilon_a_prop)-log(varsilon_a)
    #
    # post.diff <- ifelse(is.nan(post.diff),-Inf,post.diff)
    # if (post.diff > log(runif(1,0,1))){
    #   varsilon_a <- varsilon_a_prop
    #   accept2 <- accept2 + 1
    # }
    #
    # #Step V: Sample hyperparameter lambda from G(a,b)
    # lambda_a <- rgamma(1, shape = c_tau + varsilon_a * (K * (K-1) * 0.5),
    #                    rate = d_tau + varsilon_a/2*sum(psi_a))

    batchlength = 10;
    TARGACCEPT = 0.3
    if(j < inits$burnin){

      if ((accept/j)>0.3) scale1 <- 1.01*scale1
      if ((accept/j)<0.15) scale1 <- 0.99*scale1
      if ((accept2/j)>0.3)  scale2 <- 1.01*scale2
      if ((accept2/j)<0.15) scale2 <- 0.99*scale2
    }

    if(j %% 1000 == 0){
      cat(" Iteration ", j, " accept =  ", accept, " scale1 =  ", scale1, " accept2 =  ", accept2, " scale2 =  ", scale2, " \n")
    }


    # Sample w
    u <- (yt - B %*%xt)
    u_proposal <- A %*% (yt - B %*%xt)
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

recursive_seperate.HB <- function(y, t_start = 100, t_pred = 12, K, p, dist = "MST", SV = T, outname = NULL, samples = 30000, burnin = 10000, thin = 1){
  t_max = nrow(y)
  if (is.null(outname)) outname = paste("Recursive_", dist, "_", SV, "_", t_start+1, ".RData", sep = "")
  time_current <- t_start # index of the longer y
  t_current <- time_current - p # length of the shorter y
  y_current <- matrix(y[c(1:time_current), ], ncol = K)
  prior <- get_prior(tail(y_current, time_current - p), p = p, dist=dist, SV = SV)
  inits <- get_init(prior, samples = samples, burnin = burnin, thin = thin)

  if (dist == "Gaussian") {
    mcmc <- BVAR.Gaussian.SV.HB(y = tail(y_current, time_current - p), K = K, p = p,
                                y0 = head(y_current, p), prior = prior, inits = inits)
  }
  if (dist == "MST") {
    mcmc <- BVAR.MST.SV.HB(y = tail(y_current, time_current - p), K = K, p = p,
                           y0 = head(y_current, p), prior = prior, inits = inits)
  }

  Chain <- list(mcmc = mcmc,
                y = tail(y_current, time_current - p),
                y0 = head(y_current, p),
                K = K,
                p = p,
                dist = dist,
                prior = prior,
                inits = inits)
  y_obs_future <- matrix(y[c((time_current+1):(time_current+t_pred)), ],ncol = K)
  forecast_err <- forecast_density(Chain = Chain, y_obs_future = y_obs_future, t_current = t_current)
  out_recursive <- list(time_id = time_current+1,
                        forecast_err = forecast_err,
                        dist = dist,
                        SV = SV,
                        y_current = y_current,
                        y_obs_future = y_obs_future)
  save(out_recursive, file = outname)
  return( out_recursive)
}
