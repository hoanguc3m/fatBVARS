#' Bayesian inference of VAR model without SV
#'
#' Bayesian inference of VAR model without SV
#' \deqn{y_t = B x_t + SQRT(w_t) A^(-1) Sigma eps_t}
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
#' datagen <- sim.VAR.novol(dist="Gaussian")
#' y <- datagen$y
#' prior <- get_prior(y, p = 2, dist="Gaussian", SV = F)
#' inits <- get_init(prior)
#' Chain1 <- BVAR.novol(y, K = 5, p = 2, dist = "Gaussian", y0 = NULL, prior = prior, inits = inits)
#' plot(Chain1)
#' }
BVAR.novol <- function(y, K, p, dist, y0 = NULL, prior = NULL, inits = NULL){
  if (!(dist %in% c("Gaussian","Student","Skew.Student",
                    "MT","MST",
                    "OT","OST",
                    "dynSkew.Student", "OrthSkewNorm", "dynOST") ))
    stop("dist is not implemented.")

  if (prior$SV == FALSE){
    Start = Sys.time()
    if (dist == "Gaussian") Chain <- BVAR.Gaussian.novol(y, K, p, y0, prior, inits)
    if (dist == "Student") Chain <- BVAR.Student.novol(y, K, p, y0, prior, inits)
    # if (dist == "Skew.Student") Chain <- BVAR.Skew.Student.novol(y, K, p, y0, prior, inits)
    if (dist == "Skew.Student") Chain <- BVAR.Skew.Student.novol(y, K, p, y0, prior, inits)
    if (dist == "MT") Chain <- BVAR.MT.novol(y, K, p, y0, prior, inits)
    if (dist == "MST") Chain <- BVAR.MST.novol(y, K, p, y0, prior, inits)
    if (dist == "OT") Chain <- BVAR.OT.novol(y, K, p, y0, prior, inits)
    if (dist == "OST") Chain <- BVAR.OST.novol(y, K, p, y0, prior, inits)

    if (dist == "dynSkew.Student") Chain <- BVAR.dynSkew.Student.novol(y, K, p, y0, prior, inits)
    if (dist == "OrthSkewNorm") Chain <- BVAR.OrthSkewNorm.novol(y, K, p, y0, prior, inits)
    if (dist == "dynOST") Chain <- BVAR.dynOST.novol(y, K, p, y0, prior, inits)
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

  } else {
    warning("prior$SV is TRUE")
  }
  return(out)
}

###########################################################################
#' @export
BVAR.Gaussian.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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
  samples <- inits$samples
  A <- inits$A0
  B <- inits$B0
  sigma <- inits$sigma
  # Sigma <- solve(inits$A0) %*% diag(sigma, nrow = length(sigma))
  # Sigma2 <- Sigma %*% t(Sigma)
  # Sigma2_inv <- solve(Sigma2)
  V_b_prior_inv <- solve(V_b_prior)
  # new precompute
  theta.prior.precmean <- V_b_prior_inv %*% b_prior

  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # # Sample B
    # V_b_post_inv <- V_b_prior_inv + kronecker(xt %*% t(xt),Sigma2_inv)
    # b_post <- kronecker(xt, Sigma2_inv) %*% vec(yt)
    #
    # V_b_post <- solve(V_b_post_inv)
    # b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    # b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    # B <- Vec_to_Mat(b_sample, K,p)

    # sample B
    A.tmp <- diag(1/sigma, K) %*% A
    y.tilde <- as.vector( A.tmp %*% yt )
    x.tilde <- kronecker( t(xt), A.tmp )
    theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
    b_sample <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T )
                         + rnorm(K*m) )
    B <- matrix(b_sample,K,m)

    # Sample sigma
    sigma2 <- rep(0,K)
    u_ort <- A %*% (yt - B %*%xt)
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
#    sigma_post_b <- sigma0_S0 + apply(u_ort^2, 1, sum)
    sigma_post_b <- sigma0_S0 + rowSums(u_ort^2)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)

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

    # Sigma <- solve(A) %*% diag(sigma, nrow = length(sigma))
    # Sigma2 <- Sigma %*% t(Sigma)
    # Sigma2_inv <- solve(Sigma2)
    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma)
    if (j %% 1000 == 0) { cat(" Iteration ", j, " \n")}
  }
  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("sigma",c(1:K), sep = ""))
  return(as.mcmc(t(mcmc)))
}

###########################################################################
#' @export
BVAR.Student.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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

  samples <- inits$samples
  A <- inits$A0
  B <- Vec_to_Mat(inits$B0, K, p)
  sigma <- inits$sigma
  # Sigma <- solve(inits$A0) %*% diag(inits$sigma, nrow = K)
  # Sigma2 <- Sigma %*% t(Sigma)
  # Sigma2_inv <- solve(Sigma2)
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

  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + 1 + t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # # Sample B
    # xt_G <- xt / reprow(sqrt(w_sample), m)
    # yt_G <- yt / w_sqrt
    # V_b_post_inv <- V_b_prior_inv + kronecker(xt_G %*% t(xt_G),Sigma2_inv)
    # b_post <- kronecker(xt_G, Sigma2_inv) %*% vec(yt_G)
    #
    # V_b_post <- solve(V_b_post_inv)
    # b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    # b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    # B <- Vec_to_Mat(b_sample, K,p)

    # sample B
    A.tmp <- diag(1/sigma, K) %*% A
    wt <- as.vector(1/w_sqrt)
    y.tilde <- as.vector( A.tmp %*% yt ) * wt
    x.tilde <- kronecker( t(xt), A.tmp ) * wt
    theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
    b_sample <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T )
                           + rnorm(K*m) )
    B <- matrix(b_sample,K,m)

    # Sample sigma
    sigma2 <- rep(0,K)
    u_ort <- A %*% ((yt - B %*%xt)/ w_sqrt)    # change from Gaussian
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
#    sigma_post_b <- sigma0_S0 + apply(u_ort^2, 1, sum)
    sigma_post_b <- sigma0_S0 + rowSums(u_ort^2)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)

    # Sample A0
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

    # Sigma <- solve(A) %*% diag(sigma, nrow = length(sigma))
    # Sigma2 <- Sigma %*% t(Sigma)
    # Sigma2_inv <- solve(Sigma2)
    # Sigma2_inv <- (Sigma2_inv + t(Sigma2_inv))*0.5

    # Sample w
    u <- (yt - B %*%xt)
    shape_w <- nu*0.5 + K*0.5
    # rate_w <- as.numeric(nu*0.5 + 0.5 * apply((A %*% u / sigma)^2, MARGIN = 2, FUN = sum))
    # for (i in c(1:t_max)){
    #   #w_sample[i] <- rinvgamma(1, shape = shape_w, rate = rate_w[i])
    #   w_sample[i] <- 1/rgamma(1, shape = shape_w, rate = rate_w[i])
    # }
    rate_w <- as.numeric(nu*0.5 + 0.5 * colSums((A %*% u / sigma)^2))
    w_sample <- rinvgamma( n=t_max, shape = shape_w, rate = rate_w )
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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, nu, as.numeric(w_sample))
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
                        paste("sigma",c(1:K), sep = ""),
                        paste("nu"),
                        paste("w",c(1:t_max), sep = ""))

  return(as.mcmc(t(mcmc)))
}


#############################################################################################
#' @export
BVAR.Skew.Student.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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
  gamma <- inits$gamma
  D <- diag(gamma, nrow = length(gamma))

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
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + 1 + t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # # Sample B
    # xt_G <- xt / reprow(sqrt(w_sample), m)
    # yt_G <- (yt - D %*% w)  / w_sqrt
    # V_b_post_inv <- V_b_prior_inv + kronecker(xt_G %*% t(xt_G),Sigma2_inv)
    # b_post <- kronecker(xt_G, Sigma2_inv) %*% vec(yt_G)
    #
    # V_b_post <- solve(V_b_post_inv)
    # b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    # b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    # B <- Vec_to_Mat(b_sample, K,p)
    #
    # # Sample gamma
    # V_gamma_post_inv = solve(V_gamma_prior) + sum(w_sample) * Sigma2_inv
    # gamma_post = apply( Sigma2_inv %*% (yt - B %*% xt), MARGIN = 1, FUN = sum)
    # V_gamma_post <- solve(V_gamma_post_inv)
    # gamma_post <- V_gamma_post %*% ( solve(V_gamma_prior) %*% gamma_prior + gamma_post)
    # gamma <- as.numeric(gamma_post + t(chol(V_gamma_post)) %*% rnorm(K))
    # D <- diag(gamma, nrow = length(gamma))

    # Sample B and gamma
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
    D <- diag(gamma, K)

    # Sample sigma
    sigma2 <- rep(0,K)
    u_ort <- A %*% ((yt - B %*% xt - D %*% (w-mu.xi))/ w_sqrt)    # change from Gaussian
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
#    sigma_post_b <- sigma0_S0 + apply(u_ort^2, 1, sum)
    sigma_post_b <- sigma0_S0 + rowSums(u_ort^2)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)

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


    # Sample w
    q2 <- colSums( ( ( A %*% ( yt - B%*%xt + mu.xi * gamma) ) / sigma )^2 )
    p2 <- sum( ( c( A %*% gamma ) / sigma )^2 )
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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, gamma, nu, as.numeric(w_sample))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2), " ", round(gamma,2), " \n")
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
                        paste("nu"),
                        paste("w",c(1:t_max), sep = ""))

  return(as.mcmc(t(mcmc)))
}

###########################################################################
#' @export
BVAR.MT.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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
  samples <- inits$samples
  A <- inits$A0
  B <- Vec_to_Mat(inits$B0, K, p)
  sigma <- inits$sigma
  # Sigma <- solve(inits$A0) %*% diag(sigma, nrow = length(sigma))
  # Sigma2 <- Sigma %*% t(Sigma)
  # Sigma2_inv <- solve(Sigma2)
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

  # new precompute
  theta.prior.precmean <- V_b_prior_inv %*% b_prior

  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K*t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # # Sample B
    # V_b_post_inv <- V_b_prior_inv +
    #   Reduce( f = "+",
    #           x = lapply(1:t_max, function(i) kronecker(xt[,i] %*% t(xt[,i]), (w_sqrt_inv[,i] %*% t(w_sqrt_inv[,i]))) ),
    #           accumulate = FALSE) * (kronecker(matrix(1,m,m), Sigma2_inv))
    #
    # b_post <- Reduce(f = "+",
    #                   x = lapply(1:t_max, function(i) kronecker(xt[,i], ( (w_sqrt_inv[,i] %*% t(w_sqrt_inv[,i])) * Sigma2_inv ) %*% yt[,i])),
    #                   accumulate = FALSE)
    #
    # V_b_post <- solve(V_b_post_inv)
    # b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    # b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    # B <- Vec_to_Mat(b_sample, K,p)

    # Sample B
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

    # Sample sigma
    sigma2 <- rep(0,K)
    u_ort <- A %*% ((yt - B %*%xt)/ w_sqrt)    # change from Gaussian
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
    # sigma_post_b <- sigma0_S0 + apply(u_ort^2, 1, sum)
    sigma_post_b <- sigma0_S0 + rowSums(u_ort^2)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)

    # Sample A0
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
#    log_w_temp <- log(w_temp)
#    w_temp_sqrt_inv = 1 / sqrt(w_temp)
    # prior_num <- apply(dinvgamma(w_temp, shape = nu*0.5, rate = nu*0.5, log = T), MARGIN = 2, FUN = sum)
    # prior_denum <- apply(dinvgamma(w, shape = nu*0.5, rate = nu*0.5, log = T), MARGIN = 2, FUN = sum)
    # proposal_num <- apply(dinvgamma(w_temp, shape = a_target, rate = b_target, log = T), MARGIN = 2, FUN = sum)
    # proposal_denum <- apply(dinvgamma(w, shape = a_target, rate = b_target, log = T), MARGIN = 2, FUN = sum)
    prior_num <- colSums(dinvgamma(w_temp, shape = nu*0.5, rate = nu*0.5, log = T))
    prior_denum <- colSums(dinvgamma(w, shape = nu*0.5, rate = nu*0.5, log = T))
    proposal_num <- colSums(dinvgamma(w_temp, shape = a_target, rate = b_target, log = T))
    proposal_denum <- colSums(dinvgamma(w, shape = a_target, rate = b_target, log = T))

    # num_mh <- prior_num - proposal_num -
    #   0.5 * apply(log_w_temp, MARGIN = 2, FUN = sum) -
    #   0.5 * sapply(1:t_max, function(i) t(u[,i]) %*% ((w_temp_sqrt_inv[,i] %*% t(w_temp_sqrt_inv[,i])) * Sigma2_inv) %*% u[,i])
    A.w.u.s.prop <- ( A %*% ( u/w_temp_sqrt ) )*( 1/ sigma )
    num_mh <- prior_num - proposal_num - 0.5 * colSums(log(w_temp)) - 0.5 * colSums(A.w.u.s.prop^2)
    A.w.u.s.curr <- ( A %*% ( u/w_sqrt ) )*( 1/ sigma )
    denum_mh <- prior_denum - proposal_denum - 0.5 * colSums(log(w)) - 0.5 * colSums(A.w.u.s.curr^2)
    # denum_mh <- prior_denum - proposal_denum -
    #   0.5 * apply(log(w), MARGIN = 2, FUN = sum) -
    #   0.5 * sapply(1:t_max, function(i) t(u[,i]) %*% ((w_sqrt_inv[,i] %*% t(w_sqrt_inv[,i])) * Sigma2_inv)  %*% u[,i])
    alpha = num_mh - denum_mh
    temp = log(runif(t_max))
    acre = alpha > temp
    w[,acre] <- w_temp[,acre]
    w_sqrt[,acre] <- w_temp_sqrt[,acre]
#    w_sqrt_inv <- 1/w_sqrt
    acount_w[acre] <- acount_w[acre] + 1

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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, nu, as.numeric(w))
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
                        paste("sigma",c(1:K), sep = ""),
                        paste("nu",c(1:K), sep = ""),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K))
  )

  return(as.mcmc(t(mcmc)))
}

###########################################################################

#' @export
BVAR.MST.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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
  samples <- inits$samples
  A <- inits$A0
  B <- Vec_to_Mat(inits$B0, K, p)
  sigma <- inits$sigma
  # Sigma <- solve(inits$A0) %*% diag(sigma, nrow = length(sigma))
  # Sigma2 <- Sigma %*% t(Sigma)
  # Sigma2_inv <- solve(Sigma2)
  V_b_prior_inv <- solve(V_b_prior)

  # Multi degrees of freedom
  nu <- inits$nu
  mu.xi <- nu/( nu - 2 )
  logsigma_nu <- rep(0,K)
  acount_nu <- rep(0,K)
  acount_w <- rep(0, t_max)
  gamma <- inits$gamma
  D <- diag(gamma, K)

  # Init w as Gaussian
  w_sample <- rep(1, t_max)
  w <- reprow(w_sample, K)
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
  for (j in c(1:samples)){

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
    D <- diag(gamma, K)


    # Sample sigma
    sigma2 <- rep(0,K)
    u_ort <- A %*% ((yt - B %*% xt - D %*% ( w - mu.xi ))/ w_sqrt)    # change from Gaussian
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
    sigma_post_b <- sigma0_S0 + rowSums(u_ort^2)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)

    # Sample A0
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

    # Sigma <- solve(A) %*% diag(sigma, nrow = length(sigma))
    # Sigma2 <- Sigma %*% t(Sigma)
    # Sigma2_inv <- solve(Sigma2)
    # Sigma2_inv <- (Sigma2_inv + t(Sigma2_inv))*0.5
    Sigma2_inv <- t(A)%*%diag(1/sigma2)%*%A

    # Sample w
    u <- (yt - B %*%xt + mu.xi * gamma)
    u_proposal <- A %*% (yt - B %*%xt) + mu.xi * gamma

    a_target <- (nu*0.5 + 1*0.5) * 0.75 # adjust by 0.75
    b_target <- (nu*0.5 + 0.5 * u_proposal^2 / sigma^2) * 0.75 # adjust by 0.75
    w_temp <- matrix(rinvgamma( n = K*t_max, shape = a_target, rate = b_target), nrow = K)
    log_w_temp <- log(w_temp)
    w_temp_sqrt <- sqrt(w_temp)
    #    w_temp_sqrt_inv = 1 / sqrt(w_temp)
    # prior_num <- colSums(dinvgamma(w_temp, shape = nu*0.5, rate = nu*0.5, log = T))
    # prior_denum <- colSums(dinvgamma(w, shape = nu*0.5, rate = nu*0.5, log = T))
    # proposal_num <- colSums(dinvgamma(w_temp, shape = a_target, rate = b_target, log = T))
    # proposal_denum <- colSums(dinvgamma(w, shape = a_target, rate = b_target, log = T))
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


    #     num_mh <- prior_num - proposal_num -
    #       0.5 * colSums(log_w_temp) -
    #       0.5 * sapply(1:t_max, function(i) { u.gw <- w_temp_sqrt_inv[,i]*(u[,i]- gamma*w_temp[,i]); t(u.gw) %*% Sigma2_inv %*% u.gw } )
    #     #  0.5 * sapply(1:t_max, function(i) t(u[,i]- gamma*w_temp[,i]) %*% ((w_temp_sqrt_inv[,i] %*% t(w_temp_sqrt_inv[,i])) * Sigma2_inv) %*% (u[,i]- gamma*w_temp[,i]))
    #     denum_mh <- prior_denum - proposal_denum -
    #       0.5 * colSums(log(w)) -
    #       0.5 * sapply(1:t_max, function(i) { u.gw <- (u[,i] - gamma*w[,i]) / w_sqrt[,i]; t(u.gw) %*% Sigma2_inv %*%  u.gw } )
    #     #  0.5 * sapply(1:t_max, function(i) t(u[,i] - gamma*w[,i]) %*% ((w_sqrt_inv[,i] %*% t(w_sqrt_inv[,i])) * Sigma2_inv) %*% (u[,i] - gamma*w[,i]))
    #     alpha = num_mh - denum_mh
    #     temp = log(runif(t_max))
    #     acre = alpha > temp
    #     w[,acre] <- w_temp[,acre]
    #     w_sqrt <- sqrt(w)
    # #    w_sqrt_inv <- 1/w_sqrt
    #     acount_w[acre] <- acount_w[acre] + 1

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
    mu.xi <- nu/( nu - 2 )

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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, gamma, nu, as.numeric(w))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2) , " ", round(gamma,2),  " \n")
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

  return(as.mcmc(t(mcmc)))
}

###########################################################################
#' @export
BVAR.OT.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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
  samples <- inits$samples
  A <- inits$A0
  B <- Vec_to_Mat(inits$B0, K, p)
  sigma <- inits$sigma
  sigma2 <- sigma^2

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
  w_sqrt_inv <- 1/w_sqrt

  # new precompute
  theta.prior.precmean <- V_b_prior_inv %*% b_prior

  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K*t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # # Sample B
    # b_post = rep(0, m*K)
    # V_b_post_inv = V_b_prior_inv
    # for (i in c(1:t_max)){
    #   V_b_post_inv <- V_b_post_inv + kronecker(xt[,i] %*% t(xt[,i]), t(A)%*% diag(1/sigma2 / w[,i] ) %*% A )
    #   b_post <- b_post + kronecker(xt[,i], (t(A)%*% diag(1/sigma2 / w[,i] ) %*% A) %*% yt[,i])
    # }
    #
    # V_b_post <- solve(V_b_post_inv)
    # b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    # b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    # B <- Vec_to_Mat(b_sample, K,p)

    # Sample B
    wt <- rep( 1/sigma, t_max ) / as.numeric(w_sqrt)
    y.tilde <- as.vector( A %*% yt ) * wt
    x.tilde <- kronecker( t(xt), A ) * wt
    theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
    b_sample <- backsolve( theta.prec.chol,
                           backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                      upper.tri = T, transpose = T )
                           + rnorm(K*m) )
    B <- matrix(b_sample,K,m)

    # Sample sigma
    u_ort <- (A %*% (yt - B %*%xt))/ w_sqrt    # change from Gaussian
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
    sigma_post_b <- sigma0_S0 + apply(u_ort^2, 1, sum)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)

    # Sample A0
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

    # Sample w
    u <-  A %*% (yt - B %*%xt)
    w <- matrix( rinvgamma( n = K*t_max, shape = (nu+1)*0.5, rate = 0.5*( nu + (u^2)/sigma2 ) ), K, t_max )
    # for (i in c(1:t_max)){
    #   w[,i] <- mapply(rinvgamma, n = 1, shape = nu*0.5 + 1*0.5, rate = nu*0.5 + 0.5 * (u[,i]^2) / sigma2)
    # }
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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, nu, as.numeric(w))
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
                        paste("sigma",c(1:K), sep = ""),
                        paste("nu",c(1:K), sep = ""),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K))
  )

  return(as.mcmc(t(mcmc)))
}

###########################################################################
#' @export
BVAR.OST.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
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
  samples <- inits$samples
  A <- inits$A0
  B <- Vec_to_Mat(inits$B0, K, p)
  sigma <- inits$sigma
  sigma2 <- sigma^2

  V_b_prior_inv <- solve(V_b_prior)

  # Multi degrees of freedom
  nu <- inits$nu
  mu.xi <- nu/( nu - 2 )
  logsigma_nu <- rep(0,K)
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
  w_sample <- rep(1, t_max)
  w <- reprow(w_sample, K)
  w_sqrt <- sqrt(w)
  w_sqrt_inv <- 1/w_sqrt
  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K*t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # # Sample B
    # b_post = rep(0, m*K)
    # V_b_post_inv = V_b_prior_inv
    # inv_A <- solve(A)
    # for (i in c(1:t_max)){
    #   V_b_post_inv <- V_b_post_inv + kronecker(xt[,i] %*% t(xt[,i]), t(A)%*% diag(1/sigma2 / w[,i] ) %*% A )
    #   b_post <- b_post + kronecker(xt[,i], (t(A)%*% diag(1/sigma2 / w[,i] ) %*% A) %*% (yt[,i] - inv_A %*% (gamma * w[,i])) )
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
    # V_gamma_post_inv <- solve(V_gamma_prior) + diag(apply(w, MARGIN = 1, FUN = sum)/ sigma2)
    # gamma_post <- apply(u, MARGIN = 1, FUN = sum) / sigma2
    #
    # V_gamma_post <- solve(V_gamma_post_inv)
    # gamma_post <- V_gamma_post %*% ( solve(V_gamma_prior) %*% gamma_prior + gamma_post)
    # gamma <- as.numeric(gamma_post + t(chol(V_gamma_post)) %*% rnorm(K))
    # D <- diag(gamma)

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
    theta <- backsolve( theta.prec.chol,
                        backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                   upper.tri = T, transpose = T )
                        + rnorm(K*(m+1)) )
    B <- matrix(theta[1:(m*K)],K,m)
    b_sample <- as.vector(B)
    gamma <- theta[(m*K+1):((m+1)*K)]
    D <- diag(gamma)

    # Sample sigma
    u_ort <- (A %*% (yt - B %*%xt) - D %*% ( w - mu.xi ))/ w_sqrt    # change from Gaussian
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
    sigma_post_b <- sigma0_S0 + apply(u_ort^2, 1, sum)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)

    # Sample A0
    u_std <- (yt - B %*%xt) # change from Gaussian
    u_neg <- - u_std
    a_sample <- rep(0, K * (K - 1) /2)
    if (K > 1) {
      for (i in c(2:K)){
        id_end <- i*(i-1)/2
        id_start <- id_end - i + 2
        a_sub <- a_prior[id_start:id_end]
        V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
        a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = (u_std[i,] - (w[i,]  - mu.xi[i] )  * gamma[i]) / sigma[i]/ w_sqrt[i,],
                                                     xsub = matrix(u_neg[1:(i-1),] / sigma[i] / reprow(w_sqrt[i,], i-1), nrow = i-1),
                                                     a_sub = a_sub,
                                                     V_a_sub = V_a_sub)
      }
    }
    A_post <- matrix(0, nrow = K, ncol = K)
    A_post[upper.tri(A)] <- a_sample
    A <- t(A_post)
    diag(A) <- 1

    # Sample w
    u <-  A %*% (yt - B %*%xt) + gamma * mu.xi
    w <- matrix( mapply( GIGrvg::rgig, n = 1, lambda = -(nu+1)*0.5, chi = nu + (u^2)/ sigma2,
                         psi = gamma^2/sigma2 ), K, t_max )
    # w <- matrix( rinvgamma( n = K*t_max, shape = (nu+1)*0.5, rate = 0.5*( nu + (u^2)/sigma2 ) ), K, t_max )

    # # Sample w
    # inv_sigma2 = 1/sigma2
    # u <- A %*% (yt - B %*%xt)
    # for (i in c(1:t_max)){
    #   ran_stu = rt(K, df = 4)
    #   a_eq <- -0.5* (gamma^2 * inv_sigma2)
    #   b_eq <- -((0.5*1 + nu/2+1))
    #   c_eq =  0.5 * (u[,i]^2 * inv_sigma2) + nu/2
    #   w_mode <- (- b_eq - sqrt(b_eq^2-4*a_eq*c_eq))/(2*a_eq)
    #   mode_target <- log( w_mode )
    #   cov_target <- 1/( c_eq / w_mode - a_eq * w_mode) * 1.2 # scale by 1.2
    #
    #   log_w_temp = mode_target + sqrt(cov_target)*ran_stu
    #   w_temp = exp(log_w_temp)
    #
    #   num_mh =  mapply( dinvgamma, x = w_temp, shape = nu*0.5, rate = nu*0.5, log = T) +  # prior
    #     ( - 0.5 * 1 * log_w_temp - 0.5 * (u[,i] - gamma*w_temp)^2 * inv_sigma2 / w_temp) - # posterior
    #     dt( (log_w_temp - mode_target)/sqrt(cov_target), df = 4, log = T) + log_w_temp # proposal
    #   denum_mh = mapply( dinvgamma, x = w[,i], shape = nu*0.5, rate = nu*0.5, log = T) +  # prior
    #     ( - 0.5 * 1 * log(w[,i]) - 0.5 * (u[,i] - gamma*w[,i])^2 * inv_sigma2 / w[,i]) - # posterior
    #     dt( (log(w[,i]) - mode_target)/sqrt(cov_target), df = 4, log = T) + log(w[,i]) # proposal
    #   alpha = num_mh - denum_mh
    #   temp = log(runif(K))
    #   acre = alpha > temp
    #   w[,i] <- ifelse(acre, w_temp, w[,i])
    #   w_sqrt[,i] <- sqrt(w[,i])
    #   acount_w[i] <- acount_w[i] + mean(acre)
    # }
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
    mu.xi <- nu/( nu - 2 )

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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, gamma, nu, as.numeric(w))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2), " ", round(gamma,2), " \n")
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

  return(as.mcmc(t(mcmc)))
}

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
  Sigma <- solve(inits$A0) %*% diag(inits$sigma, nrow = K)
  Sigma2 <- Sigma %*% t(Sigma)
  Sigma2_inv <- solve(Sigma2)
  V_b_prior_inv <- solve(V_b_prior)

  nu <- inits$nu
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

  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + 1 + K + K * t_max + t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # Sample B
    xt_G <- xt / reprow(sqrt(w_sample), m)
    yt_G <- (yt - Gamma_T * w)  / w_sqrt
    V_b_post_inv <- V_b_prior_inv + kronecker(xt_G %*% t(xt_G),Sigma2_inv)
    b_post <- kronecker(xt_G, Sigma2_inv) %*% vec(yt_G)

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
    A_inv <- solve(A)
    for (i in c(1:t_max)){
      Sigma_T[ ((i-1)*K + 1):(i*K) , ] <- A_inv %*% diag(sigma^2, nrow = K) %*% t(A_inv) * w_sample[i]
    }
    Gamma_filter <- fatBVARS::carterkohn(yt - B %*% xt, w_sample %x% diag(K),
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
    u_ort <- A %*% ((yt - B %*% xt - Gamma_T * w)/ w_sqrt)    # change from Gaussian
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
    sigma_post_b <- sigma0_S0 + apply(u_ort^2, 1, sum)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)

    # Sample A0
    u_std <- (yt - B %*% xt - Gamma_T * w)/ w_sqrt # change from Gaussian
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
    u <- (yt - B %*%xt)
    ran_stu = rt(t_max, df = 4)
    a_eq <- as.numeric(-0.5* sapply(1:t_max, FUN = function(i) (t(Gamma_T[,i]) %*% Sigma2_inv %*% Gamma_T[,i]) ))
    b_eq <- -((0.5*K + nu/2+1))
    c_eq =  as.numeric(nu*0.5 + 0.5 * apply((A %*% u / sigma)^2, MARGIN = 2, FUN = sum))
    w_mode <- (- b_eq - sqrt(b_eq^2-4*a_eq*c_eq))/(2*a_eq)
    w_mode <- ifelse(abs(w_mode)<1e-5, - c_eq/b_eq, w_mode)
    mode_target <- log( w_mode )
    cov_target <- 1/( c_eq / w_mode - a_eq * w_mode)
    log_w_temp = mode_target + sqrt(cov_target)*ran_stu
    w_temp = exp(log_w_temp)

    num_mh =  dinvgamma(w_temp, shape = nu*0.5, rate = nu*0.5, log = T) +  # prior
      ( - 0.5 * K * log_w_temp - 0.5 * apply((A %*% (u - Gamma_T * reprow(w_temp,K)) / sigma)^2, MARGIN = 2, FUN = sum) / w_temp) - # posterior
      dt( (log_w_temp - mode_target)/sqrt(cov_target), df = 4, log = T) + log_w_temp # proposal

    denum_mh = dinvgamma(w_sample, shape = nu*0.5, rate = nu*0.5, log = T) +  # prior
      ( - 0.5 * K * log(w_sample) - 0.5 * apply((A %*% (u - Gamma_T * w) / sigma)^2, MARGIN = 2, FUN = sum) / w_sample) - # posterior
      dt( (log(w_sample) - mode_target)/sqrt(cov_target), df = 4, log = T) + log(w_sample) # proposal
    alpha = as.numeric(num_mh - denum_mh)
    temp = log(runif(t_max))
    acre = alpha > temp
    w_sample <- ifelse(acre, w_temp, w_sample)
    w <- reprow(w_sample, K)
    w_sqrt <- reprow(sqrt(w_sample), K)
    acount_w <- acount_w + acre

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

###########################################################################
#' @export
BVAR.OrthSkewNorm1.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)

  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "OrthSkewNorm", SV = FALSE)
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
  sigma2 <- sigma^2

  V_b_prior_inv <- solve(V_b_prior)

  # Multi degrees of freedom
  nu <- inits$nu
  logsigma_nu <- rep(0,K)
  acount_nu <- rep(0,K)
  acount_w <- rep(0, t_max)
  gamma <- inits$gamma
  D <- diag(gamma)

  # Init w as Gaussian
  w_sample <- rep(0, t_max)
  w <- reprow(w_sample, K)
  # Output
  mcmc <- matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K + K*t_max,
                 ncol = (samples - inits$burnin)%/% inits$thin)
  for (j in c(1:samples)){
    # Sample B
    b_post = rep(0, m*K)
    V_b_post_inv = V_b_prior_inv
    inv_A <- solve(A)
    for (i in c(1:t_max)){
      V_b_post_inv <- V_b_post_inv + kronecker(xt[,i] %*% t(xt[,i]), t(A) %*% diag(1/sigma2) %*% A )
      b_post <- b_post + kronecker(xt[,i], (t(A)%*% diag(1/sigma2) %*% A) %*% (yt[,i] - inv_A %*% diag(1/sqrt(sigma2)) %*% (gamma * w[,i])) )
    }

    V_b_post <- solve(V_b_post_inv)
    b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    B <- Vec_to_Mat(b_sample, K,p)

    # Sample gamma
    gamma_post = rep(0, K)
    V_gamma_post_inv = solve(V_gamma_prior)
    u <- diag(1/sqrt(sigma2)) %*% A %*% (yt - B %*%xt)

    V_gamma_post_inv <- solve(V_gamma_prior) + diag(apply(w, MARGIN = 1, FUN = sum))
    gamma_post <- apply(u * w, MARGIN = 1, FUN = sum)

    V_gamma_post <- solve(V_gamma_post_inv)
    gamma_post <- V_gamma_post %*% ( solve(V_gamma_prior) %*% gamma_prior + gamma_post)
    gamma <- as.numeric(gamma_post + t(chol(V_gamma_post)) %*% rnorm(K))
    D <- diag(gamma)

    # Sample sigma
    u_ort <- (A %*% (yt - B %*%xt) - D %*% w)/ w_sqrt    # change from Gaussian
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
    sigma_post_b <- sigma0_S0 + apply(u_ort^2, 1, sum)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)

    # Sample A0
    u_std <- (yt - B %*%xt) # change from Gaussian
    u_neg <- - u_std
    a_sample <- rep(0, K * (K - 1) /2)
    if (K > 1) {
      for (i in c(2:K)){
        id_end <- i*(i-1)/2
        id_start <- id_end - i + 2
        a_sub <- a_prior[id_start:id_end]
        V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
        a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = (u_std[i,] - w[i,] * gamma[i]) / sigma[i]/ w_sqrt[i,],
                                                     xsub = matrix(u_neg[1:(i-1),] / sigma[i] / reprow(w_sqrt[i,], i-1), nrow = i-1),
                                                     a_sub = a_sub,
                                                     V_a_sub = V_a_sub)
      }
    }
    A_post <- matrix(0, nrow = K, ncol = K)
    A_post[upper.tri(A)] <- a_sample
    A <- t(A_post)
    diag(A) <- 1

    # Sample w
    inv_sigma2 = 1/sigma2
    u <- A %*% (yt - B %*%xt)
    for (i in c(1:t_max)){
      ran_stu = rt(K, df = 4)
      a_eq <- -0.5* (gamma^2 * inv_sigma2)
      b_eq <- -((0.5*1 + nu/2+1))
      c_eq =  0.5 * (u[,i]^2 * inv_sigma2) + nu/2
      w_mode <- (- b_eq - sqrt(b_eq^2-4*a_eq*c_eq))/(2*a_eq)
      mode_target <- log( w_mode )
      cov_target <- 1/( c_eq / w_mode - a_eq * w_mode) * 1.2 # scale by 1.2

      log_w_temp = mode_target + sqrt(cov_target)*ran_stu
      w_temp = exp(log_w_temp)

      num_mh =  mapply( dinvgamma, x = w_temp, shape = nu*0.5, rate = nu*0.5, log = T) +  # prior
        ( - 0.5 * 1 * log_w_temp - 0.5 * (u[,i] - gamma*w_temp)^2 * inv_sigma2 / w_temp) - # posterior
        dt( (log_w_temp - mode_target)/sqrt(cov_target), df = 4, log = T) + log_w_temp # proposal
      denum_mh = mapply( dinvgamma, x = w[,i], shape = nu*0.5, rate = nu*0.5, log = T) +  # prior
        ( - 0.5 * 1 * log(w[,i]) - 0.5 * (u[,i] - gamma*w[,i])^2 * inv_sigma2 / w[,i]) - # posterior
        dt( (log(w[,i]) - mode_target)/sqrt(cov_target), df = 4, log = T) + log(w[,i]) # proposal
      alpha = num_mh - denum_mh
      temp = log(runif(K))
      acre = alpha > temp
      w[,i] <- ifelse(acre, w_temp, w[,i])
      w_sqrt[,i] <- sqrt(w[,i])
      acount_w[i] <- acount_w[i] + mean(acre)
    }
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
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, gamma, nu, as.numeric(w))
    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", logsigma_nu," ", min(acount_w)," ", max(acount_w)," ", mean(acount_w), " ", round(nu,2), " ", round(gamma,2), " \n")
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

  return(as.mcmc(t(mcmc)))
}

############################################################################################## library(Matrix)
library(Matrix)
library(magic)
library(tmvnsim)
#' @export
BVAR.OrthSkewNorm.novol <- function(y, K, p, y0 = NULL, prior = NULL, inits = NULL){
  # Init regressors in the right hand side
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)

  # Init prior and initial values
  m = K * p + 1
  if (is.null(prior)){
    prior <- get_prior(y, p, priorStyle = "Minnesota", dist = "OrthSkewNorm", SV = FALSE)
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
  Sigma <- solve(inits$A0) %*% diag(sigma, nrow = length(sigma))
  Sigma2 <- Sigma %*% t(Sigma)
  Sigma2_inv <- solve(Sigma2)
  V_b_prior_inv <- solve(V_b_prior)

  # Multi degrees of freedom
  gamma <- inits$gamma
  D <- diag(gamma)

  # Init z as Truncated Gaussian
  z <- matrix(abs(rnorm(K * t_max)), ncol = t_max, nrow = K)

  # Output
  mcmc <-  matrix(NA, nrow = m*K + 0.5*K*(K-1) + K + K + K*t_max,
                  ncol = (samples - inits$burnin)%/% inits$thin)

  for (j in c(1:samples)){

    # Sample B
    yt_G <- (yt - solve(A) %*% (D %*% z) )

    V_b_post_inv <- V_b_prior_inv + kronecker(xt %*% t(xt),Sigma2_inv)
    b_post <- kronecker(xt, Sigma2_inv) %*% vec(yt_G)

    V_b_post <- solve(V_b_post_inv)
    b_post <- V_b_post %*% ( solve(V_b_prior) %*% b_prior + b_post)
    b_sample <- b_post + t(chol(V_b_post)) %*% rnorm(m*K)
    B <- Vec_to_Mat(b_sample, K,p)

    # Sample sigma
    sigma2 <- rep(0,K)
    u_ort <- (A %*% (yt - B %*%xt) - D %*% z)    # change from Gaussian
    sigma_post_a <- sigma0_T0 + rep(t_max,K)
    sigma_post_b <- sigma0_S0 + apply(u_ort^2, 1, sum)
    for (i in c(1:K)){
      sigma2[i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
    }
    sigma <- sqrt(sigma2)

    # Sample A0
    u_std <- (yt - B %*%xt) # change from Gaussian
    u_neg <- - u_std
    a_sample <- rep(0, K * (K - 1) /2)
    if (K > 1) {
      for (i in c(2:K)){
        id_end <- i*(i-1)/2
        id_start <- id_end - i + 2
        a_sub <- a_prior[id_start:id_end]
        V_a_sub <- V_a_prior[id_start:id_end, id_start:id_end]
        a_sample[c(id_start:id_end)] <- sample_A_ele(ysub = (u_std[i,] - z[i,] * gamma[i]) / sigma[i],
                                                     xsub = matrix(u_neg[1:(i-1),] / sigma[i], nrow = i-1),
                                                     a_sub = a_sub,
                                                     V_a_sub = V_a_sub)
      }
    }
    A_post <- matrix(0, nrow = K, ncol = K)
    A_post[upper.tri(A)] <- a_sample
    A <- t(A_post)
    diag(A) <- 1
    invA <- solve(A)

    Sigma <- invA %*% diag(sigma, nrow = length(sigma))
    Sigma2 <- Sigma %*% t(Sigma)
    Sigma2_inv <- solve(Sigma2)
    Sigma2_inv <- (Sigma2_inv + t(Sigma2_inv))*0.5

    # Sample gamma
    V_gamma_post_inv = solve(V_gamma_prior) + diag(apply(z^2, MARGIN = 1, FUN = sum) / sigma^2)
    V_gamma_post <- solve(V_gamma_post_inv)

    gamma_post <- apply(z  * (A %*% (yt - B %*% xt)), MARGIN = 1, FUN = sum) / sigma^2

    gamma_post <- V_gamma_post %*% ( solve(V_gamma_prior) %*% gamma_prior + gamma_post)
    gamma <- as.numeric(gamma_post + t(chol(V_gamma_post)) %*% rnorm(K))
    D <- diag(gamma)

    # Sample Z
    lb <- rep(0, K)
    for (i in c(1:t_max)){
      V_Zt_inv <- 1 + gamma^2 / sigma^2
      V_Zt <- 1 / V_Zt_inv
      mu_Zt <- as.numeric(V_Zt * gamma  / sigma^2 * (A %*% (yt[,i] - B %*% xt[,i])) )
      #z[,i] <- mapply(FUN =  msm::rtnorm, n = 1, mean = mu_Zt, sd = sqrt(V_Zt), lower=0) # library(msm)
      z[,i] <- mapply(FUN =  TruncatedNormal::rtnorm, n = 1, mu = mu_Zt, sd = sqrt(V_Zt), lb=0, ub = Inf, method = "fast") # library(TruncatedNormal)
      #z[,i] <- TruncatedNormal::rtmvnorm(n = 1, mu = mu_Zt, sigma = V_Zt, lb = lb) # library(TruncatedNormal)
      # z[,i] <- tmvtnorm::rtmvnorm(1, mean = mu_Zt, sigma = V_Zt, lower=rep(0, length = K)) # library(tmvtnorm)
      # z[,i] <- tmvnsim(1, K, means = mu_Zt, sigma = V_Zt, lower=rep(0, length = K))$samp # library(tmvnsim)
    }


    if ((j > inits$burnin) & (j %% inits$thin == 0))
      mcmc[, (j - inits$burnin) %/% inits$thin] <- c(b_sample, a_sample, sigma, gamma, as.numeric(z))

    if (j %% 1000 == 0) {
      cat(" Iteration ", j, " ", gamma ,  " \n")
    }
  }
  nameA <- matrix(paste("a", reprow(c(1:K),K), repcol(c(1:K),K), sep = "_"), ncol = K)
  nameA <- nameA[upper.tri(nameA, diag = F)]
  row.names(mcmc) <- c( paste("B0",c(1:K), sep = ""),
                        sprintf("B%d_%d_%d",reprow(c(1:p),K*K), rep(repcol(c(1:K),K), p), rep(reprow(c(1:K),K), p)),
                        nameA,
                        paste("sigma",c(1:K), sep = ""),
                        paste("gamma",c(1:K), sep = ""),
                        sprintf("w_%d_%d", repcol(c(1:K),t_max), reprow(c(1:t_max),K)))
  return(as.mcmc(t(mcmc)))
}

#############################################################################################
