#' Simulate data from Gaussian BVAR model with SV
#'
#' This function simulate data from a BVAR model with SV.
#' \deqn{y_t = B x_t + SQRT(w_t) A^(-1) Sigma_t eps_t}
#' @param dist The variable specifies the BVAR error distribution. It should be one of
#' c("Gaussian","Student","Skew.Student","Skew.Student", "MT","Skew.MT","MST").
#' @param K The number of variables in BVAR model.
#' @param p The number of lags in BVAR model.
#' @param t_max The number of observations
#' @param b0 The coefficients of B matrix.
#' @param a0 The coefficients of A matrix.
#' @param h The initial log diag(Sigma_t) matrix. sigma_h is set at diag(0.03, K)
#' @param nu The degree of freedom.
#' @param gamma The skewness from \eqn{[-gamma, gamma]}.
#' @param seednum The default seed.
#' @param burn_in The discarded observations.
#' @return A list of simulated data with its specification.
#' @export
#' @examples
#' \dontrun{
#' datagen <- sim.VAR.SV(dist="Gaussian")
#' y <- datagen$y
#' }
sim.VAR.SV <- function(dist, K = 5, p = 2, t_max = 1000,
                          b0 = 0.5, a0 = 0.5, h = 0, sigma_h = NULL,
                          y0 = matrix(0, ncol = K, nrow = p),
                          nu = 6, gamma = 0.5, sigma_G = NULL,
                          seednum = 0, burn_in = 0){
  if (!(dist %in% c("Gaussian","Student","Skew.Student",
                    "MT","MST",
                    "OT","OST",
                    "dynSkew.Student", "dynMST", "dynOST") ))
    stop("dist is not implemented.")

  if (dist == "Gaussian") datagen <- sim.VAR.Gaussian.SV(K, p, t_max, b0, a0, h, sigma_h, y0, seednum, burn_in)
  if (dist == "Student") datagen <- sim.VAR.Student.SV(K, p, t_max, b0, a0, h, sigma_h, y0, nu, seednum, burn_in)
  # if (dist == "Skew.Student") datagen <- sim.VAR.Skew.Student.SV(K, p, t_max, b0, a0, h, nu, gamma, seednum, burn_in)
  if (dist == "Skew.Student") datagen <- sim.VAR.Skew.Student.SV(K, p, t_max, b0, a0, h, sigma_h, y0, nu, gamma, seednum, burn_in)
  if (dist == "MT") datagen <- sim.VAR.MT.SV(K, p, t_max, b0, a0, h, sigma_h, y0, nu, seednum, burn_in)
  if (dist == "MST") datagen <- sim.VAR.MST.SV(K, p, t_max, b0, a0, h, sigma_h, y0, nu, gamma, seednum, burn_in)
  if (dist == "OT") datagen <- sim.VAR.OT.SV(K, p, t_max, b0, a0, h, sigma_h, y0, nu, seednum, burn_in)
  if (dist == "OST") datagen <- sim.VAR.OST.SV(K, p, t_max, b0, a0, h, sigma_h, y0, nu, gamma, seednum, burn_in)

  if (dist == "dynSkew.Student") datagen <- sim.VAR.dynSkew.Student.SV(K, p, t_max, b0, a0, h, sigma_h, y0, nu, gamma, sigma_G, seednum, burn_in)
  if (dist == "dynMST") datagen <- sim.VAR.dynMST.SV(K, p, t_max, b0, a0, h, sigma_h, y0, nu, gamma, sigma_G, seednum, burn_in)
  if (dist == "dynOST") datagen <- sim.VAR.dynOST.SV(K, p, t_max, b0, a0, h, sigma_h, y0, nu, gamma, sigma_G, seednum, burn_in)

  return(datagen)
}

#' @export
sim.VAR.Gaussian.SV <- function(K = 5, p = 2, t_max = 1000,
                                b0 = 0.5, a0 = 0.5, h = 0, sigma_h = NULL,
                                y0 = matrix(0, ncol = K, nrow = p),
                                seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }


  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  # No skew
  # No tail
  # Volatility volatility
  if (is.null(sigma_h)){
    Vh <- seq(3e-2, 3e-2, length.out = K)
  } else {
    Vh <-  sigma_h
  }

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  inv_A0 <- solve(A0)

  for (i in c(1:t_max)){
    h <- h +  Vh * rnorm(K)
    Sigma_t <-  inv_A0 %*% diag(exp(0.5*h), nrow = K)
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt
    ysim <-  B0 %*% xt + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K),
       Vh = Vh,
       dist = "Gaussian", SV = TRUE)
}

#' @export
sim.VAR.Student.SV <- function(K = 5, p = 2, t_max = 1000,
                               b0 = 0.5, a0 = 0.5, h = 0, sigma_h = NULL,
                               y0 = matrix(0, ncol = K, nrow = p),
                               nu = 6, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }


  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  # No skew
  # Tail of student
  w_t <- rinvgamma(t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # Volatility volatility
  if (is.null(sigma_h)){
    Vh <- seq(3e-2, 3e-2, length.out = K)
  } else {
    Vh <-  sigma_h
  }

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  inv_A0 <- solve(A0)

  for (i in c(1:t_max)){
    h <- h +  Vh * rnorm(K)
    Sigma_t <-  diag(w_sqrt_t[i], nrow = K) %*% inv_A0 %*% diag(exp(0.5*h), nrow = K)
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt
    ysim <-  B0 %*% xt + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K),
       nu = nu, w = w_t[(burn_in+1):(burn_in+t_max)],
       Vh = Vh,
       dist = "Student", SV = TRUE)
}

#' @export
sim.VAR.Skew.Student.SV <- function(K = 5, p = 2, t_max = 1000,
                                     b0 = 0.5, a0 = 0.5, h = 0, sigma_h = NULL,
                                     y0 = matrix(0, ncol = K, nrow = p),
                                     nu = 6, gamma = 0.5, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }


  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }

  # Skewness
  if(K > 1){
    if (length(gamma) == 1){
      gamma <- seq(-gamma, gamma, length.out = K)
    }
  }

  D <- diag(gamma, nrow = K)
  # Tail of student
  w_t <- rinvgamma(t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # Volatility volatility
  if (is.null(sigma_h)){
    Vh <- seq(3e-2, 3e-2, length.out = K)
  } else {
    Vh <-  sigma_h
  }

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  inv_A0 <- solve(A0)

  mu_xi <- nu / (nu - 2)
  for (i in c(1:t_max)){
    h <- h +  Vh * rnorm(K)
    Sigma_t <-  diag(w_sqrt_t[i], nrow = K) %*% inv_A0 %*% diag(exp(0.5*h), nrow = K)
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt + gamma * (w_t[i] - mu_xi)
    ysim <-  B0 %*% xt + gamma * (w_t[i] - mu_xi) + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K),
       nu = nu, gamma = gamma, w = w_t[(burn_in+1):(burn_in+t_max)],
       Vh = Vh,
       dist = "Skew.Student", SV = TRUE)
}
#' @export
sim.VAR.MT.SV <- function(K = 5, p = 2, t_max = 1000,
                                    b0 = 0.5, a0 = 0.5, h = 0, sigma_h = NULL,
                                    y0 = matrix(0, ncol = K, nrow = p),
                                    nu = 6, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }


  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }

  # No skew
  # Tail of student
  if (length(nu) == 1) nu = rep(nu,K)
  w_t <- mapply(rinvgamma, n = t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # Volatility volatility
  if (is.null(sigma_h)){
    Vh <- seq(3e-2, 3e-2, length.out = K)
  } else {
    Vh <-  sigma_h
  }

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  inv_A0 <- solve(A0)

  for (i in c(1:t_max)){
    h <- h +  Vh * rnorm(K)
    Sigma_t <-  diag(w_sqrt_t[i,]) %*% inv_A0 %*% diag(exp(0.5*h), nrow = K)
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt
    ysim <-  B0 %*% xt + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K),
       nu = nu, w = w_t[(burn_in+1):(burn_in+t_max),],
       Vh = Vh,
       dist = "MT", SV = TRUE)
}
#' @export
sim.VAR.MST.SV <- function(K = 5, p = 2, t_max = 1000,
                                          b0 = 0.5, a0 = 0.5, h = 0, sigma_h = NULL,
                                          y0 = matrix(0, ncol = K, nrow = p),
                                          nu = 6, gamma = 0.5, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }


  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }

  # Skewness
  if (length(gamma) == 1){
    gamma <- seq(-gamma, gamma, length.out = K)
  }
  D <- diag(gamma, nrow = K)
  # Tail of student
  if (length(nu) == 1) nu = rep(nu,K)
  mu_xi <- nu / (nu - 2)

  w_t <- mapply(rinvgamma, n = t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # Volatility volatility
  if (is.null(sigma_h)){
    Vh <- seq(3e-2, 3e-2, length.out = K)
  } else {
    Vh <-  sigma_h
  }

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  inv_A0 <- solve(A0)

  for (i in c(1:t_max)){
    h <- h +  Vh * rnorm(K)
    Sigma_t <-  diag(w_sqrt_t[i,]) %*% inv_A0 %*% diag(exp(0.5*h), nrow = K)
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt + gamma * (w_t[i,] - mu_xi)
    ysim <-  B0 %*% xt + gamma * (w_t[i,] - mu_xi) + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K), h = h,
       nu = nu, gamma = gamma,
       w = w_t[(burn_in+1):(burn_in+t_max),],
       Vh = Vh,
       dist = "MST", SV = TRUE)
}

#' @export
sim.VAR.OT.SV <- function(K = 5, p = 2, t_max = 1000,
                                        b0 = 0.5, a0 = 0.5, h = 0, sigma_h = NULL,
                                        y0 = matrix(0, ncol = K, nrow = p),
                                        nu = 6, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }


  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }

  # No skew
  # Tail of student
  if (length(nu) == 1) nu = rep(nu,K)
  w_t <- mapply(rinvgamma, n = t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # Volatility volatility
  if (is.null(sigma_h)){
    Vh <- seq(3e-2, 3e-2, length.out = K)
  } else {
    Vh <-  sigma_h
  }

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  inv_A0 <- solve(A0)

  for (i in c(1:t_max)){
    h <- h +  Vh * rnorm(K)
    Sigma_t <-  inv_A0 %*% diag(w_sqrt_t[i,] * exp(0.5*h))
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt
    ysim <-  B0 %*% xt + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K),
       nu = nu, w = w_t[(burn_in+1):(burn_in+t_max),],
       Vh = Vh,
       dist = "OT", SV = TRUE)
}
#' @export
sim.VAR.OST.SV <- function(K = 5, p = 2, t_max = 1000,
                                              b0 = 0.5, a0 = 0.5, h = 0, sigma_h = NULL,
                                              y0 = matrix(0, ncol = K, nrow = p),
                                              nu = 6, gamma = 0.5, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }


  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }

  # Skewness
  if (length(gamma) == 1){
    gamma <- seq(-gamma, gamma, length.out = K)
  }
  D <- diag(gamma, nrow = K)
  # Tail of student
  if (length(nu) == 1) nu = rep(nu,K)
  mu_xi <- nu / (nu - 2)

  w_t <- mapply(rinvgamma, n = t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # Volatility volatility
  if (is.null(sigma_h)){
    Vh <- seq(3e-2, 3e-2, length.out = K)
  } else {
    Vh <-  sigma_h
  }

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  inv_A0 <- solve(A0)

  for (i in c(1:t_max)){
    h <- h +  Vh * rnorm(K)
    Sigma_t <-  inv_A0 %*% diag(w_sqrt_t[i,] * exp(0.5*h))
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt + inv_A0 %*% (gamma * (w_t[i,] - mu_xi))
    ysim <-  as.numeric(y_mean[i,]) + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K), h = h,
       nu = nu, gamma = gamma,
       w = w_t[(burn_in+1):(burn_in+t_max),],
       Vh = Vh,
       dist = "OST", SV = TRUE)
}
#' @export
sim.VAR.dynSkew.Student.SV <- function(K = 5, p = 2, t_max = 1000,
                                       b0 = 0.5, a0 = 0.5, h = 0, sigma_h = NULL,
                                       y0 = matrix(0, ncol = K, nrow = p),
                                       nu = 6, gamma = 0.5, sigma_G = NULL, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }


  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }

  # Skewness
  if(K > 1){
    if (length(gamma) == 1){
      gamma <- seq(-gamma, gamma, length.out = K)
    }
  }

  D <- diag(gamma, nrow = K)
  # Tail of student
  w_t <- rinvgamma(t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # Volatility volatility
  if (is.null(sigma_h)){
    Vh <- seq(3e-2, 3e-2, length.out = K)
  } else {
    Vh <-  sigma_h
  }

  # Sigma_Gamma volatility
  if (is.null(sigma_G)){
    sigma_G <- seq(2e-2, 2e-2, length.out = K)
  }

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  Gamma_T <- matrix(NA, nrow = t_max, ncol = K)
  h_true <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  inv_A0 <- solve(A0)

  mu_xi <- nu / (nu - 2)
  for (i in c(1:t_max)){
    h <- h +  Vh * rnorm(K)
    h_true[i,] <- h
    Sigma_t <-  diag(w_sqrt_t[i], nrow = K) %*% inv_A0 %*% diag(exp(0.5*h), nrow = K)
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    gamma <- gamma +  sigma_G * rnorm(K)
    Gamma_T[i,] <- gamma

    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt + gamma * (w_t[i] - mu_xi)
    ysim <-  B0 %*% xt + gamma * (w_t[i] - mu_xi) + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       h_true = h_true,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K),
       nu = nu, gamma = gamma, w = w_t[(burn_in+1):(burn_in+t_max)],
       Vh = Vh,
       Gamma_T = Gamma_T,
       sigma_G = sigma_G,
       dist = "dynSkew.Student", SV = TRUE)
}

#' @export
sim.VAR.dynMST.SV <- function(K = 5, p = 2, t_max = 1000,
                                             b0 = 0.5, a0 = 0.5, h = 0, sigma_h = NULL,
                                             y0 = matrix(0, ncol = K, nrow = p),
                                             nu = 6, gamma = 0.5, sigma_G = NULL, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }


  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  lh0 <- h

  # Skewness
  if (length(gamma) == 1){
    gamma <- seq(-gamma, gamma, length.out = K)
  }
  D <- diag(gamma, nrow = K)
  # Tail of student
  if (length(nu) == 1) nu = rep(nu,K)
  w_t <- mapply(rinvgamma, n = t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # Volatility volatility
  if (is.null(sigma_h)){
    Vh <- seq(3e-2, 3e-2, length.out = K)
  } else {
    Vh <-  sigma_h
  }

  # Sigma_Gamma volatility
  if (is.null(sigma_G)){
    sigma_G <- seq(2e-2, 2e-2, length.out = K)
  }

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  Gamma_T <- matrix(NA, nrow = t_max, ncol = K)
  h_true <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  inv_A0 <- solve(A0)

  for (i in c(1:t_max)){
    h <- h +  Vh * rnorm(K)
    h_true[i,] <- h
    gamma <- gamma + sigma_G * rnorm(K)
    Gamma_T[i,] <- gamma
    Sigma_t <-  diag(w_sqrt_t[i,]) %*% inv_A0 %*% diag(exp(0.5*h), nrow = K)
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt + gamma * w_t[i,]
    ysim <-  B0 %*% xt + gamma * w_t[i,] + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K), h = h, lh0 = lh0,
       h_true = h_true,
       nu = nu, gamma = gamma,
       w = w_t[(burn_in+1):(burn_in+t_max),],
       Vh = Vh,
       Gamma_T = Gamma_T,
       sigma_G = sigma_G,
       dist = "dynMST", SV = TRUE)
}

#' @export
sim.VAR.dynOST.SV <- function(K = 5, p = 2, t_max = 1000,
                                                 b0 = 0.5, a0 = 0.5, h = 0, sigma_h = NULL,
                                                 y0 = matrix(0, ncol = K, nrow = p),
                                                 nu = 6, gamma = 0.5, sigma_G = NULL, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (p > 0){
    if (length(b0) == 1) {
      for (i in c(1:p)){
        B0 <- cbind(B0, b0^i*diag(K))
      }
    } else {
      B0 <- matrix(b0, nrow = K)
    }
  } else {
    B0 <- matrix(b0, nrow = K)
  }


  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  lh0 <- h

  # Skewness
  if (length(gamma) == 1){
    gamma <- seq(-gamma, gamma, length.out = K)
  }
  D <- diag(gamma, nrow = K)
  # Tail of student
  if (length(nu) == 1) nu = rep(nu,K)
  w_t <- mapply(rinvgamma, n = t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # Volatility volatility
  if (is.null(sigma_h)){
    Vh <- seq(3e-2, 3e-2, length.out = K)
  } else {
    Vh <-  sigma_h
  }

  # Sigma_Gamma volatility
  if (is.null(sigma_G)){
    sigma_G <- seq(2e-2, 2e-2, length.out = K)
  }

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  Gamma_T <- matrix(NA, nrow = t_max, ncol = K)
  h_true <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  inv_A0 <- solve(A0)

  for (i in c(1:t_max)){
    h <- h +  Vh * rnorm(K)
    h_true[i,] <- h
    Sigma_t <-  inv_A0 %*% diag(w_sqrt_t[i,] * exp(0.5*h))
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    gamma <- gamma +  sigma_G * rnorm(K)
    Gamma_T[i,] <- gamma

    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    if (p > 0){
      xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    } else {
      xt <- matrix(1,nrow = K)
    }
    y_mean[i,] <- B0 %*% xt + inv_A0 %*% (gamma * w_t[i,])
    ysim <-  as.numeric(y_mean[i,]) + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K), lh0 = lh0,
       h_true = h_true,
       nu = nu, gamma = gamma,
       w = w_t[(burn_in+1):(burn_in+t_max),],
       Vh = Vh,
       Gamma_T = Gamma_T,
       sigma_G = sigma_G,
       dist = "dynOST", SV = TRUE)
}

#' #' @export
#' sim.VAR.Skew.Student.SV <- function(K = 5, p = 2, t_max = 1000,
#'                                     b0 = 0.5, a0 = 0.5, h = 0,
#'                                     nu = 6, gamma = 0.5, seednum = 0, burn_in = 0){
#'   t_max = t_max + burn_in
#'   set.seed(seednum)
#'   # Sample matrix coefficient B
#'   B0 <- cbind(rep(0,K))
#'   for (i in c(1:p)){
#'     B0 <- cbind(B0, b0^i*diag(K))
#'   }
#'   # B0 <- vec(B0)
#'
#'   # Sample matrix corr A0
#'   A0 <- matrix(rep(a0, K*K) , K )
#'   diag(A0) <- 1
#'   A0[upper.tri(A0)] <- 0
#'
#'
#'   # Sample matrix variance h
#'   if (length(h) == 1){
#'     h <- rep(h,K)
#'   }
#'   Vh <- diag(seq(3e-2, 3e-2, length.out = K))
#'   # Skewness
#'   if (length(gamma) == 1){
#'     gamma <- seq(-gamma, gamma, length.out = K)
#'   }
#'   D <- diag(gamma, nrow = K)
#'   z_t <- matrix(abs(rnorm(K*t_max)), ncol = K)
#'   # Tail of student
#'   w_t <- rinvgamma(t_max, shape = nu/2, rate = nu/2)
#'   w_sqrt_t <- sqrt(w_t)
#'
#'   # Volatility volatility
#'   ystar <- matrix(0, ncol = K, nrow = p)
#'   volatility = NULL
#'
#'   for (i in c(1:t_max)){
#'     h <- h +  Vh %*% rnorm(K)
#'     Sigma <- solve(A0) %*% diag(as.vector(exp(0.5*h)))
#'     Sigma2 <- Sigma %*% t(Sigma)
#'     xt <- rbind(1, vec( t(ystar)[,(p+i-1):i]))
#'     ysim <-  B0 %*% xt  + D %*% z_t[i,] + w_sqrt_t[i] * (Sigma %*% rnorm(K))
#'     ystar <- rbind(ystar, t(ysim))
#'     volatility <- rbind(volatility, t(h))
#'   }
#'
#'   t_max = t_max - burn_in
#'   list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
#'        K = K, p = p, t_max = t_max,
#'        A0 = A0, B0 = matrix(B0, nrow = K),
#'        nu = nu, gamma = gamma,
#'        w = w_t[(burn_in+1):(burn_in+t_max)],
#'        z = z_t,
#'        dist = "Skew.Student", SV = TRUE)
#' }
