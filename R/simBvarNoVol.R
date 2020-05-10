#' Simulate data from Gaussian BVAR model without SV
#'
#' This function simulate data from a BVAR model without SV.
#' \deqn{y_t = B x_t + SQRT(w_t) A^(-1) Sigma eps_t}
#' @param dist The variable specifies the BVAR error distribution. It should be one of
#' c("Gaussian","Student","Skew.Student","Hyper.Student", "multiStudent","Skew.multiStudent","Hyper.multiStudent").
#' @param K The number of variables in BVAR model.
#' @param p The number of lags in BVAR model.
#' @param t_max The number of observations
#' @param b0 The coefficients of B matrix.
#' @param a0 The coefficients of A matrix.
#' @param h The log diag(Sigma) matrix.
#' @param nu The degree of freedom.
#' @param gamma The skewness from \eqn{[-gamma, gamma]}.
#' @param seednum The default seed.
#' @param burn_in The discarded observations.
#' @return A list of simulated data with its specification.
#' @export
#' @examples
#' \dontrun{
#' datagen <- sim.VAR.novol(dist="Gaussian")
#' y <- datagen$y
#' }
sim.VAR.novol <- function(dist, K = 5, p = 2, t_max = 1000,
                                   b0 = 0.6, a0 = 0.5, h = 0, nu = 6, gamma = 0.5,
                                   y0 = matrix(0, ncol = K, nrow = p),
                                   seednum = 0, burn_in = 0){
  if (dist == "Gaussian") datagen <- sim.VAR.Gaussian.novol(K, p, t_max, b0, a0, h, y0, seednum, burn_in)
  if (dist == "Student") datagen <- sim.VAR.Student.novol(K, p, t_max, b0, a0, h, y0, nu, seednum, burn_in)
  # if (dist == "Skew.Student") datagen <- sim.VAR.Skew.Student.novol(K, p, t_max, b0, a0, h, nu, gamma, seednum, burn_in)
  if (dist == "Hyper.Student") datagen <- sim.VAR.Hyper.Student.novol(K, p, t_max, b0, a0, h, y0, nu, gamma, seednum, burn_in)
  if (dist == "multiStudent") datagen <- sim.VAR.multiStudent.novol(K, p, t_max, b0, a0, h, y0, nu, seednum, burn_in)
  if (dist == "Hyper.multiStudent") datagen <- sim.VAR.Hyper.multiStudent.novol(K, p, t_max, b0, a0, h, y0, nu, gamma, seednum, burn_in)
  if (dist == "multiOrthStudent") datagen <- sim.VAR.multiOrthStudent.novol(K, p, t_max, b0, a0, h, y0, nu, seednum, burn_in)
  if (dist == "Hyper.multiOrthStudent") datagen <- sim.VAR.Hyper.multiOrthStudent.novol(K, p, t_max, b0, a0, h, y0, nu, gamma, seednum, burn_in)

  return(datagen)
}

#' @export
sim.VAR.Gaussian.novol <- function(K = 5, p = 2, t_max = 1000,
                                   b0 = 0.6, a0 = 0.5, h = 0,
                                   y0 = matrix(0, ncol = K, nrow = p),
                                   seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (length(b0) == 1) {
    for (i in c(1:p)){
      B0 <- cbind(B0, b0^i*diag(K))
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
  # No volatility

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)

  Sigma <- solve(A0) %*% diag(as.vector(exp(0.5*h)))
  Sigma2 <- Sigma %*% t(Sigma)

  y_var <- reprow(Sigma[lower.tri(Sigma, diag = T)], t_max) # Sigma[lower.tri(Sigma, diag = T)] <- y_var[1,]
  volatility <- reprow(diag(Sigma2), t_max)

  for (i in c(1:t_max)){
    xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    y_mean[i,] <- B0 %*% xt
    ysim <-  B0 %*% xt + Sigma %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 =B0, h = h,
       dist = "Gaussian", SV = FALSE)
}
#' @export
sim.VAR.Student.novol <- function(K = 5, p = 2, t_max = 1000,
                                  b0 = 0.6, a0 = 0.5, h = 0,
                                  y0 = matrix(0, ncol = K, nrow = p),
                                  nu = 6, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (length(b0) == 1) {
    for (i in c(1:p)){
      B0 <- cbind(B0, b0^i*diag(K))
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

  # No volatility
  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)

  Sigma <- solve(A0) %*% diag(as.vector(exp(0.5*h)))
  Sigma2 <- Sigma %*% t(Sigma)

  for (i in c(1:t_max)){
    xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    y_mean[i,] <- B0 %*% xt
    Sigma_t <- diag(w_sqrt_t[i], nrow = K) %*% Sigma
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma_t %*% t(Sigma_t))
    ysim <-  B0 %*% xt + (Sigma_t %*% eps[i,])
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 =B0, h = h,
       nu = nu, w = w_t[(burn_in+1):(burn_in+t_max)],
       dist = "Student", SV = FALSE)
}
#' @export
sim.VAR.Hyper.Student.novol <- function(K = 5, p = 2, t_max = 1000,
                                        b0 = 0.6, a0 = 0.5, h = 0,
                                        y0 = matrix(0, ncol = K, nrow = p),
                                        nu = 6, gamma = 0, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (length(b0) == 1) {
    for (i in c(1:p)){
      B0 <- cbind(B0, b0^i*diag(K))
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
  D <- diag(gamma)
  # Tail of student
  w_t <- rinvgamma(t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # No volatility
  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)

  Sigma <- solve(A0) %*% diag(as.vector(exp(0.5*h)))
  Sigma2 <- Sigma %*% t(Sigma)

  for (i in c(1:t_max)){
    xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    y_mean[i,] <- B0 %*% xt + gamma * w_t[i]
    Sigma_t <- diag(w_sqrt_t[i], nrow = K) %*% Sigma
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma_t %*% t(Sigma_t))
    ysim <-  B0 %*% xt  + gamma * w_t[i] + (Sigma_t %*% eps[i,])
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 =B0, h = h,
       nu = nu, gamma = gamma,
       w = w_t[(burn_in+1):(burn_in+t_max)],
       dist = "Hyper.Student", SV = FALSE)
}
#' @export
sim.VAR.multiStudent.novol <- function(K = 5, p = 2, t_max = 1000,
                                       b0 = 0.6, a0 = 0.5, h = 0,
                                       y0 = matrix(0, ncol = K, nrow = p),
                                       nu = 6, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (length(b0) == 1) {
    for (i in c(1:p)){
      B0 <- cbind(B0, b0^i*diag(K))
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

  # No volatility
  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)

  Sigma <- solve(A0) %*% diag(as.vector(exp(0.5*h)))
  Sigma2 <- Sigma %*% t(Sigma)

  for (i in c(1:t_max)){
    xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    y_mean[i,] <- B0 %*% xt
    Sigma_t <- diag(w_sqrt_t[i,]) %*% Sigma
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma_t %*% t(Sigma_t))
    ysim <-  B0 %*% xt + (Sigma_t %*% eps[i,])
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 =B0, h = h,
       nu = nu, w = w_t[(burn_in+1):(burn_in+t_max),],
       dist = "multiStudent", SV = FALSE)
}
#' @export
sim.VAR.Hyper.multiStudent.novol <- function(K = 5, p = 2, t_max = 1000,
                                             b0 = 0.6, a0 = 0.5, h = 0,
                                             y0 = matrix(0, ncol = K, nrow = p),
                                             nu = 6, gamma = 0, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (length(b0) == 1) {
    for (i in c(1:p)){
      B0 <- cbind(B0, b0^i*diag(K))
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
  D <- diag(gamma)
  # Tail of student
  if (length(nu) == 1) nu = rep(nu,K)
  w_t <- mapply(rinvgamma, n = t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # No volatility
  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)

  Sigma <- solve(A0) %*% diag(exp(0.5*h))
  Sigma2 <- Sigma %*% t(Sigma)

  for (i in c(1:t_max)){
    xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    y_mean[i,] <- B0 %*% xt + gamma * w_t[i,]
    Sigma_t <- diag(w_sqrt_t[i,]) %*% Sigma
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma_t %*% t(Sigma_t))
    ysim <-  B0 %*% xt  + gamma * w_t[i,] + (Sigma_t %*% eps[i,])
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 =B0, h = h,
       nu = nu, gamma = gamma,
       w = w_t[(burn_in+1):(burn_in+t_max),],
       dist = "Hyper.multiStudent", SV = FALSE)
}
#' @export
sim.VAR.multiOrthStudent.novol <- function(K = 5, p = 2, t_max = 1000,
                                           b0 = 0.6, a0 = 0.5, h = 0,
                                           y0 = matrix(0, ncol = K, nrow = p),
                                           nu = 6, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (length(b0) == 1) {
    for (i in c(1:p)){
      B0 <- cbind(B0, b0^i*diag(K))
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

  # No volatility
  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)

  inv_A0 <- solve(A0)

  for (i in c(1:t_max)){
    xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    y_mean[i,] <- B0 %*% xt
    Sigma_t <-  inv_A0 %*% diag(w_sqrt_t[i,] * exp(0.5*h))
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma_t %*% t(Sigma_t))
    ysim <-  B0 %*% xt + (Sigma_t %*% eps[i,])
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 =B0, h = h,
       nu = nu, w = w_t[(burn_in+1):(burn_in+t_max),],
       dist = "multiOrthStudent", SV = FALSE)
}
#' @export
sim.VAR.Hyper.multiOrthStudent.novol <- function(K = 5, p = 2, t_max = 1000,
                                                 b0 = 0.6, a0 = 0.5, h = 0,
                                                 y0 = matrix(0, ncol = K, nrow = p),
                                                 nu = 6, gamma = 0, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (length(b0) == 1) {
    for (i in c(1:p)){
      B0 <- cbind(B0, b0^i*diag(K))
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
  D <- diag(gamma)
  # Tail of student
  if (length(nu) == 1) nu = rep(nu,K)
  w_t <- mapply(rinvgamma, n = t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # No volatility
  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)

  inv_A0 <- solve(A0)

  for (i in c(1:t_max)){
    xt <- rbind(1, vec( t(ystar[(p+i-1):i,])))
    y_mean[i,] <- B0 %*% xt + inv_A0 %*% (gamma * w_t[i,])
    Sigma_t <-  inv_A0 %*% diag(w_sqrt_t[i,] * exp(0.5*h))
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma_t %*% t(Sigma_t))
    ysim <-  as.numeric(y_mean[i,]) + (Sigma_t %*% eps[i,])
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 =B0, h = h,
       nu = nu, gamma = gamma,
       w = w_t[(burn_in+1):(burn_in+t_max),],
       dist = "Hyper.multiOrthStudent", SV = FALSE)
}
# #' @export
# sim.VAR.Skew.Student.novol <- function(K = 5, p = 2, t_max = 1000,
#                                        b0 = 0.6, a0 = 0.5, h = 0,
#                                        nu = 6, gamma = 0, seednum = 0, burn_in = 0){
#   t_max = t_max + burn_in
#   set.seed(seednum)
#   # Sample matrix coefficient B
#   B0 <- cbind(rep(0,K))
#   for (i in c(1:p)){
#     B0 <- cbind(B0, b0^i*diag(K))
#   }
#   # B0 <- vec(B0)
#
#   # Sample matrix corr A0
#   A0 <- matrix(rep(a0, K*K) , K )
#   diag(A0) <- 1
#   A0[upper.tri(A0)] <- 0
#
#
#   # Sample matrix variance h
#   if (length(h) == 1){
#     h <- rep(h,K)
#   }
#   # Skewness
#   if (length(gamma) == 1){
#     gamma <- seq(-gamma, gamma, length.out = K)
#   }
#   D <- diag(gamma)
#   z_t <- matrix(abs(rnorm(K*t_max)), ncol = K)
#   # Tail of student
#   w_t <- rinvgamma(t_max, shape = nu/2, rate = nu/2)
#   w_sqrt_t <- sqrt(w_t)
#
#   # No volatility
#   ystar <- matrix(0, ncol = K, nrow = p)
#
#   for (i in c(1:t_max)){
#     Sigma <- solve(A0) %*% diag(as.vector(exp(0.5*h)))
#     Sigma2 <- Sigma %*% t(Sigma)
#     xt <- rbind(1, vec( t(ystar)[,(p+i-1):i]))
#     ysim <-  B0 %*% xt  + D %*% z_t[i,] + w_sqrt_t[i] * (Sigma %*% rnorm(K))
#     ystar <- rbind(ystar, t(ysim))
#   }
#
#   t_max = t_max - burn_in
#   list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
#        K = K, p = p, t_max = t_max,
#        A0 = A0, B0 =B0, h = h,
#        nu = nu, gamma = gamma,
#        w = w_t[(burn_in+1):(burn_in+t_max)],
#        z = z_t,
#        dist = "Skew.Student", SV = FALSE)
# }
