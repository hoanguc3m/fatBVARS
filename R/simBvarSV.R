#' Simulate data from Gaussian BVAR model with SV
#'
#' This function simulate data from a BVAR model with SV.
#' \deqn{y_t = B x_t + SQRT(w_t) A^(-1) Sigma_t eps_t}
#' @param dist The variable specifies the BVAR error distribution. It should be one of
#' c("Gaussian","Student","Skew.Student","Hyper.Student", "multiStudent","Skew.multiStudent","Hyper.multiStudent").
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
                          b0 = 0.6, a0 = 0.5, h = 0, nu = 6, gamma = 0.5,
                          seednum = 0, burn_in = 0){
  if (dist == "Gaussian") datagen <- sim.VAR.Gaussian.SV(K, p, t_max, b0, a0, h, seednum, burn_in)
  if (dist == "Student") datagen <- sim.VAR.Student.SV(K, p, t_max, b0, a0, h, nu, seednum, burn_in)
  if (dist == "Skew.Student") datagen <- sim.VAR.Skew.Student.SV(K, p, t_max, b0, a0, h, nu, gamma, seednum, burn_in)
  if (dist == "Hyper.Student") datagen <- sim.VAR.Hyper.Student.SV(K, p, t_max, b0, a0, h, nu, gamma, seednum, burn_in)
  if (dist == "multiStudent") datagen <- sim.VAR.multiStudent.SV(K, p, t_max, b0, a0, h, nu, seednum, burn_in)
  if (dist == "Hyper.multiStudent") datagen <- sim.VAR.Hyper.multiStudent.SV(K, p, t_max, b0, a0, h, nu, gamma, seednum, burn_in)

  return(datagen)
}

#' @export
sim.VAR.Gaussian.SV <- function(K = 5, p = 2, t_max = 1000,
                                b0 = 0.6, a0 = 0.5, h = 0,
                                seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  for (i in c(1:p)){
    B0 <- cbind(B0, b0^i*diag(K))
  }
  #B0 <- vec(B0)

  # Sample matrix corr A0
  A0 <- matrix(rep(a0, K*K) , K )
  diag(A0) <- 1
  A0[upper.tri(A0)] <- 0

  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  Vh <- diag(seq(3e-2, 3e-2, length.out = K))
  # No skew
  # No tail
  # Volatility logvol

  ystar <- matrix(0, ncol = K, nrow = p)
  eps <- matrix(rnorm(t_max*K), ncol = K)
  logvol = NULL
  for (i in c(1:t_max)){
    h <- h +  Vh %*% rnorm(K)
    Sigma <- solve(A0) %*% diag(as.vector(exp(0.5*h)))
    Sigma2 <- Sigma %*% t(Sigma)
    xt <- rbind(1, vec( t(ystar)[,(p+i-1):i]))
    ysim <-  B0 %*% xt  + Sigma %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
    logvol <- rbind(logvol, t(h))
  }

  t_max = t_max - burn_in
  list(y = ystar[(p+burn_in+1):(p+burn_in+t_max),],
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K), Sigma = Sigma,
       logvol = logvol, Vh = Vh,
       dist = "Gaussian", SV = TRUE)
}

#' @export
sim.VAR.Student.SV <- function(K = 5, p = 2, t_max = 1000,
                               b0 = 0.6, a0 = 0.5, h = 0,
                               nu = 6, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  for (i in c(1:p)){
    B0 <- cbind(B0, b0^i*diag(K))
  }
  # B0 <- vec(B0)

  # Sample matrix corr A0
  A0 <- matrix(rep(a0, K*K) , K )
  diag(A0) <- 1
  A0[upper.tri(A0)] <- 0

  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  Vh <- diag(seq(3e-2, 3e-2, length.out = K))
  # No skew
  # Tail of student
  w_t <- rinvgamma(t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # Volatility logvol
  ystar <- matrix(0, ncol = K, nrow = p)
  logvol = NULL

  for (i in c(1:t_max)){
    h <- h +  Vh %*% rnorm(K)
    Sigma <- solve(A0) %*% diag(as.vector(exp(0.5*h)))
    Sigma2 <- Sigma %*% t(Sigma)
    xt <- rbind(1, vec( t(ystar)[,(p+i-1):i]))
    ysim <-  B0 %*% xt + w_sqrt_t[i] * (Sigma %*% rnorm(K))
    ystar <- rbind(ystar, t(ysim))
    logvol <- rbind(logvol, t(h))
  }

  t_max = t_max - burn_in
  list(y = ystar[(p+burn_in+1):(p+burn_in+t_max),],
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K), Sigma = Sigma,
       nu = nu, w = w_t[(burn_in+1):(burn_in+t_max)],
       logvol = logvol, Vh = Vh,
       dist = "Student", SV = TRUE)
}

#' @export
sim.VAR.Skew.Student.SV <- function(K = 5, p = 2, t_max = 1000,
                                    b0 = 0.6, a0 = 0.5, h = 0,
                                    nu = 6, gamma = 0.5, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  for (i in c(1:p)){
    B0 <- cbind(B0, b0^i*diag(K))
  }
  # B0 <- vec(B0)

  # Sample matrix corr A0
  A0 <- matrix(rep(a0, K*K) , K )
  diag(A0) <- 1
  A0[upper.tri(A0)] <- 0


  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  Vh <- diag(seq(3e-2, 3e-2, length.out = K))
  # Skewness
  if (length(gamma) == 1){
    gamma <- seq(-gamma, gamma, length.out = K)
  }
  D <- diag(gamma)
  z_t <- matrix(abs(rnorm(K*t_max)), ncol = K)
  # Tail of student
  w_t <- rinvgamma(t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # Volatility logvol
  ystar <- matrix(0, ncol = K, nrow = p)
  logvol = NULL

  for (i in c(1:t_max)){
    h <- h +  Vh %*% rnorm(K)
    Sigma <- solve(A0) %*% diag(as.vector(exp(0.5*h)))
    Sigma2 <- Sigma %*% t(Sigma)
    xt <- rbind(1, vec( t(ystar)[,(p+i-1):i]))
    ysim <-  B0 %*% xt  + D %*% z_t[i,] + w_sqrt_t[i] * (Sigma %*% rnorm(K))
    ystar <- rbind(ystar, t(ysim))
    logvol <- rbind(logvol, t(h))
  }

  t_max = t_max - burn_in
  list(y = ystar[(p+burn_in+1):(p+burn_in+t_max),],
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K), Sigma = Sigma,
       nu = nu, gamma = gamma,
       w = w_t[(burn_in+1):(burn_in+t_max)],
       z = z_t,
       dist = "Skew.Student", SV = TRUE)
}

#' @export
sim.VAR.Hyper.Student.SV <- function(K = 5, p = 2, t_max = 1000,
                                     b0 = 0.6, a0 = 0.5, h = 0,
                                     nu = 6, gamma = 0.5, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  for (i in c(1:p)){
    B0 <- cbind(B0, b0^i*diag(K))
  }
  # B0 <- vec(B0)

  # Sample matrix corr A0
  A0 <- matrix(rep(a0, K*K) , K )
  diag(A0) <- 1
  A0[upper.tri(A0)] <- 0


  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  Vh <- diag(seq(3e-2, 3e-2, length.out = K))
  # Skewness
  if (length(gamma) == 1){
    gamma <- seq(-gamma, gamma, length.out = K)
  }
  D <- diag(gamma)
  # Tail of student
  w_t <- rinvgamma(t_max, shape = nu/2, rate = nu/2)
  w_sqrt_t <- sqrt(w_t)

  # Volatility logvol
  ystar <- matrix(0, ncol = K, nrow = p)
  logvol = NULL

  for (i in c(1:t_max)){
    h <- h +  Vh %*% rnorm(K)
    Sigma <- solve(A0) %*% diag(as.vector(exp(0.5*h)))
    Sigma2 <- Sigma %*% t(Sigma)
    xt <- rbind(1, vec( t(ystar)[,(p+i-1):i]))
    ysim <-  B0 %*% xt  + gamma * w_t[i] + w_sqrt_t[i] * (Sigma %*% rnorm(K))
    ystar <- rbind(ystar, t(ysim))
    logvol <- rbind(logvol, t(h))
  }

  t_max = t_max - burn_in
  list(y = ystar[(p+burn_in+1):(p+burn_in+t_max),],
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K), Sigma = Sigma,
       nu = nu, gamma = gamma,
       w = w_t[(burn_in+1):(burn_in+t_max)],
       dist = "Hyper.Student", SV = TRUE)
}
#' @export
sim.VAR.multiStudent.SV <- function(K = 5, p = 2, t_max = 1000,
                                    b0 = 0.6, a0 = 0.5, h = 0,
                                    nu = 6, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  for (i in c(1:p)){
    B0 <- cbind(B0, b0^i*diag(K))
  }
  # B0 <- vec(B0)

  # Sample matrix corr A0
  A0 <- matrix(rep(a0, K*K) , K )
  diag(A0) <- 1
  A0[upper.tri(A0)] <- 0

  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  Vh <- diag(seq(3e-2, 3e-2, length.out = K))
  # No skew
  # Tail of student
  w_t <- matrix(rinvgamma(t_max*K, shape = nu/2, rate = nu/2), ncol = K)
  w_sqrt_t <- sqrt(w_t)

  # Volatility logvol
  ystar <- matrix(0, ncol = K, nrow = p)
  logvol = NULL

  for (i in c(1:t_max)){
    h <- h +  Vh %*% rnorm(K)
    Sigma <- solve(A0) %*% diag(as.vector(exp(0.5*h)))
    Sigma2 <- Sigma %*% t(Sigma)
    xt <- rbind(1, vec( t(ystar)[,(p+i-1):i]))
    ysim <-  B0 %*% xt + diag(w_sqrt_t[i,]) %*% (Sigma %*% rnorm(K))
    ystar <- rbind(ystar, t(ysim))
    logvol <- rbind(logvol, t(h))
  }

  t_max = t_max - burn_in
  list(y = ystar[(p+burn_in+1):(p+burn_in+t_max),],
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K), Sigma = Sigma,
       nu = nu, w = w_t[(burn_in+1):(burn_in+t_max),],
       logvol = logvol, Vh = Vh,
       dist = "multiStudent", SV = TRUE)
}
#' @export
sim.VAR.Hyper.multiStudent.SV <- function(K = 5, p = 2, t_max = 1000,
                                          b0 = 0.6, a0 = 0.5, h = 0,
                                          nu = 6, gamma = 0.5, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- cbind(rep(0,K))
  for (i in c(1:p)){
    B0 <- cbind(B0, b0^i*diag(K))
  }
  # B0 <- vec(B0)

  # Sample matrix corr A0
  A0 <- matrix(rep(a0, K*K) , K )
  diag(A0) <- 1
  A0[upper.tri(A0)] <- 0


  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  Vh <- diag(seq(3e-2, 3e-2, length.out = K))
  # Skewness
  if (length(gamma) == 1){
    gamma <- seq(-gamma, gamma, length.out = K)
  }
  D <- diag(gamma)
  # Tail of student
  w_t <- matrix(rinvgamma(t_max*K, shape = nu/2, rate = nu/2), ncol = K)
  w_sqrt_t <- sqrt(w_t)

  # No volatility
  ystar <- matrix(0, ncol = K, nrow = p)
  logvol = NULL

  for (i in c(1:t_max)){
    h <- h +  Vh %*% rnorm(K)
    Sigma <- solve(A0) %*% diag(as.vector(exp(0.5*h)))
    Sigma2 <- Sigma %*% t(Sigma)
    xt <- rbind(1, vec( t(ystar)[,(p+i-1):i]))
    ysim <-  B0 %*% xt  + gamma * w_t[i,] + diag(w_sqrt_t[i,]) %*% (Sigma %*% rnorm(K))
    ystar <- rbind(ystar, t(ysim))
    logvol <- rbind(logvol, t(h))
  }

  t_max = t_max - burn_in
  list(y = ystar[(p+burn_in+1):(p+burn_in+t_max),],
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K), h = h, Sigma = Sigma,
       nu = nu, gamma = gamma,
       w = w_t[(burn_in+1):(burn_in+t_max),],
       logvol = logvol, Vh = Vh,
       dist = "Hyper.multiStudent", SV = TRUE)
}
