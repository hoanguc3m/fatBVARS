# Minnesota prior from bvars Package
# URL: https://github.com/joergrieger/bvars

#' @export
get_prior_minnesota <- function(y, p, intercept=TRUE, ...){
  # lambda1=0.2, lambda2=0.5, lambda3=1, lambda4=2
  arguments <- eval(substitute(alist(...)))
  lambda1 <- ifelse(is.null(arguments$lambda1), 0.2, arguments$lambda1)
  lambda2 <- ifelse(is.null(arguments$lambda2), 0.5, arguments$lambda2)
  lambda3 <- ifelse(is.null(arguments$lambda3), 1, arguments$lambda3)
  lambda4 <- ifelse(is.null(arguments$lambda4), 2, arguments$lambda4)


  t_max <- nrow(y)
  K   <- ncol(y)

  # Prior b0 mean
  constant = 0
  B0 <- NULL
  if(intercept == TRUE) {
    constant = 1
    B0 <- cbind(rep(0,K))
  }

  if (p > 0) B0 <- cbind(B0, diag(K))
  if (p > 1){
    for (i in c(2:p)){
      B0 <- cbind(B0, matrix(0,K,K))
    }
  }
  b0 <- as.numeric(B0)

  #
  # Prior for covariance matrix
  #

  sigmasq <- rep(0,K)

  if (p > 0){
    for(ii in 1:K){
      Ylagi         <- stats::embed(y[,ii],dimension = p + 1)[,-1]
      Yi            <- y[(p + 1):t_max,ii]
      arest         <- stats::lm(Yi~Ylagi)
      sigmasq[ii] <- summary(arest)$sigma^2 # the square root of the estimated variance of the random error sigma^2 = 1/(n-p) Sum(w[i] R[i]^2)
    }
  } else {
    sigmasq <- apply(y, MARGIN = 2, FUN = var)
  }



  #print(sigmasq)

  #
  # Covariance matrix for the prior
  #

  M <- K * p + constant
  Vi <- rep(0, K * M)

  #
  # Add Covariance coefficients for intercepts
  #

  # Covariance for the intercept
  if (intercept==TRUE){
    for(ii in 1:K){
      Vi[ii]    <- lambda3 * sigmasq[ii]
    }
  }

  if (p > 0){

  # without intercept
  for(jj in 1:p){ #loop over the jj-th lag

    for(kk in 1:K){ #kk-th variable

      for(ii in 1:K){ # loop over the ii-th equation

        indx <- (kk - 1) * K  + (jj - 1) * K * K + ii + constant * K

        if(ii==kk){
          Vi[indx] <- lambda1/(jj^lambda4)
        } else{
          Vi[indx] <- (lambda1 * lambda2)/(jj ^ lambda4) *(sigmasq[ii] / sigmasq[kk])
        }
      }
    }
  }

  }
  # Vec_to_Mat(Vi, K, p)

  V_b_prior = diag(Vi, nrow = K * M, ncol = K * M)

  pr <- list(type = "Minnesota",
             p = p,
             intercept = intercept,
             b0 = b0,
             V_b_prior = V_b_prior,
             sigma = sqrt(sigmasq))
  return(pr)
}

######################################################################################
get_prior_OLS <- function(y, p, tau = nrow(y) - p, scale_factor = 4 ){
  priordat = y[1:(tau+p),]
  K <- ncol(priordat)
  M <- K + p*(K^2) # nr of beta parameters
  numa <- 0.5*K*(K-1) # nr of VCV elements

  xt <- NULL
  t_max <- nrow(priordat)

  yt = t(y)
  for (i in c( p:(t_max-1))){
    xt <- cbind(xt, rbind(1, vec(yt[,i:(i-p+1)])))
  }

  #t_max <- t_max - p
  yt <- t(priordat[(p+1):(dim(priordat)[1]),])  # LHS variable
  tau <- dim(yt)[2] # New sample size

  # B0 mean and var
  B0 = yt %*% t(xt) %*% solve( xt %*% t(xt) )
  b0 = matrixcalc::vec(B0)
  H = t(xt) %*% solve( xt %*% t(xt) ) %*% xt
  SSE = yt %*% ( diag(rep(1, tau)) - H ) %*% t(yt)
  Sigma2 = SSE / (tau - K * p -1)
  V_b_prior = kronecker( solve(xt %*% t(xt)), Sigma2  ) * scale_factor # scale factor

  # A0 mean and var
  A0 <- t(chol(Sigma2))
  sigma <- diag(A0)

  pr <- list(type = "OLS",
             p = p,
             intercept = T,
             b0 = b0,
             V_b_prior = V_b_prior,
             sigma = sigma)

}

######################################################################################

#' Prior distribution of BVAR model
#'
#' This function returns a prior specification of BVAR-SV-fatTail model.
#' @param y The input data as matrix T x K where T is the number of observations and K is the number of variables
#' @param p The number of lags in BVAR model.
#' @param priorStyle The prior style in BVAR model should be c("Minnesota", "OLS")
#' @param dist The variable specifies the BVAR error distribution. It should be one of
#' c("Gaussian","Student", "Skew.Student", "MT","MST","orthoStudent").
#' @param SV The indicator if this BVAR model has a Stochastic volatility part.
#' @return A list of prior specifications. \eqn{b \sim N(b0, V_b_prior)}, \eqn{a \sim N(0, 1000I)}, \eqn{sigmaSq \sim IG(0.5*sigma_T0, 0.5*sigma_S0)},
#' \eqn{nu \sim G(a,b)}, \eqn{gamma \sim N(0, I)}.
#' @export
#' @examples
#' \dontrun{
#' prior <- get_prior(y, p = 2, dist = c("Student"), SV = FALSE)
#' }
get_prior <- function(y, p, priorStyle = c("Minnesota"),
                      dist = c("Gaussian"),
                      SV = FALSE, ...){
  if (!(dist %in% c("Gaussian","Student","Skew.Student",
                    "MT","MST",
                    "OT","OST",
                    "dynSkew.Student", "OrthSkewNorm", "dynOST") ))
    stop("dist is not implemented.")

  arguments <- eval(substitute(alist(...)))

  K <- ncol(y)
  M <- K + p*(K^2) # nr of beta parameters
  numa <- 0.5*K*(K-1) # nr of VCV elements

  # return b0, V_b_prior, sigma
  if (priorStyle == "Minnesota"){
    prior_sub <- get_prior_minnesota(y, p, intercept = TRUE, ...)
  }
  if (priorStyle == "OLS"){

  }


  # B mean and var
  b0 = prior_sub$b0
  V_b_prior = prior_sub$V_b_prior

  # A mean and var
  a0 <- rep(0, numa) # A \sim N(0,10I)
  a_Vprior <- ifelse(is.null(arguments$a_Vprior), 10, arguments$a_Vprior)
  V_a_prior <- diag(x = rep(a_Vprior,numa), nrow = numa, ncol = numa)

  sigma <- prior_sub$sigma

  #Gaussian
  prior_collect <- list(b_prior=b0, V_b_prior = V_b_prior,
                        a_prior = a0, V_a_prior = V_a_prior,
                        sigma = sigma,
                        sigma_T0 = 1, sigma_S0 = 1)
  #Student
  if (dist !="Gaussian"){
    prior_collect$nu_gam_a = 2
    prior_collect$nu_gam_b = 0.1
  }
  #Skew.Student
  if (dist =="Skew.Student" |
      dist =="MST" | dist =="OST"|
      dist =="dynSkew.Student" | dist =="OrthSkewNorm"| dist =="dynOST"){
    prior_collect$gamma_prior = rep(0, K)
    prior_collect$V_gamma_prior = diag(K)
  }
  prior_collect$t_max <- nrow(y)
  prior_collect$dist <- dist
  prior_collect$SV <- SV
  prior_collect$priorStyle <- priorStyle

  if (SV) {
    prior_collect$sigma_h <- 1;
  }
  return(prior_collect)
}

######################################################################################
#' Initial values from the prior of BVAR model
#'
#' This function returns initial values of BVAR-SV-fatTail model from its prior.
#' @param prior The prior distribution.
#' @return A list of initial values.
#' @export
#' @examples
#' \dontrun{
#' inits <- get_init(prior)
#' }
#'
get_init <- function(prior, samples = 1100, burnin = 100, thin = 1){
  dist = prior$dist
  if (!(dist %in% c("Gaussian","Student","Skew.Student",
                    "MT","MST",
                    "OT","OST",
                    "dynSkew.Student", "OrthSkewNorm", "dynOST") ))
    stop("dist is not implemented.")

  SV = prior$SV
  K = length(prior$sigma)
  p = (length(prior$b_prior) / K - 1) / K

  A = matrix(0, K, K)
  A[upper.tri(A, diag=FALSE)] = prior$a_prior + sqrt(diag(prior$V_a_prior)) * rnorm(K*(K-1)*0.5)
  A <- t(A)
  diag(A) <- 1

  B = matrix(prior$b_prior, nrow = K) + sqrt(diag(prior$V_b_prior)) * rnorm(K + K*K*p)
  sigma = prior$sigma + rnorm(K, 0,1)^2

  #Gaussian
  inits <- list(samples = samples,
                burnin = burnin,
                thin = thin,
                A0 = A,
                B0 = B,
                sigma = sigma
  )

  #Student
  if (dist !="Gaussian"){
    inits$nu = 4 + rgamma(1, shape = 2, rate = 0.1)
  }
  #Skew.Student
  if (dist =="Skew.Student" |
      dist =="MST" | dist =="OST"|
      dist =="dynSkew.Student" | dist =="OrthSkewNorm"| dist =="dynOST"){
    inits$gamma = rep(0.001,K) + 0.1 *rnorm(K)
  }
  #MT
  if (dist == "MT" | dist =="OT" |
      dist == "MST" | dist == "OST" |
      dist == "dynMST"| dist == "dynOST"){
    inits$nu = 4 + rgamma(K, shape = 2, rate = 0.1)
  }
  inits$nu <- ifelse(inits$nu > 100, 50, inits$nu)

  #Stochastic vol
  if(SV){
    h <- matrix(log(prior$sigma^2), ncol = prior$t_max, nrow = K) + 0.0001 * rnorm(prior$t_max * K)
    inits$h <- h
    inits$sigma_h <- 0.0001*diag(K)
  }
  return(inits)
}
