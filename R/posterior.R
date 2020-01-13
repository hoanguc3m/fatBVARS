
#' Posterior distribution of BVAR model
#'
#' This function returns a posterior distribution of BVAR-SV-fatTail model.
#' @param param The one posterior sample of the parameters obtained by fatBVARSVobj$mcmc
#' @param data The fatBVARSV object from command BVAR.
#' @return The unnormalized log posterior density of the fatBVARSVobj at param.
#' @export
#' @examples
#' \dontrun{
#' posterior <- get_posterior(Chain$mcmc[1,], Chain)
#' }
log_posterior <- function(param, data){

  K <- data$K
  p <- data$p

  y <- data$y
  y0 <- data$y0
  dist <- data$dist
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)

  prior <- data$prior
  SV <- prior$SV


  M <- K + p*(K^2) # nr of beta parameters
  numa <- 0.5*K*(K-1) # nr of VCV elements



  # B
  b0 = get_post(param, element = "B")
  B = matrix(b0, nrow = K)

  # A mean and var
  a0 <- get_post(param, element = "a")
  A = matrix(0, K, K)
  A[upper.tri(A, diag=FALSE)] = a0
  A <- t(A)
  diag(A) <- 1


  if (SV){
    # sigma_h
    sigma_h <- get_post(param, element = "sigma_h")
    h <- matrix(get_post(param, element = "h"), nrow = K)

    #Gaussian
    if (dist == "Gaussian"){
      log_prior <- mvnfast::dmvn(X = b0, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
        mvnfast::dmvn(X = a0, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
        sum(invgamma::dinvgamma(sigma_h^2, shape = 1* 0.5, rate = 0.0001 * 0.5, log = T)) +
        sum(mapply(dnorm, h[,c(2:t_max)], mean = h[,c(1:(t_max-1))], sd = sqrt(sigma_h), log = T)) + # apply ele in matrix using seq ele in array nu
        sum(mapply(dnorm, h[,1], mean = 0, sd = 1, log = T))

      log_likelihood <- sum(mvnfast::dmvn(X = (y - t(xt) %*% t(B)) %*% t(A) / exp(t(h*0.5)),
                                          mu = rep(0,K),
                                          sigma = diag(1, K), log = T, isChol = T)) - 0.5 * sum(h)
    }
    #Student
    if (dist == "Student"){

      # nu
      nu <- as.numeric(get_post(param, element = "nu"))
      w <- get_post(param, element = "w")

      log_prior <- mvnfast::dmvn(X = b0, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
        mvnfast::dmvn(X = a0, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
        sum(invgamma::dinvgamma(sigma_h^2, shape = 1* 0.5, rate = 0.0001 * 0.5, log = T)) +
        sum(mapply(dnorm, h[,c(2:t_max)], mean = h[,c(1:(t_max-1))], sd = sqrt(sigma_h), log = T)) + # apply ele in matrix using seq ele in array nu
        sum(mapply(dnorm, h[,1], mean = 0, sd = 1, log = T)) +
        sum(invgamma::dinvgamma(w, shape = nu*0.5, rate = nu*0.5, log = T)) +
        dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)


      log_likelihood <- sum(mvnfast::dmvn(X = ((y - t(xt) %*% t(B))/repcol(sqrt(w),K))  %*% t(A) / exp(t(h*0.5)),
                                          mu = rep(0,K),
                                          sigma = diag(1, K), log = T, isChol = T)) - 0.5 * sum(h) - 0.5 * K * sum(log(w)) # Jacobian
    }


    #Hyper.Student
    if (dist == "Hyper.Student"){

      # nu
      nu <- as.numeric(get_post(param, element = "nu"))
      w <- get_post(param, element = "w")
      gamma <- get_post(param, element = "gamma")

      log_prior <- mvnfast::dmvn(X = b0, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
        mvnfast::dmvn(X = a0, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
        sum(invgamma::dinvgamma(sigma_h^2, shape = 1* 0.5, rate = 0.0001 * 0.5, log = T)) +
        sum(mapply(dnorm, h[,c(2:t_max)], mean = h[,c(1:(t_max-1))], sd = sqrt(sigma_h), log = T)) + # apply ele in matrix using seq ele in array nu
        sum(mapply(dnorm, h[,1], mean = 0, sd = 1, log = T)) +
        sum(invgamma::dinvgamma(w, shape = nu*0.5, rate = nu*0.5, log = T)) +
        dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T) +
        mvnfast::dmvn(X = gamma, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)


      log_likelihood <- sum(mvnfast::dmvn(X = ((y - t(xt) %*% t(B) - w %*% t(gamma))/repcol(sqrt(w),K))  %*% t(A) / exp(t(h*0.5)),
                                          mu = rep(0,K),
                                          sigma = diag(1, K), log = T, isChol = T)) - 0.5 * sum(h) - 0.5 * K * sum(log(w)) # Jacobian
    }

    #multiStudent
    if (dist == "multiStudent"){

      # nu
      nu <- get_post(param, element = "nu")
      w <- matrix(get_post(param, element = "w"), nrow = K)


      log_prior <- mvnfast::dmvn(X = b0, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
        mvnfast::dmvn(X = a0, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
        sum(invgamma::dinvgamma(sigma_h^2, shape = 1* 0.5, rate = 0.0001 * 0.5, log = T)) +
        sum(mapply(dnorm, h[,c(2:t_max)], mean = h[,c(1:(t_max-1))], sd = sqrt(sigma_h), log = T)) + # apply ele in matrix using seq ele in array nu
        sum(mapply(dnorm, h[,1], mean = 0, sd = 1, log = T)) +
        sum(invgamma::dinvgamma(w, shape = nu*0.5, rate = nu*0.5, log = T)) +
        sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T))



      log_likelihood <- sum(mvnfast::dmvn(X = ((y - t(xt) %*% t(B))/sqrt(t(w)))  %*% t(A) / exp(t(h*0.5)),
                                          mu = rep(0,K),
                                          sigma = diag(1, K), log = T, isChol = T)) - 0.5 * sum(h) - 0.5 * sum(log(w)) # Jacobian
    }

    #Hyper.multiStudent
    if (dist == "Hyper.multiStudent"){

      # nu
      nu <- get_post(param, element = "nu")
      w <- matrix(get_post(param, element = "w"), nrow = K)
      gamma <- get_post(param, element = "gamma")

      log_prior <- mvnfast::dmvn(X = b0, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
        mvnfast::dmvn(X = a0, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
        sum(invgamma::dinvgamma(sigma_h^2, shape = 1* 0.5, rate = 0.0001 * 0.5, log = T)) +
        sum(mapply(dnorm, h[,c(2:t_max)], mean = h[,c(1:(t_max-1))], sd = sqrt(sigma_h), log = T)) + # apply ele in matrix using seq ele in array nu
        sum(mapply(dnorm, h[,1], mean = 0, sd = 1, log = T)) +
        sum(invgamma::dinvgamma(w, shape = nu*0.5, rate = nu*0.5, log = T)) +
        sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
        mvnfast::dmvn(X = gamma, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)


      log_likelihood <- sum(mvnfast::dmvn(X = ((y - t(xt) %*% t(B) - t(w * gamma))/sqrt(t(w))) %*% t(A) / exp(t(h*0.5)),
                                          mu = rep(0,K),
                                          sigma = diag(1, K), log = T, isChol = T)) - 0.5 * sum(h) - 0.5 * sum(log(w)) # Jacobian
    }
  } else {
    {
      # sigma
      sigma <- get_post(param, element = "sigma")

      #Gaussian
      if (dist == "Gaussian"){
        log_prior <- mvnfast::dmvn(X = b0, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
          mvnfast::dmvn(X = a0, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
          sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5, rate = prior$sigma_S0 * 0.5, log = T))
        log_likelihood <- sum(mvnfast::dmvn(X = (y - t(xt) %*% t(B)),
                                            mu = rep(0,K),
                                            sigma = solve(A) %*% diag(sigma), log = T, isChol = T))
      }

      #Student
      if (dist == "Student"){

        # nu
        nu <- as.numeric(get_post(param, element = "nu"))
        w <- get_post(param, element = "w")

        log_prior <- mvnfast::dmvn(X = b0, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
          mvnfast::dmvn(X = a0, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
          sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5, rate = prior$sigma_S0 * 0.5, log = T)) +
          sum(invgamma::dinvgamma(w, shape = nu*0.5, rate = nu*0.5, log = T)) +
          dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)


        log_likelihood <- sum(mvnfast::dmvn(X = (y - t(xt) %*% t(B))/repcol(sqrt(w),K),
                                            mu = rep(0,K),
                                            sigma = solve(A) %*% diag(sigma), log = T, isChol = T)) -
          0.5 * K * sum(log(w)) # Jacobian
      }

      #Hyper.Student
      if (dist == "Hyper.Student"){

        # nu
        nu <- as.numeric(get_post(param, element = "nu"))
        w <- get_post(param, element = "w")
        gamma <- get_post(param, element = "gamma")

        log_prior <- mvnfast::dmvn(X = b0, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
          mvnfast::dmvn(X = a0, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
          sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5, rate = prior$sigma_S0 * 0.5, log = T)) +
          sum(invgamma::dinvgamma(w, shape = nu*0.5, rate = nu*0.5, log = T)) +
          dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T) +
          mvnfast::dmvn(X = gamma, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
        # D = diag(gamma)
        log_likelihood <- sum(mvnfast::dmvn(X = (y - t(xt) %*% t(B) - w %*% t(gamma))/repcol(sqrt(w),K),
                                            mu = rep(0,K),
                                            sigma = solve(A) %*% diag(sigma), log = T, isChol = T)) -
          0.5 * K * sum(log(w)) # Jacobian
      }


      #multiStudent
      if (dist == "multiStudent"){

        # nu
        nu <- get_post(param, element = "nu")
        w <- matrix(get_post(param, element = "w"), nrow = K)

        log_prior <- mvnfast::dmvn(X = b0, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
          mvnfast::dmvn(X = a0, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
          sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5, rate = prior$sigma_S0 * 0.5, log = T)) +
          sum(mapply(invgamma::dinvgamma, t(w), shape = nu*0.5, rate = nu*0.5, log = T)) + # apply ele in matrix using seq ele in array nu
          sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T))

        log_likelihood <- sum(mvnfast::dmvn(X = (y - t(xt) %*% t(B))/sqrt(t(w)),
                                            mu = rep(0,K),
                                            sigma = solve(A) %*% diag(sigma), log = T, isChol = T)) -
          0.5 * sum(log(w)) # Jacobian
      }

      #Hyper.multiStudent
      if (dist == "Hyper.multiStudent"){

        # nu
        nu <- get_post(param, element = "nu")
        w <- matrix(get_post(param, element = "w"), nrow = K)
        gamma <- get_post(param, element = "gamma")

        log_prior <- mvnfast::dmvn(X = b0, mu = prior$b_prior, sigma = prior$V_b_prior, log = T) +
          mvnfast::dmvn(X = a0, mu = prior$a_prior, sigma = prior$V_a_prior, log = T) +
          sum(invgamma::dinvgamma(sigma^2, shape = prior$sigma_T0 * 0.5, rate = prior$sigma_S0 * 0.5, log = T)) +
          sum(mapply(invgamma::dinvgamma, t(w), shape = nu*0.5, rate = nu*0.5, log = T)) + # apply ele in matrix using seq ele in array nu
          sum(dgamma(nu, shape = prior$nu_gam_a, rate = prior$nu_gam_b, log = T)) +
          mvnfast::dmvn(X = gamma, mu = prior$gamma_prior, sigma = prior$V_gamma_prior, log = T)
        # D = diag(gamma)
        log_likelihood <- sum(mvnfast::dmvn(X = (y - t(xt) %*% t(B) - t(w * gamma))/sqrt(t(w)),
                                            mu = rep(0,K),
                                            sigma = solve(A) %*% diag(sigma), log = T, isChol = T)) -
          0.5 * sum(log(w)) # Jacobian
      }
    }
  }
  return(log_prior + log_likelihood)
}


#' Log Likelihood of BVAR model
#'
#' This function returns a log likelihood of BVAR-SV-fatTail model.
#' @param param The one sample of the parameters obtained by fatBVARSVobj$mcmc
#' @param data The fatBVARSV object from command BVAR.
#' @return The log likelihood density of the fatBVARSVobj at param.
#' @export
#' @examples
#' \dontrun{
#' eval_ll <- log_likelihood(Chain$mcmc[1,], Chain)
#' }
log_likelihood <- function(param, data, aggregate = TRUE){

  K <- data$K
  p <- data$p

  y <- data$y
  y0 <- data$y0
  dist <- data$dist
  t_max <- nrow(y)
  yt = t(y)
  xt <- makeRegressor(y, y0, t_max, K, p)

  prior <- data$prior
  SV <- prior$SV


  M <- K + p*(K^2) # nr of beta parameters
  numa <- 0.5*K*(K-1) # nr of VCV elements



  # B
  b0 = get_post(param, element = "B")
  B = matrix(b0, nrow = K)

  # A mean and var
  a0 <- get_post(param, element = "a")
  A = matrix(0, K, K)
  A[upper.tri(A, diag=FALSE)] = a0
  A <- t(A)
  diag(A) <- 1


  if (SV){
    # sigma_h
    sigma_h <- get_post(param, element = "sigma_h")
    h <- matrix(get_post(param, element = "h"), nrow = K)

    #Gaussian
    if (dist == "Gaussian"){
      log_likelihood <- mvnfast::dmvn(X = (y - t(xt) %*% t(B)) %*% t(A) / exp(t(h*0.5)),
                                          mu = rep(0,K),
                                          sigma = diag(1, K), log = T, isChol = T) - 0.5 * apply(h, MARGIN = 2, sum)

    }
    #Student
    if (dist == "Student"){

      # nu
      nu <- as.numeric(get_post(param, element = "nu"))
      w <- as.numeric(get_post(param, element = "w"))

      log_likelihood <- mvnfast::dmvn(X = ((y - t(xt) %*% t(B))/repcol(sqrt(w),K))  %*% t(A) / exp(t(h*0.5)),
                                          mu = rep(0,K),
                                          sigma = diag(1, K), log = T, isChol = T) - 0.5 * apply(h, MARGIN = 2, sum) - 0.5 * K * log(w) # Jacobian
    }


    #Hyper.Student
    if (dist == "Hyper.Student"){

      # nu
      nu <- as.numeric(get_post(param, element = "nu"))
      w <- get_post(param, element = "w")
      gamma <- get_post(param, element = "gamma")

      log_likelihood <- mvnfast::dmvn(X = ((y - t(xt) %*% t(B) - w %*% t(gamma))/repcol(sqrt(w),K))  %*% t(A) / exp(t(h*0.5)),
                                          mu = rep(0,K),
                                          sigma = diag(1, K), log = T, isChol = T) - 0.5 * apply(h, MARGIN = 2, sum) - 0.5 * K * log(w) # Jacobian
    }

    #multiStudent
    if (dist == "multiStudent"){

      # nu
      nu <- get_post(param, element = "nu")
      w <- matrix(get_post(param, element = "w"), nrow = K)

      log_likelihood <- mvnfast::dmvn(X = ((y - t(xt) %*% t(B))/sqrt(t(w)))  %*% t(A) / exp(t(h*0.5)),
                                          mu = rep(0,K),
                                          sigma = diag(1, K), log = T, isChol = T) - 0.5 * apply(h, MARGIN = 2, sum) - 0.5 * apply(log(w), MARGIN = 2, sum) # Jacobian
    }

    #Hyper.multiStudent
    if (dist == "Hyper.multiStudent"){

      # nu
      nu <- get_post(param, element = "nu")
      w <- matrix(get_post(param, element = "w"), nrow = K)
      gamma <- get_post(param, element = "gamma")

      log_likelihood <- mvnfast::dmvn(X = ((y - t(xt) %*% t(B) - t(w * gamma))/sqrt(t(w))) %*% t(A) / exp(t(h*0.5)),
                                          mu = rep(0,K),
                                          sigma = diag(1, K), log = T, isChol = T) - 0.5 * apply(h, MARGIN = 2, sum) - 0.5 * apply(log(w), MARGIN = 2, sum) # Jacobian
    }
  } else {
    {
      # sigma
      sigma <- get_post(param, element = "sigma")

      #Gaussian
      if (dist == "Gaussian"){
        log_likelihood <- mvnfast::dmvn(X = (y - t(xt) %*% t(B)),
                                            mu = rep(0,K),
                                            sigma = solve(A) %*% diag(sigma), log = T, isChol = T)
      }

      #Student
      if (dist == "Student"){

        # nu
        nu <- as.numeric(get_post(param, element = "nu"))
        w <- get_post(param, element = "w")

        log_likelihood <- mvnfast::dmvn(X = (y - t(xt) %*% t(B))/repcol(sqrt(w),K),
                                            mu = rep(0,K),
                                            sigma = solve(A) %*% diag(sigma), log = T, isChol = T) -
          0.5 * K * log(w) # Jacobian
      }

      #Hyper.Student
      if (dist == "Hyper.Student"){

        # nu
        nu <- as.numeric(get_post(param, element = "nu"))
        w <- get_post(param, element = "w")
        gamma <- get_post(param, element = "gamma")

        log_likelihood <- mvnfast::dmvn(X = (y - t(xt) %*% t(B) - w %*% t(gamma))/repcol(sqrt(w),K),
                                            mu = rep(0,K),
                                            sigma = solve(A) %*% diag(sigma), log = T, isChol = T) -
          0.5 * K * log(w) # Jacobian
      }


      #multiStudent
      if (dist == "multiStudent"){

        # nu
        nu <- get_post(param, element = "nu")
        w <- matrix(get_post(param, element = "w"), nrow = K)

        log_likelihood <- mvnfast::dmvn(X = (y - t(xt) %*% t(B))/sqrt(t(w)),
                                            mu = rep(0,K),
                                            sigma = solve(A) %*% diag(sigma), log = T, isChol = T) -
          0.5 * apply(log(w), MARGIN = 2, sum)  # Jacobian
      }

      #Hyper.multiStudent
      if (dist == "Hyper.multiStudent"){

        # nu
        nu <- get_post(param, element = "nu")
        w <- matrix(get_post(param, element = "w"), nrow = K)
        gamma <- get_post(param, element = "gamma")

        log_likelihood <- mvnfast::dmvn(X = (y - t(xt) %*% t(B) - t(w * gamma))/sqrt(t(w)),
                                            mu = rep(0,K),
                                            sigma = solve(A) %*% diag(sigma), log = T, isChol = T) -
          0.5 * apply(log(w), MARGIN = 2, sum) # Jacobian
      }
    }
  }
  if (aggregate) log_likelihood <- sum(log_likelihood)
  return(log_likelihood)
}


#' WAIC of BVAR model
#'
#' This function returns a WAIC of BVAR-SV-fatTail model.
#' @param data The fatBVARSV object from command BVAR.
#' @return The WAIC of the fatBVARSVobj at param.
#' @export
#' @examples
#' \dontrun{
#' WAIC_chain <- WAIC(Chain)
#' }
WAIC <- function(Chain){
  ndraws <- nrow(Chain$mcmc)
  t_max <- nrow(Chain$y)

  log_ll <- matrix(NA, nrow = ndraws, ncol = t_max)
  for (j in c(1:ndraws)){
    log_ll[j,] <- log_likelihood(param = Chain$mcmc[j,], data = Chain, aggregate = F)
    if(anyNA(log_ll[j,])) { cat("Error ", j, " \t")}
  }
  lpd = sum(log(apply(exp(log_ll), MARGIN = 2, mean))) # log predicted density
  p_waic = sum(apply(log_ll, MARGIN = 2, FUN = var)) # Effective parameters
  waic_deviance <- -2 * ( lpd - p_waic )
  return(list(waic = waic_deviance,
              lpd = lpd,
              p_waic = p_waic))
}


