# importance_check: input importance weights


#' @export
importance_check <- function(sum_log){
  N <- length(sum_log)
  impt_w <- exp(sum_log - max(sum_log))
  # Threshold u
  prop_vec <- c(0.55, 0.5, 0.4, 0.3, 0.1, 0.01)

  test1_results <- rep(NA, length(prop_vec)) # Wald
  test2_results <- rep(NA, length(prop_vec)) # Score
  test3_results <- rep(NA, length(prop_vec)) # LR

  for (proportion in prop_vec){
    j <- (1:length(prop_vec))[proportion == prop_vec]
    u <- sort(impt_w)[(1 - proportion) * N]
    z_all <-ifelse(impt_w > u, impt_w - u, 0)
    z <- z_all[z_all >0]
    n <- length(z)


    # Wald test
    f_tau <- function(tau){
      xi <- sum(log(1+ tau * z)) / n
      return(n/ tau - (1+1/xi) * sum ( z/(1+tau*z) ))
    }

    #tau_init <- 1/ mean(z)
    tau <- uniroot(f = f_tau, interval = c(0.1/ mean(z),100/mean(z) ))$root
    #tau <- uniroot(f = f_tau, interval = c(1,1e6))$root
    # cat("tau = ", tau)

    xi <- sum(log(1+ tau * z)) / n
    beta <- xi/tau
    t_stat <- sqrt(n/3/beta^2) *(xi - 0.5)
    t_norm <- qnorm(0.95)
    test1_results[j] <- t_stat # (t_stat > t_norm) # Wald

    # Score test

    score_beta_restrict <- function(beta){
      return(n / beta * ( 3/ n * sum(z/(2 * beta + z)) - 1  ))
    }


    beta_restrict <- uniroot(f = score_beta_restrict, interval = c(.Machine$double.eps^0.5,5),
                             tol = .Machine$double.eps)$root
    # cat("beta = ", beta_restrict)
    #score_beta_restrict(beta_restrict)
    score_xi <- function(beta){
      4 * sum(log(1 + 0.5 * z/beta)) - 6 * sum(z / (2 * beta + z))
    }

    test2_results[j] <- score_xi(beta_restrict) / sqrt(2 * n) # ((score_xi(beta_restrict) / sqrt(2 * n)) > t_norm)

    # LR test

    pareto_fztheta <- function(xi, beta){
       return(- n * log(beta) - (1+1/xi) *sum(log (1+xi/beta*z)))
    }
    # pareto_fztheta <- function(par){
    #   xi = par[1]; beta = par[2]
    #   return( n * log(beta) + (1+1/xi) *sum(log (1+xi/beta*z))) # Negative likelihood
    # }
    # theta_res <- optim(par = c(1,beta_restrict), fn = pareto_fztheta, method = "L-BFGS-B", lower = c(0.5,0.0001))
    # cat(c(theta_res$par, beta_restrict) , "\t")
    # LR <- 2 * ( - pareto_fztheta( c(xi2,beta2)) + # Negative
    #               pareto_fztheta( c(0.5,beta_restrict) ) )

    if (xi > 0.5) {
      xi2 <- xi
      beta2 <- beta
    } else {
      xi2 <- 0.5
      beta2 <- beta_restrict
    }

    LR <- 2 * ( pareto_fztheta(xi = xi2, beta = beta2) -
                  pareto_fztheta(xi = 0.5, beta = beta_restrict) )

    test3_results[j] <- LR # (LR > 0.5*(qchisq(p = 0.95, df = 0) + qchisq(p = 0.95, df = 1)  )  )
  }

  return(list(Wald_test = test1_results,
         Score_test = test2_results,
         LR_test = test3_results,
         norm_critical = t_norm,
         chi_critical = 2.69 # 0.5*(qchisq(p = 0.95, df = 0) + qchisq(p = 0.95, df = 1))
         ))
}
