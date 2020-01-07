
#' Forecast future observations of BVAR model
#'
#' This function returns a forecast future observations of BVAR-SV-fatTail model.
#' @param Chain The fatBVARSV object from command BVAR.
#' @param t_pred The time prediction horizon.
#' @return The list of forecast future observations of BVAR model
#' @export
#' @examples
#' \dontrun{
#' forecast1 <- get_forecast(Chain1)
#' }
#'
get_forecast <- function(Chain, y0 = NULL, t_pred = 12, Nfsample = NULL){
  K <- Chain$K
  p <- Chain$p

  if( is.null(y0)) {
    # y <- Chain$y
    y0 <- tail(Chain$y, p) # Get last p obs
  }

  dist <- Chain$dist
  SV <- Chain$prior$SV

  mcmc <- Chain$mcmc
  if (is.null(Nfsample)){
    Nfsample <- nrow(mcmc)
    frange <- c(1:Nfsample)
  } else {
    frange <- sample(nrow(mcmc), Nfsample)
  }



  y_pred <- array(NA, dim = c(t_pred, K, Nfsample))
  ymean_pred <- array(NA, dim = c(t_pred, K, Nfsample))
  yvar_pred <- array(NA, dim = c(t_pred, K, Nfsample))
  vol_pred <- array(NA, dim = c(t_pred, K, Nfsample))

  for (i in frange){
    b0 <- get_post(mcmc[i,], element = "B")
    a0 <- get_post(mcmc[i,], element = "a")

    if (SV) {
        sigma_h <- get_post(mcmc[i,], element = "sigma_h")
        h <- as.numeric(tail(t(matrix(get_post(mcmc[i,], element = "h"), nrow = K)) ,1))
        if (dist == "Gaussian") {
          pred_tmp <- sim.VAR.Gaussian.SV(K = K, p = p, t_max = t_pred,
                                             b0 = b0, a0 = a0, h = h, sigma_h = sigma_h,
                                             y0 = y0,
                                             seednum = sample(1:100000,1), burn_in = 0)
        }

        if (dist == "Student") {
          nu <- get_post(mcmc[i,], element = "nu")
          pred_tmp <- sim.VAR.Student.SV(K = K, p = p, t_max = t_pred,
                                            b0 = b0, a0 = a0, h = h, sigma_h = sigma_h,
                                            y0 = y0, nu = nu,
                                            seednum = sample(1:100000,1), burn_in = 0)
        }
        if (dist == "Hyper.Student") {
          nu <- get_post(mcmc[i,], element = "nu")
          gamma <- get_post(mcmc[i,], element = "gamma")
          pred_tmp <- sim.VAR.Hyper.Student.SV(K = K, p = p, t_max = t_pred,
                                                  b0 = b0, a0 = a0, h = h, sigma_h = sigma_h,
                                                  y0 = y0, nu = nu, gamma = gamma,
                                                  seednum = sample(1:100000,1), burn_in = 0)
        }
        if (dist == "multiStudent") {
          nu <- get_post(mcmc[i,], element = "nu")
          pred_tmp <- sim.VAR.multiStudent.SV(K = K, p = p, t_max = t_pred,
                                                 b0 = b0, a0 = a0, h = h, sigma_h = sigma_h,
                                                 y0 = y0, nu = nu,
                                                 seednum = sample(1:100000,1), burn_in = 0)
        }
        if (dist == "Hyper.multiStudent") {
          nu <- get_post(mcmc[i,], element = "nu")
          gamma <- get_post(mcmc[i,], element = "gamma")
          pred_tmp <- sim.VAR.Hyper.multiStudent.SV(K = K, p = p, t_max = t_pred,
                                                       b0 = b0, a0 = a0, h = h, sigma_h = sigma_h,
                                                       y0 = y0, nu = nu, gamma = gamma,
                                                       seednum = sample(1:100000,1), burn_in = 0)
        }
      } else {
        sigma <- get_post(mcmc[i,], element = "sigma")
        h <- log(sigma)
        if (dist == "Gaussian") {
            pred_tmp <- sim.VAR.Gaussian.novol(K = K, p = p, t_max = t_pred,
                                               b0 = b0, a0 = a0, h = h,
                                               y0 = y0,
                                               seednum = sample(1:100000,1), burn_in = 0)
          }
        if (dist == "Student") {
          nu <- get_post(mcmc[i,], element = "nu")
          pred_tmp <- sim.VAR.Student.novol(K = K, p = p, t_max = t_pred,
                                             b0 = b0, a0 = a0, h = h,
                                             y0 = y0, nu = nu,
                                             seednum = sample(1:100000,1), burn_in = 0)
        }
        if (dist == "Hyper.Student") {
          nu <- get_post(mcmc[i,], element = "nu")
          gamma <- get_post(mcmc[i,], element = "gamma")
          pred_tmp <- sim.VAR.Hyper.Student.novol(K = K, p = p, t_max = t_pred,
                                            b0 = b0, a0 = a0, h = h,
                                            y0 = y0, nu = nu, gamma = gamma,
                                            seednum = sample(1:100000,1), burn_in = 0)
        }
        if (dist == "multiStudent") {
          nu <- get_post(mcmc[i,], element = "nu")
          pred_tmp <- sim.VAR.multiStudent.novol(K = K, p = p, t_max = t_pred,
                                            b0 = b0, a0 = a0, h = h,
                                            y0 = y0, nu = nu,
                                            seednum = sample(1:100000,1), burn_in = 0)
        }
        if (dist == "Hyper.multiStudent") {
          nu <- get_post(mcmc[i,], element = "nu")
          gamma <- get_post(mcmc[i,], element = "gamma")
          pred_tmp <- sim.VAR.Hyper.multiStudent.novol(K = K, p = p, t_max = t_pred,
                                                 b0 = b0, a0 = a0, h = h,
                                                 y0 = y0, nu = nu, gamma = gamma,
                                                 seednum = sample(1:100000,1), burn_in = 0)
        }
      }

    y_pred[,,i] <- pred_tmp$y
    ymean_pred[,,i] <- pred_tmp$y_mean
    yvar_pred[,,i] <- pred_tmp$y_var
    vol_pred[,,i] <- pred_tmp$logvol

  }
  return(list(y_pred = y_pred,
              ymean_pred = ymean_pred,
              yvar_pred = yvar_pred,
              vol_pred = vol_pred))
}

#' @export
forecast_density <- function(Chain, y_current = NULL, y_obs_future){
  K <- Chain$K
  p <- Chain$p

  if (! ncol(y_obs_future) == K) { stop("ncol(y_obs_future) != K") }

  t_pred = nrow(y_obs_future)
  predictive_samples <- get_forecast(Chain = Chain, y0 = y_current, t_pred = t_pred) # Nfsample

  log_pred<- matrix(NA, nrow = t_pred, ncol = K)
  emp_CDF<- matrix(NA, nrow = t_pred, ncol = K)
  log_biv_pred <- array(NA, dim = c(K*(K-1)*0.5,t_pred))

  for (i in c(1:K)){
    for (j in c(1:t_pred)){
      #predict_den <- stats::density(predictive_samples$y_pred[j,i,])
      predict_den <- KernSmooth::bkde(x = predictive_samples$y_pred[j,i,]) # fast and accurate
      y_min <- min(predict_den$x)
      y_max <- max(predict_den$x)
      y_obs <- as.numeric(y_obs_future[j,i])
      appximate_density <- approxfun(predict_den)
      if(y_obs < y_min) {y_obs = y_min}
      if(y_obs > y_max) {y_obs = y_max}
      log_pred[j,i] <- log(appximate_density(y_obs))
      emp_CDF[j,i] <- ecdf(x = predictive_samples$y_pred[j,i,])(y_obs)
    }
  }
  # for (i in c(1: (K-1))){
  #   for (k in c((i+1):K)){
  #     for (j in c(1:t_pred)){
  #       # predict_den <- pmpp::kde2D(data = cbind(predictive_samples$y_pred[j,i,],
  #       #                                             predictive_samples$y_pred[j,k,])
  #       #                                   )
  #       # KernSmooth::bkde2D(x = cbind(predictive_samples$y_pred[j,i,],
  #       #                              predictive_samples$y_pred[j,k,]), bandwidth = ?)
  #       predict_den <- MASS::kde2d(x = predictive_samples$y_pred[j,i,],
  #                                 y = predictive_samples$y_pred[j,k,], n = 100)
  #
  #       log_biv_pred[(i-1)*K+k-i,j] <- log(fields::interp.surface(predict_den, cbind(y_obs_future[j,i], y_obs_future[j,k])))
  #     }
  #
  #   }
  # }

  return(list(log_pred = log_pred,
              MSFE = (apply(predictive_samples$y_pred, MARGIN = c(1,2), FUN = mean) - y_obs_future)^2,
              MAFE = abs(apply(predictive_samples$y_pred, MARGIN = c(1,2), FUN = mean) - y_obs_future),
              # log_biv_pred = log_biv_pred,
              emp_CDF = emp_CDF))
}

#' @export
recursive_forecast <- function(Chain, y0 = NULL, y_future, t_pred = 12, reestimated = F){
  K <- Chain$K
  p <- Chain$p

  if (! ncol(y_future) == K) { stop("ncol(y_future) != K") }

  if (is.null(y0)){
    y0 <- Chain$y
  }

  t_length <- nrow(y_future)
  if (t_length < t_pred)  { stop("ncol(y_future) < t_pred)") }
  max_rolling <- t_length - t_pred

  y_combine <- rbind( tail(y0, p) , y_future)
  sum_log_pred <- matrix(0, nrow = t_pred, ncol = K)
  sum_MSFE <- matrix(0, nrow = t_pred, ncol = K)
  sum_MAFE <- matrix(0, nrow = t_pred, ncol = K)
  PIT <- array(0, dim = c(t_pred, K, max_rolling))
  for (time_id in c(1:max_rolling)){
    if (time_id %% 100 == 0 ) cat("At rolling ", time_id)
    y_current <- y_combine[c(time_id:(time_id+p-1)), ]
    y_obs_future <- y_combine[c((time_id+p):(time_id+p+t_pred-1)), ]
    forecast_err <- forecast_density(Chain, y_current, y_obs_future)
    sum_log_pred <- sum_log_pred + forecast_err$log_pred
    sum_MSFE <- sum_MSFE + forecast_err$MSFE
    sum_MAFE <- sum_MAFE + forecast_err$MAFE
    PIT[,,time_id] <- forecast_err$emp_CDF
  }

  return( list( LP = sum_log_pred / max_rolling,
                MSFE = sum_MSFE / max_rolling,
                MAFE = sum_MAFE / max_rolling,
                PIT = PIT))
}
