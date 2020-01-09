
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
get_forecast <- function(Chain, y0 = NULL, t_pred = 12, t_current = NULL, Nfsample = NULL){
  K <- Chain$K
  p <- Chain$p

  if( is.null(y0)) {
    # y <- Chain$y
    y0 <- tail(Chain$y, p) # Get last p obs
  }
  if( is.null(t_current)) {
    t_current <- which(Chain$y == y0[p,1])
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

  B_mat <- get_post(mcmc, element = "B")
  A_mat <- get_post(mcmc, element = "a")
  Sigma_mat <- get_post(mcmc, element = "sigma") # No SV
  Sigma_H <- sqrt(get_post(mcmc, element = "sigma_h")) # SV
  H_mat <- get_post(mcmc, element = "h")
  if (ncol(H_mat)>0) H_mat <- H_mat[,( (t_current - 1) * K +1):(t_current * K)]
  Nu_mat <- get_post(mcmc, element = "nu")
  Gamma_mat <- get_post(mcmc, element = "gamma")
  for (i in frange){
    b0 <- B_mat[i,]
    a0 <- A_mat[i,]

    if (SV) {
        sigma_h <- Sigma_H[i,]
        h <- H_mat[i,]
        if (dist == "Gaussian") {
          pred_tmp <- sim.VAR.Gaussian.SV(K = K, p = p, t_max = t_pred,
                                             b0 = b0, a0 = a0, h = h, sigma_h = sigma_h,
                                             y0 = y0,
                                             seednum = sample(1:100000,1), burn_in = 0)
        }

        if (dist == "Student") {
          nu <- Nu_mat[i]
          pred_tmp <- sim.VAR.Student.SV(K = K, p = p, t_max = t_pred,
                                            b0 = b0, a0 = a0, h = h, sigma_h = sigma_h,
                                            y0 = y0, nu = nu,
                                            seednum = sample(1:100000,1), burn_in = 0)
        }
        if (dist == "Hyper.Student") {
          nu <- Nu_mat[i]
          gamma <- Gamma_mat[i,]
          pred_tmp <- sim.VAR.Hyper.Student.SV(K = K, p = p, t_max = t_pred,
                                                  b0 = b0, a0 = a0, h = h, sigma_h = sigma_h,
                                                  y0 = y0, nu = nu, gamma = gamma,
                                                  seednum = sample(1:100000,1), burn_in = 0)
        }
        if (dist == "multiStudent") {
          nu <- Nu_mat[i,]
          pred_tmp <- sim.VAR.multiStudent.SV(K = K, p = p, t_max = t_pred,
                                                 b0 = b0, a0 = a0, h = h, sigma_h = sigma_h,
                                                 y0 = y0, nu = nu,
                                                 seednum = sample(1:100000,1), burn_in = 0)
        }
        if (dist == "Hyper.multiStudent") {
          nu <- Nu_mat[i,]
          gamma <- Gamma_mat[i,]
          pred_tmp <- sim.VAR.Hyper.multiStudent.SV(K = K, p = p, t_max = t_pred,
                                                       b0 = b0, a0 = a0, h = h, sigma_h = sigma_h,
                                                       y0 = y0, nu = nu, gamma = gamma,
                                                       seednum = sample(1:100000,1), burn_in = 0)
        }
      } else {
        sigma <- Sigma_mat[i,]
        h <- log(sigma)
        if (dist == "Gaussian") {
            pred_tmp <- sim.VAR.Gaussian.novol(K = K, p = p, t_max = t_pred,
                                               b0 = b0, a0 = a0, h = h,
                                               y0 = y0,
                                               seednum = sample(1:100000,1), burn_in = 0)
          }
        if (dist == "Student") {
          nu <- Nu_mat[i]
          pred_tmp <- sim.VAR.Student.novol(K = K, p = p, t_max = t_pred,
                                             b0 = b0, a0 = a0, h = h,
                                             y0 = y0, nu = nu,
                                             seednum = sample(1:100000,1), burn_in = 0)
        }
        if (dist == "Hyper.Student") {
          nu <- Nu_mat[i]
          gamma <- Gamma_mat[i,]
          pred_tmp <- sim.VAR.Hyper.Student.novol(K = K, p = p, t_max = t_pred,
                                            b0 = b0, a0 = a0, h = h,
                                            y0 = y0, nu = nu, gamma = gamma,
                                            seednum = sample(1:100000,1), burn_in = 0)
        }
        if (dist == "multiStudent") {
          nu <- Nu_mat[i,]
          pred_tmp <- sim.VAR.multiStudent.novol(K = K, p = p, t_max = t_pred,
                                            b0 = b0, a0 = a0, h = h,
                                            y0 = y0, nu = nu,
                                            seednum = sample(1:100000,1), burn_in = 0)
        }
        if (dist == "Hyper.multiStudent") {
          nu <- Nu_mat[i,]
          gamma <- Gamma_mat[i,]
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

#' Forecast density of future observations of BVAR model
#'
#' This function returns a forecast density of future observations of BVAR-SV-fatTail model.
#' @param Chain The fatBVARSV object from command BVAR.
#' @param y_current The current values of y.
#' @param y_obs_future The future observable values of y.
#' @param t_pred The time prediction horizon.
#' @return The forecast density of future observations of BVAR model
#' @export
#' @examples
#' \dontrun{
#' forecast_dens1 <- forecast_density(Chain1, y_obs_future)
#' }
#'
#' @export
forecast_density <- function(Chain, y_current = NULL, y_obs_future, t_current = NULL){
  K <- Chain$K
  p <- Chain$p

  if (! ncol(y_obs_future) == K) { stop("ncol(y_obs_future) != K") }
  if (is.null(y_current)) y_current <- Chain$y
  if (is.null(t_current)) t_current = which(Chain$y == y_current[p,1])

  t_pred = nrow(y_obs_future)
  predictive_samples <- get_forecast(Chain = Chain, y0 = y_current, t_pred = t_pred, t_current = t_current) # Nfsample

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

      # appximate_density <- approxfun(predict_den)
      # log_pred[j,i] <- log(appximate_density(y_obs))

      appximate_density <- smooth.spline(x = predict_den$x, y = predict_den$y)
      pred_dens <- predict(appximate_density, y_obs)$y
      if ( pred_dens <= 0) {
        if ( abs(y_obs - y_min) < abs(y_obs - y_max) ) {
          pred_dens <- predict(appximate_density, y_min)$y
        } else {
          pred_dens <- predict(appximate_density, y_max)$y
        }
      }
      # if(y_obs < y_min) {y_obs = y_min}
      # if(y_obs > y_max) {y_obs = y_max}
      log_pred[j,i] <- log(pred_dens)
      # if(is.na(log_pred[j,i])) {
      #   cat("t_current : ", t_current, " ", i, " ", j, " ", y_obs, " ", y_min, " ", y_max)
      #   stop(" NA log ")
      # }
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
              predictive_samples = predictive_samples,
              t_current = t_current,
              # log_biv_pred = log_biv_pred,
              emp_CDF = emp_CDF))
}

#' Forecast density of future observations of BVAR model
#'
#' This function returns a forecast density of future observations of BVAR-SV-fatTail model.
#' @param Chain The fatBVARSV object from command BVAR.
#' @param y_current The current values of y.
#' @param y_obs_future The future observable values of y.
#' @param t_pred The time prediction horizon.
#' @return The forecast density of future observations of BVAR model
#' @export
#' @examples
#' \dontrun{
#' forecast_dens1 <- forecast_density(Chain1, y_obs_future)
#' }
#'
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
  max_rolling <- t_length - t_pred + 1

  y_combine <- rbind( tail(y0, p) , y_future)
  sum_log_pred <- matrix(0, nrow = t_pred, ncol = K)
  sum_MSFE <- matrix(0, nrow = t_pred, ncol = K)
  sum_MAFE <- matrix(0, nrow = t_pred, ncol = K)
  PIT <- array(0, dim = c(t_pred, K, max_rolling))
  for (time_id in c(1:max_rolling)){
    if (time_id %% 10 == 0 ) cat("At rolling ", time_id, " \t")
    y_current <- y_combine[c(time_id:(time_id+p-1)), ]
    y_obs_future <- y_combine[c((time_id+p):(time_id+p+t_pred-1)), ]
    forecast_err <- forecast_density(Chain = Chain, y_current = y_current,
                                     y_obs_future = y_obs_future,
                                     t_current = which(Chain$y == y_current[p,1]) )
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
