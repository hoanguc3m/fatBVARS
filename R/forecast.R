
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
get_forecast <- function(Chain, y0 = NULL, t_pred = 12, t_current = NULL, Nfsample = NULL, thin = 1){
  K <- Chain$K
  p <- Chain$p

  if( is.null(y0)) {
    # y <- Chain$y
    y0 <- tail(Chain$y, p) # Get last p obs
  }
  if( is.null(t_current)) {
    t_current <- which(Chain$y == tail(y0,p)[p,1])
  }

  dist <- Chain$dist
  SV <- Chain$prior$SV

  ndraws <- nrow(Chain$mcmc)
  mcmc <- Chain$mcmc[seq(from = 1, to = ndraws, by = thin),]
  if (is.null(Nfsample)){
    Nfsample <- nrow(mcmc)
    frange <- c(1:Nfsample)
  } else {
    frange <- sample(nrow(mcmc), Nfsample, replace = T)
  }

  y_pred <- array(NA, dim = c(t_pred, K, Nfsample))
  ymean_pred <- array(NA, dim = c(t_pred, K, Nfsample)) # Predictive mean
  yvar_pred <- array(NA, dim = c(t_pred, 0.5*K*(K+1), Nfsample)) # Predictive Chol Sigma
  vol_pred <- array(NA, dim = c(t_pred, K, Nfsample)) # Predictive diag(Sigma)

  B_mat <- get_post(mcmc, element = "B")
  A_mat <- get_post(mcmc, element = "a")
  Sigma_mat <- get_post(mcmc, element = "sigma") # No SV
  Sigma_H <- sqrt(get_post(mcmc, element = "sigma_h")) # SV
  H_mat <- get_post(mcmc, element = "h")
  if (ncol(H_mat)>0) H_mat <- H_mat[,( (t_current - 1) * K +1):(t_current * K)]
  Nu_mat <- get_post(mcmc, element = "nu")
  Gamma_mat <- get_post(mcmc, element = "gamma")
  count_id <- 0
  for (i in frange){
    count_id <- count_id+1
    b0 <- B_mat[i,]
    a0 <- A_mat[i,]

    if (SV) {
        sigma_h <- Sigma_H[i,]
        if (K == 1) {
          h <- H_mat[i]
        } else {
          h <- H_mat[i,]
        }

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
        if (dist == "multiOrthStudent") {
          nu <- Nu_mat[i,]
          pred_tmp <- sim.VAR.multiOrthStudent.SV(K = K, p = p, t_max = t_pred,
                                              b0 = b0, a0 = a0, h = h, sigma_h = sigma_h,
                                              y0 = y0, nu = nu,
                                              seednum = sample(1:100000,1), burn_in = 0)
        }
        if (dist == "Hyper.multiOrthStudent") {
          nu <- Nu_mat[i,]
          gamma <- Gamma_mat[i,]
          pred_tmp <- sim.VAR.Hyper.multiOrthStudent.SV(K = K, p = p, t_max = t_pred,
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
        if (dist == "multiOrthStudent") {
          nu <- Nu_mat[i,]
          pred_tmp <- sim.VAR.multiOrthStudent.novol(K = K, p = p, t_max = t_pred,
                                                 b0 = b0, a0 = a0, h = h,
                                                 y0 = y0, nu = nu,
                                                 seednum = sample(1:100000,1), burn_in = 0)
        }
        if (dist == "Hyper.multiOrthStudent") {
          nu <- Nu_mat[i,]
          gamma <- Gamma_mat[i,]
          pred_tmp <- sim.VAR.Hyper.multiOrthStudent.novol(K = K, p = p, t_max = t_pred,
                                                       b0 = b0, a0 = a0, h = h,
                                                       y0 = y0, nu = nu, gamma = gamma,
                                                       seednum = sample(1:100000,1), burn_in = 0)
        }
      }

    y_pred[,,count_id] <- pred_tmp$y
    ymean_pred[,,count_id] <- pred_tmp$y_mean
    yvar_pred[,,count_id] <- pred_tmp$y_var
    vol_pred[,,count_id] <- pred_tmp$volatility

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
  if (is.null(t_current)) t_current = which(Chain$y == tail(y_current,p)[p,1] )

  t_pred = nrow(y_obs_future)
  Nfsample = 10000
  predictive_samples <- get_forecast(Chain = Chain, y0 = tail(y_current,p), t_pred = t_pred, t_current = t_current, thin = 10, Nfsample = Nfsample) # Nfsample

  #############################################################
  # Non cummulative
  #############################################################
  log_pred<- matrix(NA, nrow = t_pred, ncol = K)
  bilog_pred<- matrix(NA, nrow = t_pred, ncol = K*(K-1)*0.5)
  emp_CDF<- matrix(NA, nrow = t_pred, ncol = K)

  for (i in c(1:K)){
    for (j in c(1:t_pred)){
      y_obs <- as.numeric(y_obs_future[j,i])
      log_pred[j,i] <- log(mean(dnorm(y_obs, mean = predictive_samples$ymean_pred[j,i,], sd = sqrt(predictive_samples$vol_pred[j,i,]))))
      emp_CDF[j,i] <- ecdf(x = predictive_samples$y_pred[j,i,])(y_obs)
    }
  }
  if (K > 1){
    Sigma2 <- array(0, dim = c(K, K, Nfsample))
    Sigma <- matrix(0, nrow = K, ncol = K)
    for (j in c(1:t_pred)){
      for (sample_id in c(1:Nfsample)){
        Sigma[lower.tri(Sigma, diag = T)] <- predictive_samples$yvar_pred[j,,sample_id]
        Sigma2[,,sample_id] <- Sigma %*% t(Sigma)
      }

      for (i in c(1: (K-1))){
        for (k in c((i+1):K)){
          pred_bidens <- rep(NA, Nfsample)
          for (sample_id in c(1:Nfsample)){
            pred_bidens[sample_id] <- dmvn(X = c(y_obs_future[j,i], y_obs_future[j,k]),
                                           mu = predictive_samples$ymean_pred[j, c(i,k),sample_id],
                                           sigma = Sigma2[c(i,k),c(i,k), sample_id], log = FALSE, isChol = FALSE)
          }
          bilog_pred[j,(i-1)*K + k-i*(i+1)*0.5] <- log(mean(pred_bidens))
        }
      }
    }
  }

  #############################################################
  # Cummulative
  #############################################################
  clog_pred<- matrix(NA, nrow = t_pred, ncol = K)
  cbilog_pred<- matrix(NA, nrow = t_pred, ncol = K*(K-1)*0.5)
  cemp_CDF<- matrix(NA, nrow = t_pred, ncol = K)
  cy_pred <- apply(predictive_samples$y_pred, MARGIN = c(2,3), FUN = cumsum)
  cy_obs_future <- apply(y_obs_future, MARGIN = 2, FUN = cumsum)
  # for (i in c(1:K)){
  #   for (j in c(1:t_pred)){
  #     y_obs <- as.numeric(cy_obs_future[j,i])
  #     predict_den <- stats::density(cy_pred[j,i,], adjust = 3)
  #     id <- which(predict_den$y == 0);
  #     if (length(id)>0) {
  #       for (fix in id){
  #         predict_den$y[fix] <- mean(predict_den$y[(fix-1):(fix+1)])
  #       }
  #     }
  #     appximate_density <- smooth.spline(x = predict_den$x, y = log(predict_den$y))
  #     pred_dens <- predict(appximate_density, y_obs)$y
  #     clog_pred[j,i] <- pred_dens
  #     cemp_CDF[j,i] <- ecdf(x = cy_pred[j,i,])(y_obs)
  #   }
  # }

  CRPS <- matrix(NA, nrow = t_pred, ncol = K)
  SCRPS <- matrix(NA, nrow = t_pred, ncol = K)

  for (i in c(1:K)){
    for (j in c(1:t_pred)){
      y_obs <- as.numeric(y_obs_future[j,i])
      samples1 <- predictive_samples$y_pred[j,i,]
      samples2 <- sample(samples1, size = length(samples1), replace = T)
      CRPS[j,i] <- - mean(abs(samples1 - y_obs)) + 0.5*mean(abs(samples1 - samples2))
      SCRPS[j,i] <- - mean(abs(samples1 - y_obs)) / mean(abs(samples1 - samples2)) - 0.5*log(mean(abs(samples1 - samples2)))
    }
  }

  return(list(log_pred = log_pred,
              bilog_pred = bilog_pred,
              emp_CDF = emp_CDF,
              MSFE = (apply(predictive_samples$y_pred, MARGIN = c(1,2), FUN = mean) - y_obs_future)^2,
              MAFE = abs(apply(predictive_samples$y_pred, MARGIN = c(1,2), FUN = mean) - y_obs_future),
              # predictive_samples = predictive_samples,
              t_current = t_current,
              clog_pred = clog_pred,
              cbilog_pred = cbilog_pred,
              cemp_CDF = cemp_CDF,
              cMSFE = (apply(cy_pred, MARGIN = c(1,2), FUN = mean) - cy_obs_future)^2,
              cMAFE = abs(apply(cy_pred, MARGIN = c(1,2), FUN = mean) - cy_obs_future),
              CRPS = CRPS,
              SCRPS = SCRPS
  ))
}

#' Recursive forecast of BVAR model
#'
#' This function returns a recursive forecast of future observations of BVAR-SV-fatTail model.
#' @param y the whole data matrix T x k.
#' @param t_start the start point for forecast.
#' @param t_pred The time prediction horizon.
#' @param K The number of variables in y.
#' @param p The number of lags in BVAR model.
#' @return The recursive forecast density of future observations of BVAR model
#' @export
#' @examples
#' \dontrun{
#' recuresive_model1 <- recursive_forecast(y, t_start, t_pred, K, p, dist = "Hyper.multiStudent", SV = T)
#' }
#'
#' @export
recursive_forecast <- function(y, t_start = 100, t_pred = 12, K, p, dist = "Hyper.multiStudent", SV = T, numCores = NULL){
    t_max = nrow(y)
    max_rolling <- t_max - (t_start + t_pred) + 1

    log_pred <- array(NA, dim = c(t_pred, K, max_rolling))
    bilog_pred <- array(NA, dim = c(t_pred, K*(K-1)*0.5, max_rolling))
    MSFE <- array(NA, dim = c(t_pred, K, max_rolling))
    MAFE <- array(NA, dim = c(t_pred, K, max_rolling))
    PIT <- array(NA, dim = c(t_pred, K, max_rolling))

    if (is.null(numCores)) numCores <- parallel::detectCores()*0.5
    out_recursive <- parallel::mclapply(1:max_rolling,
                       FUN = function(j) {
                         time_current <- j + t_start - 1
                         y_current <- y[c(1:time_current), ]
                         prior <- get_prior(tail(y_current, time_current - p), p = p, dist=dist, SV = SV)
                         inits <- get_init(prior, samples = 15000, burnin = 5000, thin = 1)
                         if (SV) {
                           Chain <- BVAR.SV(y = tail(y_current, time_current - p), K = K, p = p, dist = dist,
                                            y0 = head(y_current, p), prior = prior, inits = inits)
                         } else {
                           Chain <- BVAR.novol(y = tail(y_current, time_current - p), K = K, p = p, dist = dist,
                                               y0 = head(y_current, p), prior = prior, inits = inits)

                         }
                         y_obs_future <- y[c((time_current+1):(time_current+t_pred)), ]
                         forecast_err <- forecast_density(Chain = Chain, y_obs_future = y_obs_future)
                         list(time_id = j, forecast_err = forecast_err)
                       },
                       mc.cores = numCores)

    for (forecast_element in out_recursive){
      time_id = forecast_element$time_id
      forecast_err = forecast_element$forecast_err
      log_pred[,,time_id] <- forecast_err$log_pred
      bilog_pred[,,time_id] <- forecast_err$bilog_pred
      MSFE[,,time_id] <- forecast_err$MSFE
      MAFE[,,time_id] <- forecast_err$MAFE
      PIT[,,time_id] <- forecast_err$emp_CDF
    }

    mlog_pred <- apply(log_pred, MARGIN = c(1,2), FUN = mean)
    mbLP <- apply(bilog_pred, MARGIN = c(1,2), FUN = mean)
    mMSFE <- apply(MSFE, MARGIN = c(1,2), FUN = mean)
    mMAFE <- apply(MAFE, MARGIN = c(1,2), FUN = mean)

  return( list( mLP = mlog_pred,
                mbLP = mbLP,
                mMSFE = mMSFE,
                mMAFE = mMAFE,
                log_pred = log_pred,
                bilog_pred = bilog_pred,
                MSFE = MSFE,
                MAFE = MAFE,
                PIT = PIT))
}

#' @export
recursive_seperate <- function(y, t_start = 100, t_pred = 12, K, p, dist = "Hyper.multiStudent", SV = T, outname = NULL){
  t_max = nrow(y)
  if (is.null(outname)) outname = paste("Recursive_", dist, "_", SV, "_", t_start+1, ".RData", sep = "")
  time_current <- t_start
  y_current <- matrix(y[c(1:time_current), ], ncol = K)
  prior <- get_prior(tail(y_current, time_current - p), p = p, dist=dist, SV = SV)
  inits <- get_init(prior, samples = 15000, burnin = 5000, thin = 1)
  if (SV) {
    Chain <- BVAR.SV(y = tail(y_current, time_current - p), K = K, p = p, dist = dist,
                     y0 = head(y_current, p), prior = prior, inits = inits)
  } else {
    Chain <- BVAR.novol(y = tail(y_current, time_current - p), K = K, p = p, dist = dist,
                        y0 = head(y_current, p), prior = prior, inits = inits)
  }
  y_obs_future <- matrix(y[c((time_current+1):(time_current+t_pred)), ],ncol = K)
  forecast_err <- forecast_density(Chain = Chain, y_obs_future = y_obs_future)
  out_recursive <- list(time_id = time_current+1,
                        forecast_err = forecast_err,
                        dist = dist,
                        SV = SV,
                        y_current = y_current,
                        y_obs_future = y_obs_future)
  save(out_recursive, file = outname)
  return( out_recursive)

}
# recursive_forecast <- function(y, t_start = 100, t_pred = 12, K, p, dist = "Hyper.multiStudent", SV = T, core = NULL){
#   t_max = nrow(y)
#   max_rolling <- t_max - (t_start + t_pred) + 1
#
#   log_pred <- array(NA, dim = c(t_pred, K, max_rolling))
#   bilog_pred <- array(NA, dim = c(t_pred, K*(K-1)*0.5, max_rolling))
#   MSFE <- array(NA, dim = c(t_pred, K, max_rolling))
#   MAFE <- array(NA, dim = c(t_pred, K, max_rolling))
#   PIT <- array(NA, dim = c(t_pred, K, max_rolling))
#
#   for (time_id in c(1:max_rolling)){
#     if (time_id %% 10 == 0 ) cat("At rolling ", time_id, " \t")
#     time_current <- time_id + t_start - 1
#     y_current <- y[c(1:time_current), ]
#     prior <- get_prior(tail(y_current, time_current - p), p = p, dist=dist, SV = SV)
#     inits <- get_init(prior, samples = 15000, burnin = 5000, thin = 1)
#     if (SV) {
#       Chain <- BVAR.SV(y = tail(y_current, time_current - p), K = K, p = p, dist = dist,
#                        y0 = head(y_current, p), prior = prior, inits = inits)
#     } else {
#       Chain <- BVAR.novol(y = tail(y_current, time_current - p), K = K, p = p, dist = dist,
#                           y0 = head(y_current, p), prior = prior, inits = inits)
#
#     }
#     y_obs_future <- y[c((time_current+1):(time_current+t_pred)), ]
#     forecast_err <- forecast_density(Chain = Chain, y_obs_future = y_obs_future)
#     log_pred[,,time_id] <- forecast_err$log_pred
#     bilog_pred[,,time_id] <- forecast_err$bilog_pred
#     MSFE[,,time_id] <- forecast_err$MSFE
#     MAFE[,,time_id] <- forecast_err$MAFE
#     PIT[,,time_id] <- forecast_err$emp_CDF
#   }
#   mlog_pred <- apply(log_pred, MARGIN = c(1,2), FUN = mean)
#   mbLP <- apply(bilog_pred, MARGIN = c(1,2), FUN = mean)
#   mMSFE <- apply(MSFE, MARGIN = c(1,2), FUN = mean)
#   mMAFE <- apply(MAFE, MARGIN = c(1,2), FUN = mean)
#
#   return( list( mLP = mlog_pred,
#                 mbLP = mbLP,
#                 mMSFE = mMSFE,
#                 mMAFE = mMAFE,
#                 log_pred = log_pred,
#                 bilog_pred = bilog_pred,
#                 MSFE = MSFE,
#                 MAFE = MAFE,
#                 PIT = PIT))
# }
