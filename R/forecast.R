
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
get_forecast <- function(Chain, t_pred = 12, Nfsample = NULL){
  K <- Chain$K
  p <- Chain$p

  y <- Chain$y
  y0 <- tail(y, p) # Get last p obs

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
