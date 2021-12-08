######################################################################################

#' Impulse response function of BVAR model
#'
#' This function returns a impulse response function of BVAR-SV-fatTail model.
#' @param Chain The fatBVARSV object from command BVAR.
#' @param impulse.variable The position of impulse variable
#' @param response.variable The position of reponse variable
#' @param atT The time at which impulse response function is calculated
#' @param n.ahead Maximal time between impulse and response (defaults to 20).
#' @param draw.plot Show plot
#' @return A list of two elements: Contemporaneous impulse responses and Matrix of simulated impulse responses
#' @export
#' @examples
#' \dontrun{
#' irf_obj <- get_irf(Chain = Chain, impulse.variable = 1, response.variable = 2,
#'                    t = NULL, n.ahead = 20, draw.plot = FALSE)
#' }
#' @references bvars package https://cran.r-project.org/web/packages/bvarsv
get_irf <- function(Chain, impulse.variable = 1, response.variable = 2,
                    atT = NULL, n.ahead = 20, draw.plot = FALSE){
  # TODO: Check IRF in case of fat tail.
  K <- Chain$K
  p <- Chain$p


  ndraws <- nrow(Chain$mcmc)
  t_max <- nrow(Chain$y)

  if (is.null(atT)){
    atT = t_max
  }


  dist <- Chain$dist
  SV <- Chain$prior$SV


  if (SV) {
    sig <- exp(get_post(Chain$mcmc, element = "h")*0.5)
    sig <- sig[, c( ((atT-1) * K+1) : ((atT) * K))]
  } else {
    sig <- get_post(Chain$mcmc, element = "sigma")
  }
  if (!(dist == "Gaussian")){
    if (dist == "Student" || dist == "Skew.Student"){
      W <- sqrt(get_post(Chain$mcmc, element = "w")[, atT])
    } else {
      W <- sqrt(get_post(Chain$mcmc, element = "w")[, c( ((atT-1) * K+1) : ((atT) * K))])
    }
  }


  H.chol <- array(diag(K), dim = c(K,K, n.ahead+1))
  beta <- get_post(Chain$mcmc, element = "B")[,-c(1:K)] # Remove constant
  out <- matrix(0, ndraws, n.ahead + 1)
  out_all <- array(0, c(K,K,ndraws, n.ahead + 1))

  # Compute IR for each MC iteration
  for (j in c(1:ndraws)){

    # Cholesky of VCV matrix, depending on specification

    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- get_post(Chain$mcmc[j,], element = "a")
    A0 <- t(A0)
    diag(A0) <- 1
    A_inv <- solve(A0)

    s = n.ahead+1
    sigtemp <- sig[j,c(1:K)]
    if (dist == "Gaussian"){
      H.chol[,,s] <- A_inv %*% diag(sigtemp)
    }
    if (dist == "Student" || dist == "Skew.Student"){
      wtemp <- W[j]
      H.chol[,,s] <- wtemp * A_inv %*% diag(sigtemp)
    }

    if (dist == "MT"){
      wtemp <- W[j,c(1:K)]
      H.chol[,,s] <- diag(wtemp) %*% A_inv %*% diag(sigtemp)
    }
    if (dist == "MST"){
      wtemp <- W[j,c(1:K)]
      H.chol[,,s] <- diag(wtemp) %*% A_inv %*% diag(sigtemp)
    }
    if (dist == "OT"){
      wtemp <- W[j,c(1:K)]
      H.chol[,,s] <- A_inv %*% diag(wtemp) %*% diag(sigtemp)
    }
    if (dist == "OST"){
      wtemp <- W[j,c(1:K)]
      H.chol[,,s] <- A_inv %*% diag(wtemp) %*% diag(sigtemp)
    }
    for (s in c(1:(n.ahead))){
      H.chol[,,s] <- H.chol[,,n.ahead+1]
    }


    # Compute Impulse Responses
    aux <- IRFmats(B = matrix(beta[j,], nrow = K),
                   H.chol = H.chol)

    out[j,] <- aux[response.variable, impulse.variable, ]
    out_all[,,j,] <- aux
  }

  # Make plot
  if (draw.plot){
    pdat <- t(apply(out[,-1], 2, function(z) quantile(z, c(0.05, 0.25, 0.5, 0.75, 0.95))))
    xax <- 1:n.ahead
    matplot(x = xax, y = pdat, type = "n", ylab = "", xlab = "Horizon", bty = "n", xlim = c(1, n.ahead))
    polygon(c(xax, rev(xax)), c(pdat[,5], rev(pdat[,4])), col = "grey60", border = NA)
    polygon(c(xax, rev(xax)), c(pdat[,4], rev(pdat[,3])), col = "grey30", border = NA)
    polygon(c(xax, rev(xax)), c(pdat[,3], rev(pdat[,2])), col = "grey30", border = NA)
    polygon(c(xax, rev(xax)), c(pdat[,2], rev(pdat[,1])), col = "grey60", border = NA)
    lines(x = xax, y = pdat[,3], type = "l", col = 1, lwd = 2.5)
    abline(h = 0, lty = 2)
  }
  return(list(irf_all = out,
              mean_irf = apply(out, MARGIN = 2, FUN = mean),
              mean_all = apply(out_all, MARGIN = c(1,2,4), FUN = mean) )  )
}

#' @export
get_irf_backward <- function(Chain, impulse.variable = 1, response.variable = 2,
                    atT = NULL, n.ahead = 20, draw.plot = FALSE){
  # TODO: Check IRF in case of fat tail.
  K <- Chain$K
  p <- Chain$p

  if (is.null(atT)){
    atT = t_max
  }

  if (atT < n.ahead){
    stop("atT < n.ahead.")
  }

  ndraws <- nrow(Chain$mcmc)
  t_max <- nrow(Chain$y)

  dist <- Chain$dist
  SV <- Chain$prior$SV


  if (SV) {
    sig <- exp(get_post(Chain$mcmc, element = "h")*0.5)
    sig <- sig[, c( ((atT-1-n.ahead) * K+1) : (atT * K))]
  } else {
    sig <- get_post(Chain$mcmc, element = "sigma")
  }
  if (!(dist == "Gaussian")){
    W <- sqrt(get_post(Chain$mcmc, element = "w")[, c( ((atT-1-n.ahead) * K+1) : (atT * K))])
  }


  H.chol <- array(diag(K), dim = c(K,K, n.ahead+1))
  beta <- get_post(Chain$mcmc, element = "B")[,-c(1:K)] # Remove constant
  out <- matrix(0, ndraws, n.ahead + 1)
  out_all <- array(0, c(K,K,ndraws, n.ahead + 1))

  # Compute IR for each MC iteration
  for (j in c(1:ndraws)){

    # Cholesky of VCV matrix, depending on specification

    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- get_post(Chain$mcmc[j,], element = "a")
    A0 <- t(A0)
    diag(A0) <- 1
    A_inv <- solve(A0)

    for (s in c(1:(n.ahead+1))){
      sigtemp <- sig[j,c(((n.ahead - s+1) * K+1) : ((n.ahead - s + 2) * K)) ]
      if (dist == "Gaussian"){
        H.chol[,,s] <- A_inv %*% diag(sigtemp)
      }
      if (dist == "MT"){
        wtemp <- W[j,c(((n.ahead - s+1) * K+1) : ((n.ahead - s + 2) * K)) ]
        H.chol[,,s] <- diag(wtemp) %*% A_inv %*% diag(sigtemp)
      }
      if (dist == "MST"){
        wtemp <- W[j,c(((n.ahead - s+1) * K+1) : ((n.ahead - s + 2) * K)) ]
        H.chol[,,s] <- diag(wtemp) %*% A_inv %*% diag(sigtemp)
      }
      if (dist == "OT"){
        wtemp <- W[j,c(((n.ahead - s+1) * K+1) : ((n.ahead - s + 2) * K)) ]
        H.chol[,,s] <- A_inv %*% diag(wtemp) %*% diag(sigtemp)
      }
      if (dist == "OST"){
        wtemp <- W[j,c(((n.ahead - s+1) * K+1) : ((n.ahead - s + 2) * K)) ]
        H.chol[,,s] <- A_inv %*% diag(wtemp) %*% diag(sigtemp)
      }
    }


    # Compute Impulse Responses
    aux <- IRFmats(B = matrix(beta[j,], nrow = K),
                   H.chol = H.chol)

    out[j,] <- aux[response.variable, impulse.variable, ]
    out_all[,,j,] <- aux
  }

  # Make plot
  if (draw.plot){
    pdat <- t(apply(out[,-1], 2, function(z) quantile(z, c(0.05, 0.25, 0.5, 0.75, 0.95))))
    xax <- 1:n.ahead
    matplot(x = xax, y = pdat, type = "n", ylab = "", xlab = "Horizon", bty = "n", xlim = c(1, n.ahead))
    polygon(c(xax, rev(xax)), c(pdat[,5], rev(pdat[,4])), col = "grey60", border = NA)
    polygon(c(xax, rev(xax)), c(pdat[,4], rev(pdat[,3])), col = "grey30", border = NA)
    polygon(c(xax, rev(xax)), c(pdat[,3], rev(pdat[,2])), col = "grey30", border = NA)
    polygon(c(xax, rev(xax)), c(pdat[,2], rev(pdat[,1])), col = "grey60", border = NA)
    lines(x = xax, y = pdat[,3], type = "l", col = 1, lwd = 2.5)
    abline(h = 0, lty = 2)
  }
  return(list(irf_all = out,
              mean_irf = apply(out, MARGIN = 2, FUN = mean),
              mean_all = apply(out_all, MARGIN = c(1,2,4), FUN = mean) )  )
}
#' @export
IRFmats <- function(B, H.chol){
  n.ahead = dim(H.chol)[3] - 1
  # Dimensions
  K <- nrow(B) # nr of variables

  # AR matrices
  p <- ncol(B)/K

  # Construct companion form matrices

  if (p > 1){
    # VAR matrices
    Bc <- matrix(0, K*p, K*p)
    Bc[1:K, ] <- B
    Bc[-(1:K), 1:(K*(p-1))] <- diag(K*(p-1))

  } else {
    Bc <- B
  }

  # Matrix with impulse responses
  Phi <- array(NA, dim = c(K, K, (n.ahead+1)) )
  for (s in 1:(n.ahead+1)){
    aux <- matrixcalc::matrix.power(Bc, s-1)[1:K, 1:K]
    aux <- aux %*% H.chol[,,s]
    Phi[,,s] <- aux
  }
  return(Phi)
}

#' @export
stability <- function(b, K, p){
  # Dimensions
  B <- matrix(b, nrow = K)
  B <- B[,c(2:(K*p+1))]

  if (p > 1){
    # VAR matrices
    Bc <- matrix(0, K*p, K*p)
    Bc[1:K, ] <- B
    Bc[-(1:K), 1:(K*(p-1))] <- diag(K*(p-1))

  } else {
    Bc <- B
  }
  ee = max(abs(eigen(Bc)$values))
  return(ee<1)
}

#' @export
check_post <- function(Chain){
  b_mat <- get_post(Chain, element = "B")
  stabi <- rep(F, nrow(b_mat))
  for (i in c(1:nrow(b_mat))){
    b <- b_mat[i,]
    stabi[i] <- stability(b = b, K = K, p = p)
  }
  return(sum(stabi))
}
