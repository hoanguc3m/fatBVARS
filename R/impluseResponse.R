######################################################################################

#' Impulse response function of BVAR model
#'
#' This function returns a impulse response function of BVAR-SV-fatTail model.
#' @param Chain The fatBVARSV object from command BVAR.
#' @param impulse.variable The position of impulse variable
#' @param response.variable The position of reponse variable
#' @param time The time at which impulse response function is calculated
#' @param n.ahead Maximal time between impulse and response (defaults to 20).
#' @param scenario If 1, there is no orthogonalizaton, and the shock size corresponds to one unit of the impulse variable.
#' If scenario is either 2 (the default) or 3, the error term variance-covariance matrix is orthogonalized via Cholesky decomposition.
#' For scenario = 2, the Cholesky decomposition of the error term VCV matrix at time point t is used.
#' scenario = 3 is the variant used in Del Negro and Primiceri (2015).
#' Here, the diagonal elements are set to their averages over time, whereas the off-diagonal elements are specific to time t.
#' See the bvars package.
#' @param draw.plot Show plot
#' @return A list of two elements: Contemporaneous impulse responses and Matrix of simulated impulse responses
#' @export
#' @examples
#' \dontrun{
#' irf_obj <- get_irf(Chain = Chain, impulse.variable = 1, response.variable = 2,
#'                    t = NULL, n.ahead = 20, scenario = 2, draw.plot = TRUE)
#' }
#' @references bvars package https://cran.r-project.org/web/packages/bvarsv
get_irf <- function(Chain, impulse.variable = 1, response.variable = 2,
                    atT = NULL, n.ahead = 20, scenario = 2, draw.plot = TRUE){
  # TODO: Check IRF in case of fat tail.
  K <- Chain$K
  p <- Chain$p

  ndraws <- nrow(Chain$mcmc)
  t_max <- nrow(Chain$y)

  dist <- Chain$dist
  SV <- Chain$prior$SV

  if (is.null(atT)){
    atT = t_max
  }

  if (SV) {
    sig <- exp(get_post(Chain$mcmc, element = "h")*0.5)
    if (scenario == 3){
      sig <- matrix(apply( sig , 2, mean), nrow = K)
      sig <- apply(sig, MARGIN = 1, mean)
      sig <- diag(sig)
    } else {
      sig <- sig[, c( ((atT-1) * K+1) : (atT * K))]
    }
  } else {
    sig <- get_post(Chain$mcmc, element = "sigma")
    if (scenario == 3){
      sig <- apply( get_post(Chain$mcmc, element = "sigma") , 2, mean)
      sig <- diag(sig)
    }
  }
  H.chol <- NULL
  beta <- get_post(Chain$mcmc, element = "B")[,-c(1:K)] # Remove constant
  out <- matrix(0, ndraws, n.ahead + 1)

  # Compute IR for each MC iteration
  for (j in c(1:ndraws)){

    # Cholesky of VCV matrix, depending on specification
    if (scenario > 1){
      A0 <- matrix(0, nrow = K, ncol = K)
      A0[upper.tri(A0)] <- get_post(Chain$mcmc[j,], element = "a")
      A0 <- t(A0)
      diag(A0) <- 1

      if (scenario == 2){ # Standard orthogonalization
        H.chol <- A0 %*% diag(sig[j,])
      } else { # Orthogonalization as in Primiceri
        H.chol <- A0 %*% sig
      }
    }

    # Compute Impulse Responses
    aux <- IRFmats(B = matrix(beta[j,], nrow = K),
                   H.chol = H.chol, n.ahead = n.ahead)
    aux <- aux[response.variable, seq(from = impulse.variable, by = K, length = n.ahead + 1)]

    out[j,] <- aux
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
  return(list(contemporaneous = out[,1], irf = out[,-1])  )
}

#' @export
IRFmats <- function(B, H.chol = NULL, n.ahead = 20){
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
  Phi <- matrix(0, K, K*(n.ahead+1))
  for (s in 0:n.ahead){
    aux <- matrixcalc::matrix.power(Bc, s)[1:K, 1:K]
    if (!is.null(H.chol)) aux <- aux %*% H.chol
    Phi[,(s*K+1):((s+1)*K)] <- aux
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
