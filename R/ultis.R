##########################################################################
# Adatptive metropolis hasting #
##########################################################################
TARGACCEPT = 0.3
batchlength = 10

#' @export
adaptamount <- function(iteration){
  return( min( 0.01, 1.0 / sqrt(iteration) ) );
}

# This function is borrowed from https://github.com/FK83/bvarsv
#' @export
getmix <- function(){
  q <- c(0.00730, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.25750)      # probabilities
  m <- c(-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -1.08819) # means
  u2 <- c(5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261)    #variances
  return(list(q=q,m=m,u2=u2))
}
##########################################################################

##########################################################################
# Utility functions  #
##########################################################################
#' @export
makeRegressor <- function(y, y0, t_max, K, p){
  xt <- NULL
  yt = t(y)

  if (is.null(y0)){
    y0 <- matrix(0, ncol = K, nrow = p)

  }
  y0 <- rbind(y0, y)
  for (i in c(1:p)){
    xt <- cbind(xt, rbind(1, vec( t(y0)[,(p+i-1):i])))
  }
  for (i in c( p:(t_max-1))){
    xt <- cbind(xt, rbind(1, vec(yt[,i:(i-p+1)])))
  }
  return(xt)
}

#' @export
Vec_to_Mat <- function(b, K, p){
  return(matrix(b, nrow = K))
}

#' @export
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

#' @export
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

#' @export
sample_A_ele <- function(ysub, xsub, a_sub, V_a_sub){
  t_max = length(ysub)
  n = nrow(xsub)

  a_post = rep(0, n)
  V_a_sub_inv = solve(V_a_sub)

  V_a_post_inv <- V_a_sub_inv + xsub %*% t(xsub)
  V_a_post <- solve(V_a_post_inv)
  a_post <- V_a_post %*% ( V_a_sub_inv %*% a_sub + xsub %*% ysub)
  return(a_post  + t(chol(V_a_post)) %*% rnorm(n))
}

#' @export
sample_h_ele <- function(ytilde, sigma_h = 0.0001*diag(K),
                         h = matrix(0, nrow = t_max, ncol = K), K, t_max){
  tmp <- getmix()
  q <- tmp$q
  m_mean <- tmp$m
  u2 <- tmp$u2
  Zs <- matrix(1,t_max,1) %x% diag(K)
  # # log Gaussian prior
  # sigma_prmean <- prior$sigma_OLS in log scale
  # sigma_prvar <- k_sig*diag(M)

  sigma_prmean <- rep(0,K) # log(sigma_OLS) assumed to be 0
  sigma_prvar <- diag(K)

  k_W = 0.03 # mean of sigma_h
  W_prvar <- (K) # Wprior
  W_prmean <- ((k_W)^2)*(K)*diag(K)


  aux <- sigmahelper4(t(ytilde^2), q, m_mean, u2, h, Zs, sigma_h, sigma_prmean, sigma_prvar)
  h <- aux$Sigtdraw
  sqrtvol <- aux$sigt

  sse_2 <- h[,2:t_max] - h[,1:(t_max-1)]
  sse_2 <- sse_2 %*% t(sse_2)
  sigma_h <- solve(rWishart( 1, t_max + W_prvar, solve(sse_2 + W_prmean))[,,1])
  #sigma_h <- LaplacesDemon::rinvwishart(nu = t_max + W_prvar, S = sse_2 + W_prmean)

  aux$sigma_h <- sigma_h
  return(aux)
}
