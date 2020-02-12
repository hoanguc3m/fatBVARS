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
reprow<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

#' @export
repcol<-function(x,n){
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


  sigma_prmean <- rep(0,K) # mean h_0
  sigma_prvar <- 10*diag(K)   # variance h_0

  aux <- sigmahelper4(t(ytilde^2), q, m_mean, u2, h, Zs, sigma_h, sigma_prmean, sigma_prvar)
  h <- aux$Sigtdraw
  h0 <- as.numeric(aux$h0)
  # sqrtvol <- aux$sigt
  # [TODO] fix this
  # sse_2 <- apply( (h[,2:t_max] - h[,1:(t_max-1)])^2, MARGIN = 1, FUN = sum)
  sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
  sigma_post_a <- 1 + rep(t_max,K) # prior of sigma_h Gamma(1,0.0001)
  sigma_post_b <- 0.0001 + sse_2 # prior of sigma_h

  for (i in c(1:K)){
    sigma_h[i,i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
  }

  aux$sigma_h <- sigma_h
  return(aux)
}

#' #' @export
#' lh_post <- function(ytilde, hc, alph, HinvSH_h, t_max, K){
#'   llike = - 0.5*sum(hc) -
#'             0.5*t(ytilde) %*% sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = exp(-hc)) %*% ytilde
#'   lpost = as.numeric(llike +  -.5*Matrix::t(hc-alph)%*%HinvSH_h%*%(hc-alph))
#'   return(lpost)
#' }
#'
#' #' @export
#' sample_h_ele <- function(ytilde, sigma_h = 0.0001*diag(K),
#'                          h = matrix(0, ncol = t_max, nrow = K),
#'                          s0 = rep(0,K), K, t_max){
#'   u = as.numeric(ytilde)
#'   s2 = u^2
#'
#'   Hh = sparseMatrix(i = 1:(t_max*K),
#'                     j = 1:(t_max*K),
#'                     x = rep(1,t_max*K)) -
#'     sparseMatrix( i = (K+1):(t_max*K),
#'                   j = 1:((t_max-1)*K),
#'                   x = rep(1,(t_max-1)*K),
#'                   dims =  c(t_max*K, t_max*K))
#'   SH = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = rep(1./diag(sigma_h), t_max))
#'
#'   HinvSH_h = Matrix::t(Hh) %*% SH %*% Hh
#'   alph = Matrix::solve(Hh, sparseMatrix(i = 1:K, j = rep(1,K), x = s0, dims = c(t_max*K,1)))
#'
#'   e_h = 1
#'   # ht = rep(s0,t_max)
#'   ht = log(s2)
#'
#'   while (e_h > .001){
#'     einvhts2 = exp(-ht)*s2
#'     gh = - HinvSH_h %*% (ht-alph) - 0.5 * (1-einvhts2)
#'     Gh = - HinvSH_h -.5*sparseMatrix(i = 1:(t_max*K),j = 1:(t_max*K), x = einvhts2)
#'     newht = ht - Matrix::solve(Gh,gh)
#'     e_h = max(abs(newht-ht))
#'     ht = newht
#'   }
#'
#'
#'   Kh = -Gh
#'   CKh = Matrix::t(Matrix::chol(Kh))
#'
#'   # c_pri = -t_max*K*0.5*log(2*pi) -.5*t_max*sum(log(sigma_h))
#'   # c_IS = -t_max*K*0.5*log(2*pi) + sum(log(Matrix::diag(CKh)))
#'
#'   logc = lh_post(u, ht, alph, HinvSH_h, t_max, K) + log(3)
#'   # AR-step:
#'   flag = 0
#'   while (flag == 0){
#'     hc = ht + Matrix::solve(Matrix::t(CKh), rnorm(t_max*K))
#'     alpARc =  lh_post(u, hc, alph, HinvSH_h, t_max, K) +
#'               as.numeric(0.5*Matrix::t(hc-ht)%*%Kh%*%(hc-ht)) - logc
#'     if (alpARc > log(runif(1))){
#'       flag = 1
#'     }
#'   }
#'   # MH-step
#'   h <- as.numeric(h)
#'   alpAR = lh_post(u, h, alph, HinvSH_h, t_max, K) +
#'           as.numeric(0.5*Matrix::t(h-ht)%*%Kh%*%(h-ht)) - logc
#'   if (alpAR < 0){
#'     alpMH = 1
#'   } else if (alpARc < 0){
#'     alpMH = - alpAR
#'   } else {
#'     alpMH = alpARc - alpAR
#'   }
#'
#'   if (alpMH > log(runif(1))){
#'     h = hc
#'   }
#'   h <- matrix(h, nrow = K)
#'
#'   # sample h0
#'   V_h0 <- rep(10, K)
#'   mu_h0 <- rep(0, K)
#'   invV_h0post <- diag(1/diag(sigma_h)+1/V_h0)
#'   h0_mean <- solve(invV_h0post) %*% (mu_h0/V_h0 + h[,1]/ diag(sigma_h) )
#'   s0 <- as.numeric(h0_mean) + Matrix::solve(Matrix::t(invV_h0post), rnorm(K))
#'
#'
#'   # sample sigma_h
#'   sse_2 <- apply( (h[,2:t_max] - h[,1:(t_max-1)])^2, MARGIN = 1, FUN = sum)
#'   sigma_post_a <- 1 + rep(t_max,K)
#'   sigma_post_b <- 0.0001 + sse_2
#'
#'   for (i in c(1:K)){
#'     sigma_h[i,i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
#'   }
#'   aux <- list(sigma_h = sigma_h,
#'               Sigtdraw = h,
#'               sigt = exp(h/2),
#'               s0 = s0)
#'
#'   return(aux)
#' }
##########################################################################
# Plot functions  #
##########################################################################

#' @export
plot.fatBVARSV <- function(fatBVARSVobj, element = NULL){
  if (is.null(element)) {
    plot(fatBVARSVobj$mcmc)
  } else {
    plot(get_post.fatBVARSV(fatBVARSVobj, element))
  }
}

##########################################################################
# Get posterior functions  #
##########################################################################
#'@export
get_post <- function(obj, element = NULL, ...) {
  UseMethod("get_post", obj)
}
#' @rdname get_post
#' @export
get_post.fatBVARSV <- function(obj, element = NULL){
  if (is.null(element)) {
    return(obj$mcmc)
  } else {
    mcmc_name <- substr(colnames(obj$mcmc), start = 1, stop = nchar(element))
    mcmc_id = (mcmc_name == element)
    return(obj$mcmc[,mcmc_id])
  }
}

#' @export
get_post.numeric <- function(obj, element = NULL){
  if (is.null(element)) {
    return(obj)
  } else {
    mcmc_name <- substr(names(obj), start = 1, stop = nchar(element))
    mcmc_id = (mcmc_name == element)
    return(obj[mcmc_id])
  }
}

#' @export
get_post.matrix <- function(obj, element = NULL){
  if (is.null(element)) {
    return(obj)
  } else {
    mcmc_name <- substr(colnames(obj), start = 1, stop = nchar(element))
    mcmc_id = (mcmc_name == element)
    return(obj[,mcmc_id])
  }
}
#' @export
get_post.mcmc <- function(obj, element = NULL){
  if (is.null(element)) {
    return(obj)
  } else {
    mcmc_name <- substr(colnames(obj), start = 1, stop = nchar(element))
    mcmc_id = (mcmc_name == element)
    return(obj[,mcmc_id])
  }
}

#' @export
get_lu_param <- function(param){
  if (class(param) == "mcmc") {
    mcmc_name <- colnames(param)
  }

  if (class(param) == "numeric") {
    mcmc_name <- names(param)
  }
  num_param <- length(mcmc_name)
  lb <- rep(-Inf, num_param)
  ub <- rep(Inf, num_param)
  names(lb) <- names(ub) <- mcmc_name

  # B : no restrict
  # A : no restrict

  # sigma : positive restrict
  mcmc_id <- (substr(mcmc_name, start = 1, stop = nchar("sigma")) == "sigma")
  lb[mcmc_id] <- 0

  # w : positive restrict
  mcmc_id <- (substr(mcmc_name, start = 1, stop = nchar("nu")) == "nu")
  lb[mcmc_id] <- 2
  ub[mcmc_id] <- 100

  # w : positive restrict
  mcmc_id <- (substr(mcmc_name, start = 1, stop = nchar("w")) == "w")
  lb[mcmc_id] <- 0

  # gamma : no restrict

  return(list(lb = lb,
              ub = ub))

}

