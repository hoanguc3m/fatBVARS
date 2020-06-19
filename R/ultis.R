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
    xt <- cbind(xt, rbind(1, vec( as.matrix( t(y0)[,(p+i-1):i]))) )
  }
  for (i in c( p:(t_max-1))){
    xt <- cbind(xt, rbind(1, vec( as.matrix( yt[,i:(i-p+1)]))) )
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
sample_h_ele <- function(ytilde, sigma_h = 0.0001*diag(K), h0_mean = rep(0,K),
                         h = matrix(0, nrow = t_max, ncol = K), K, t_max){
  tmp <- getmix()
  q <- tmp$q
  m_mean <- tmp$m
  u2 <- tmp$u2
  Zs <- matrix(1,t_max,1) %x% diag(K)


  sigma_prmean <- h0_mean # mean h_0
  sigma_prvar <- 1*diag(K)   # variance h_0

  aux <- sigmahelper4(t(ytilde^2), q, m_mean, u2, h, Zs, sigma_h, sigma_prmean, sigma_prvar)
  h <- aux$Sigtdraw
  h0 <- as.numeric(aux$h0)
  # sqrtvol <- aux$sigt
  # [TODO] fix this
  # sse_2 <- apply( (h[,2:t_max] - h[,1:(t_max-1)])^2, MARGIN = 1, FUN = sum)
  if (K>1) {
    sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
  } else {
    sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
  }

  sigma_post_a <- 1 + rep(t_max,K) # prior of sigma_h Gamma(1,0.0001)
  sigma_post_b <- 0.0001 + sse_2 # prior of sigma_h

  for (i in c(1:K)){
    sigma_h[i,i] <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
  }

  aux$sigma_h <- sigma_h
  return(aux)
}

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
    return(obj$mcmc[,mcmc_id, drop = FALSE])
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
    return(as.matrix(obj[,mcmc_id], nrow = nrow(obj)))
  }
}
#' @export
get_post.mcmc <- function(obj, element = NULL){
  if (is.null(element)) {
    return(obj)
  } else {
    mcmc_name <- substr(colnames(obj), start = 1, stop = nchar(element))
    mcmc_id = (mcmc_name == element)
    return(obj[,mcmc_id, drop = FALSE])
  }
}

# get_lu_param <- function(param){
#   if (class(param) == "mcmc") {
#     mcmc_name <- colnames(param)
#   }
#
#   if (class(param) == "numeric") {
#     mcmc_name <- names(param)
#   }
#   num_param <- length(mcmc_name)
#   lb <- rep(-Inf, num_param)
#   ub <- rep(Inf, num_param)
#   names(lb) <- names(ub) <- mcmc_name
#
#   # B : no restrict
#   # A : no restrict
#
#   # sigma : positive restrict
#   mcmc_id <- (substr(mcmc_name, start = 1, stop = nchar("sigma")) == "sigma")
#   lb[mcmc_id] <- 0
#
#   # w : positive restrict
#   mcmc_id <- (substr(mcmc_name, start = 1, stop = nchar("nu")) == "nu")
#   lb[mcmc_id] <- 2
#   ub[mcmc_id] <- 100
#
#   # w : positive restrict
#   mcmc_id <- (substr(mcmc_name, start = 1, stop = nchar("w")) == "w")
#   lb[mcmc_id] <- 0
#
#   # gamma : no restrict
#
#   return(list(lb = lb,
#               ub = ub))
#
# }

#' Compute a two-dimensional kernel density estimate
#'
#' @author Michal Oleszak
#'
#' @description The kernel is assumed to be Gaussian. Bandwidth matrix is
#'              diagonal. The two bandwidth parameters are chosen optimally
#'              without ever using/assuming any parametric model for the data
#'              or any "rules of thumb". Unlike many other procedures, this one
#'              is immune to accuracy failures in the estimation of multimodal
#'              densities with widely separated modes. This function in meant to be
#'              the R implementation of the MATLAB \code{kde2d()} function written
#'              and published by Z. I. Botev at:
#'              \url{http://web.maths.unsw.edu.au/~zdravkobotev/}
#'
#' @references Z. I. Botev, J. F. Grotowski and D. P. Kroese,
#'             "Kernel Density Estimation Via Diffusion",
#'             Annals of Statistics, 2010, Volume 38, Number 5, Pages 2916-2957
#'
#' @param data   N by 2 matrix with the two variables as columns
#' @param n      size of the n by n grid over which the density is computed
#' @param limits limits of the bounding box over which the density is computed;
#'               format: c(lower_Xlim, upper_Xlim, lower_Ylim, upper_Ylim)
#'
#' @return A \code{list} with bandwidth, density and grids for the two dimensions.
#'
#' @importFrom pracma accumarray repmat fzero meshgrid histc numel zeros ones ifft
#' @importFrom stats mvfft
kde2D <- function(data, n = 2 ^ 8, limits = NULL) {
  # define auxiliary functions
  ndhist <- function(data, M) {
    bins <- pracma::zeros(nrow(data), ncol(data))
    for (i in 1:ncol(data)) {
      bins[, i] <- pracma::histc(data[, i], seq(0, 1, 1 / M))$bin
      bins[, i] <- apply(as.matrix(bins[, i]), 1, function(y) min(y, M))
    }
    out <- pracma::accumarray(
      bins[all(bins > 0), ], rep((1 / nrow(data)), nrow(data)),
      M * (pracma::ones(1, ncol(data)))
    )
    return(out)
  }
  dct2d <- function(data) {
    dct1d <- function(x) {
      x <- rbind(x[seq(1, nrow(x), 2), ], x[seq(nrow(x), 2, -2), ])
      out <- Re(weights * mvfft(x))
      return(out)
    }
    w <- as.matrix(c(1, 2 * (exp(-1i * (1:(nrow(data) - 1)) * pi / (2 * nrow(data))))))
    weights <- w[, pracma::ones(1, ncol(data))]
    data <- t(dct1d(t(dct1d(data))))
    return(data)
  }
  idct2d <- function(data) {
    idct1d <- function(x) {
      y <- Re(ifft_mat(weights * x))
      out <- pracma::zeros(nrow(data), ncol(data))
      out[seq(1, nrow(data), 2), ] <- y[1:(nrow(data) / 2), ]
      out[seq(2, nrow(data), 2), ] <- y[seq(nrow(data), nrow(data) / 2 + 1, -1), ]
      return(out)
    }
    ifft_mat <- function(x) {
      x <- apply(x, 2, function(y) pracma::ifft(y))
      return(x)
    }
    w <- as.matrix(exp(1i * (0:(nrow(data) - 1)) * pi / (2 * nrow(data))))
    weights <- w[, pracma::ones(1, ncol(data))]
    out <- idct1d(t(idct1d(data)))
    return(out)
  }
  K <- function(s) {
    if (s == 0) {
      out <- (-1) ^ s / sqrt(2 * pi)
    } else {
      out <- (-1) ^ s * prod(seq(from = 1, to = 2 * s - 1, by = 2)) / sqrt(2 * pi)
    }
    return(out)
  }
  psi <- function(s, time) {
    w <- c((exp(-I * pi ^ 2 * time))) *
      c(cbind(1, 0.5 * pracma::ones(1, length(I) - 1)))
    wx <- t(w * (I ^ s[1]))
    wy <- t(w * (I ^ s[2]))
    out <- (-1) ^ sum(s) * (wy %*% A2 %*% t(wx)) * pi ^ (2 * sum(s))
    return(out)
  }
  func <- function(s, t) {
    if (sum(s) <= 4) {
      sum_func <- func(c(s[1] + 1, s[2]), t) + func(c(s[1], s[2] + 1), t)
      const <- (1 + 1 / 2 ^ (sum(s) + 1)) / 3
      time <- c((-2 * const * K(s[1]) * K(s[2]) / N / sum_func) ^ (1 / (2 + sum(s))))
      out <- psi(s, time)
    } else {
      out <- psi(s, t)
    }
    return(out)
  }
  evolve <- function(t) {
    sum_func <- func(c(0, 2), t) + func(c(2, 0), t) + 2 * func(c(1, 1), t)
    time <- (2 * pi * N * sum_func) ^ (-1 / 3)
    out <- (t - time) / time
    return(out)
  }
  subtract_evolve <- function(t) {
    return(t - evolve(t))
  }

  # compute the kernel density estimate
  N <- nrow(data)
  # round up n to the next power of 2
  n <- 2 ^ ceiling(log2(n))
  # define default limits
  if (is.null(limits)) {
    max <- c(max(data[, 1]), max(data[, 2]))
    min <- c(min(data[, 1]), min(data[, 2]))
    range <- max - min
    max_xy <- max + range / 4
    min_xy <- min - range / 4
  } else {
    max_xy <- c(limits[2], limits[4])
    min_xy <- c(limits[1], limits[3])
  }
  # scale data
  scaling <- max_xy - min_xy
  data_trans <- (data - pracma::repmat(min_xy, N, 1)) / pracma::repmat(scaling, N, 1)
  # bin the data uniformly using regular grid
  data_init <- ndhist(data_trans, n)
  # discrete cosine transform of initial data
  data_DCT <- dct2d(data_init)
  # compute the optimal bandwidth ^ 2
  I <- (0:(n - 1)) ^ 2
  A2 <- data_DCT ^ 2
  t_star <- pracma::fzero(subtract_evolve, c(0, 0.1))[[1]]
  p_02 <- func(c(0, 2), t_star)
  p_20 <- func(c(2, 0), t_star)
  p_11 <- func(c(1, 1), t_star)
  t_y <- (p_02 ^ (3 / 4) / (4 * pi * N * p_20 ^ (3 / 4) * (p_11 + sqrt(p_20 * p_02)))) ^ (1 / 3)
  t_x <- (p_20 ^ (3 / 4) / (4 * pi * N * p_02 ^ (3 / 4) * (p_11 + sqrt(p_20 * p_02)))) ^ (1 / 3)
  # smooth the discrete cosine transform of initial data using t_star
  data_DCT_smoothed <- (t(t(exp(-(0:(n - 1)) ^ 2 * pi ^ 2 * c(t_x) / 2))) %*%
                          t(exp(-(0:(n - 1)) ^ 2 * pi ^ 2 * c(t_y) / 2))) * data_DCT
  # apply the inverse discrete cosine transform
  density <- idct2d(data_DCT_smoothed) *
    (pracma::numel(data_DCT_smoothed) / prod(scaling))
  grid <- pracma::meshgrid(
    seq(min_xy[1], max_xy[1], by = scaling[1] / (n - 1)),
    seq(min_xy[2], max_xy[2], by = scaling[2] / (n - 1))
  )
  # compute the bandwidth
  bw <- sqrt(cbind(t_x, t_y)) * scaling

  return(list(
    "bandwidth" = bw, "density" = density,
    "X" = grid[[1]], "Y" = grid[[2]]
  ))
}
