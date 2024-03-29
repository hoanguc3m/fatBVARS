% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simBvarNoVol.R
\name{sim.VAR.novol}
\alias{sim.VAR.novol}
\title{Simulate data from Gaussian BVAR model without SV}
\usage{
sim.VAR.novol(
  dist,
  K = 5,
  p = 2,
  t_max = 1000,
  b0 = 0.5,
  a0 = 0.5,
  h = 0,
  nu = 6,
  gamma = 0.5,
  y0 = matrix(0, ncol = K, nrow = p),
  sigma_G = NULL,
  seednum = 0,
  burn_in = 0
)
}
\arguments{
\item{dist}{The variable specifies the BVAR error distribution. It should be one of
c("Gaussian","Student","Skew.Student","Skew.Student", "MT","Skew.MT","MST").}

\item{K}{The number of variables in BVAR model.}

\item{p}{The number of lags in BVAR model.}

\item{t_max}{The number of observations}

\item{b0}{The coefficients of B matrix.}

\item{a0}{The coefficients of A matrix.}

\item{h}{The log diag(Sigma) matrix.}

\item{nu}{The degree of freedom.}

\item{gamma}{The skewness from \eqn{[-gamma, gamma]}.}

\item{seednum}{The default seed.}

\item{burn_in}{The discarded observations.}
}
\value{
A list of simulated data with its specification.
}
\description{
This function simulate data from a BVAR model without SV.
\deqn{y_t = B x_t + SQRT(w_t) A^(-1) Sigma eps_t}
}
\examples{
\dontrun{
datagen <- sim.VAR.novol(dist="Gaussian")
y <- datagen$y
}
}
