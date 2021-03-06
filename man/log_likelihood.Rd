% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior.R
\name{log_likelihood}
\alias{log_likelihood}
\title{Log Likelihood of BVAR model}
\usage{
log_likelihood(param, data, aggregate = TRUE)
}
\arguments{
\item{param}{The one sample of the parameters obtained by fatBVARSVobj$mcmc}

\item{data}{The fatBVARSV object from command BVAR.}
}
\value{
The log likelihood density of the fatBVARSVobj at param.
}
\description{
This function returns a log likelihood of BVAR-SV-fatTail model.
}
\examples{
\dontrun{
eval_ll <- log_likelihood(Chain$mcmc[1,], Chain)
}
}
