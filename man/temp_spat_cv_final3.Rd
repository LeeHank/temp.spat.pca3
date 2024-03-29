% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/temp_spat_pca.R
\name{temp_spat_cv_final3}
\alias{temp_spat_cv_final3}
\title{Tuning process of Temporal Spatial PCA}
\usage{
temp_spat_cv_final3(s, y_mat, sigma2_eps = 1, itermax = 30,
  tol = 0.001, tau, k)
}
\arguments{
\item{s}{location matrix}

\item{y_mat}{data matrix}

\item{sigma2_eps}{measurement error variance}

\item{itermax}{maximum iteration of EM-algorithm}

\item{tol}{stopping criterion of EM}

\item{tau}{a vector contains candidates of roughness parameter}

\item{k}{a vector contains candidates of eigen-space}
}
\value{
a list of tuning process
}
\description{
Tuning process of Temporal Spatial PCA
}
\examples{
\dontrun{
s = sim_data$location
y_mat = sim_data$y_mat[1:100,]
tau = 2^seq(-10,10,by=1)
k= c(1,2,3,4)

library(foreach)
library(doParallel)
cl = makeCluster(6)
registerDoParallel(cl)

cv_result_all = temp_spat_cv_final3(s=s, y_mat=y_mat, tau=tau, k=k)
}
}
