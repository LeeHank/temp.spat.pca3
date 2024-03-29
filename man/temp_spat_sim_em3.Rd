% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/temp_spat_pca.R
\name{temp_spat_sim_em3}
\alias{temp_spat_sim_em3}
\title{All in one function for simulation experiment use}
\usage{
temp_spat_sim_em3(s, y_mat, y_mat_new, sigma2_eps = 1, itermax = 30,
  tol = 0.001, tau, k, s_new, phi_mat_new, lambda_mat, A, new_xsi_mat)
}
\arguments{
\item{s}{location matrix}

\item{y_mat}{data matrix (training data)}

\item{y_mat_new}{data matrix (testing data)}

\item{sigma2_eps}{measurement error variance}

\item{itermax}{maximum iteration in EM-algorithm}

\item{tol}{stopping criterion}

\item{tau}{roughness parameter}

\item{k}{eigen-space}

\item{s_new}{location matrix for spatial prediction use}

\item{phi_mat_new}{true new phi matrix}

\item{lambda_mat}{true lambda matrix}

\item{A}{state transition matrix}

\item{new_xsi_mat}{true random coefficient matrix}
}
\value{
a list of parameter estimaition, prediction, and performance measure
}
\description{
All in one function for simulation experiment use
}
