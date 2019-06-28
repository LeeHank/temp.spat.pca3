% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/temp_spat_pca.R
\name{data_gen_1d_2}
\alias{data_gen_1d_2}
\title{Generate simulation data in thesis p23~24}
\usage{
data_gen_1d_2(s_old_mat, s_new_mat, k, n, lambda_mat, sigma_sq, A, V,
  normalize = FALSE)
}
\arguments{
\item{s_old_mat}{locations for model fitting}

\item{s_new_mat}{locations for spatial prediction}

\item{k}{dim of eigen-space}

\item{n}{time length}

\item{lambda_mat}{state coavariance matrix}

\item{sigma_sq}{measurement error variance}

\item{A}{transition matrix}

\item{V}{Ca}

\item{normalize}{logical value, if TRUE, data matrix will be standardized by column}
}
\value{
a list object contains data matrix, signal, and errors
\describe{
  \item{y_mat}{data matrix (time T by locations N matrix)}
  \item{xsi_mat}{state matrix (time T by eigen-space K matrix)}
  \item{epsilon_mat}{measurement error matrix (time T by locations N matrix)}
  \item{location}{locations for fitting}
  \item{phi_mat}{phi mat}
  \item{location_new}{locations for spatial prediction}
  \item{phi_mat_new}{phi mat new}
}
}
\description{
Generate simulation data in thesis p23~24
}