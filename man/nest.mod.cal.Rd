% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nesting_atmos.R
\name{nest.mod.cal}
\alias{nest.mod.cal}
\title{Function: Nested bias correction for calibration - atmospheric variables}
\usage{
nest.mod.cal(n, y, obs.sta, start.yr, end.yr, fun = "mean", tol = 0.1)
}
\arguments{
\item{n}{The number of columns (e.g., stations)}

\item{y}{A matrix of time series atmospheric variables}

\item{obs.sta}{Observed (Target) data statistics}

\item{start.yr}{Start year}

\item{end.yr}{End year}

\item{fun}{A function to compute the summary statistics, "mean" or "sum".}

\item{tol}{Tolerance for small values, default 0.1}
}
\value{

}
\description{
Function: Nested bias correction for calibration - atmospheric variables
}
\examples{

}
