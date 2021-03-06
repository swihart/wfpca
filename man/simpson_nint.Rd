% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/simpson_nint.R
\name{simpson_nint}
\alias{simpson_nint}
\title{numerical integration based on Simpsons rule}
\usage{
simpson_nint(x, fx, n.pts = max(256, length(n.x)))
}
\arguments{
\item{x}{a sequence of x values.}

\item{fx}{the value of the function to be integrated at x or a function}

\item{n.pts}{the number of points to be used in the integration. If x contains more than n.pts then n.pts will be set to length(x)}
}
\value{
results a data frame with rmse for each approach for that simulated dataset
}
\description{
based on Bolstad:sintegral(); but needed to remove the message.
}
\examples{
---
}

