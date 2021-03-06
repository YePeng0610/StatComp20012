% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatCompR.R
\name{sigmA}
\alias{sigmA}
\title{Estimation of Confidence interval of population variance}
\usage{
sigmA(alpha, x, mu = NA)
}
\arguments{
\item{alpha}{Significance level}

\item{x}{A set of data from a normal population}

\item{mu}{mean of population distribution}
}
\value{
Two critical values and an interval
}
\description{
Estimate the confidence interval of variance from the normal population(The population mean is known)
}
\examples{
\dontrun{
x<-rnorm(10,mean=0,sd=1)
m<-sigmA(alpha=0.05,x,mu=0)
m
}
}
