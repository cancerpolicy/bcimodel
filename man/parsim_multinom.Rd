% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parallelized.R
\name{parsim_multinom}
\alias{parsim_multinom}
\title{Simulate from mutinomial, parallelized}
\usage{
parsim_multinom(nsims, nreps, probs, names, ncores = 4)
}
\arguments{
\item{nsims}{Number of size=1 sims from the multinomials}

\item{nreps}{Number of times to replicate the nsims. Put the smaller number here, the bigger one as nsims}

\item{probs}{Vector of probabilities for each of the K categories}

\item{names}{Vector of names associated with each of the categories. Must be of same length as probs. Can be character or numeric. Use numeric row IDs if you're referring to a group defined by multiple variables whose groupings are defined in a separate data frame}
}
\description{
Simulate matrix of draws from multinomial
}
