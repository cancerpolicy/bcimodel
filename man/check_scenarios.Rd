% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simpolicies.R
\name{check_scenarios}
\alias{check_scenarios}
\title{Check accurate specification of 'scenarios' parameter to simpolicies()}
\usage{
check_scenarios(scenarios, treatinfo)
}
\arguments{
\item{scenarios}{Data frame of scenarios to simulate (see example below). First scenario should be the base case.}

\item{treatinfo}{Data frame with treatment hazard ratios, and for each scenario, treatment proportions, by stage-subgroups. See the example and the vignette}
}
\value{
NULL if scenarios are paired properly, or vector listing the problematic pairs
}
\description{
The 'scenarios' parameter must contain paired policies, i.e. if one scenario has early detection, a scenario must exist that has the same treatments but no early detection. This function checks for the correct pairing
}
\examples{
data(ex1)
# See scenarios
ex1$pol
# See treatment
ex1$tx
# Create version that have the "tam" scenario deleted and specify
# "base" and "tamandshift" scenarios as paired. This should fail.
fail <- ex1
fail$pol <- fail$pol[-2,]
fail$pol$num[2] <- 2
fail$pol$pairnum[2] <- 1
fail$tx$tam <- NULL

check_scenarios(ex1$pol, ex1$tx)
check_scenarios(fail$pol, fail$tx)
}
