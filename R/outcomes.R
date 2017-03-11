
################################################################################
# Functions to simulate and tabulate outcomes
################################################################################

#-------------------------------------------------------------------------------
# rexp_matrix
#-------------------------------------------------------------------------------
#' Sim from exponential given rates in matrix form
#'
#' This is simple code but it's a separate function to confirm that
#' it performs as expected
#'
#' @export

rexp_matrix <- function(ratematrix) {
    return(
           apply(ratematrix, 2, 
                 FUN=function(x) { rexp(rep(1, length(x)), rate=x) })
    )

}

#-------------------------------------------------------------------------------
# timetocancerdeath_by_policy
#-------------------------------------------------------------------------------
#' Sim time to cancer death based on policy guidelines
#' 
#' Sim non-stage-shifted times, to be updated for shifted cases using the 
#' update_time_stageshift function
#' 
#' @param policies A "scenarios" data frame containing an 'id' for the policies
#' and a 'pairnum' column indicating either NA or the paired policy, for
#  strategies with early detection
#'
#' @export

timetocancerdeath_by_policy <- function(policies, rates, pop_size, nsim) {

    # Two possible policy types: if policies$pairnum is NA, then 
    # sim times for everyone. If it's not, that indicates that
    # we will handle the paired sim for stage-shifted cases in 
    # update_time_stageshift 
    times <- lapply(policies$id, 
                     function(x) {
                         policy <- policies$num[which(policies$id==x)]
                         pairnum <- policies$pairnum[which(policies$id==x)]
                         if (is.na(pairnum)) {
                             t <- rexp_matrix(rates[[policy]])
                         } else {
                            t <- matrix(NA, nrow=pop_size, ncol=nsim)
                         }
                         return(t)
                     })
    return(times)

}

#-------------------------------------------------------------------------------
# update_time_stageshift 
#-------------------------------------------------------------------------------
#' For stage-shifted cases in an early detection scenario, update treatment
#' 
#' Takes treatments from a paired scenario with no early detection and 
#' updates treatments only for those advanced-stage cases who were
#' stage-shifted in the early detection scenario
#'
#' @export

update_time_stageshift <- function(policies, shifts, rates, times) {
    for (i in policies$num) {
        pairnum <- policies$pairnum[i]
        if (!is.na(pairnum)) { 
            # Start with times from paired scenario
            temp <- times[[pairnum]]
            # Simulate quantile-correlated times
            correlated <- sim_same_qexp(oldtime=temp, 
                                        oldrate=rates[[pairnum]],
                                        newrate=rates[[i]],
                                        prefix='c')
            # For stage-shifted cases, insert quantile-correlated
            # new times
            temp[shifts[[i]]] <- correlated[shifts[[i]]]
            # Update times
            times[[i]] <- temp
        }
    }
    return(times)
}

#-------------------------------------------------------------------------------
# summarize
#-------------------------------------------------------------------------------
                                    
cuminc <- function(futimes, times) {
    sapply(futimes, function(x) colSums(times<=x))
}

cummort <- function(futimes, times, cancerdeath) {
    sapply(futimes, function(x) colSums(times<=x & cancerdeath))
}

cumyears <- function(futimes, times) {
    sapply(futimes, function(x) colSums(ifelse(times<=x, times, x)))
}

calcmrr <- function(cummortlist, controlindex, reverse=FALSE, perc=FALSE) {
    return(
           lapply(1:length(cummortlist),
           function(x) {
               r <- cummortlist[[x]]/cummortlist[[controlindex]]
               if (reverse) r <- 1-r
               if (perc) r <- 100*r
               return(r)
           })
    )
}

calcarr <- function(cummortlist, controlindex, reverse=FALSE) {
    return(
           lapply(1:length(cummortlist),
           function(x) {
               diff <- cummortlist[[controlindex]]-cummortlist[[x]]
               if (reverse) return(-1*diff) else return(diff)
           })
    )
}

# Small function to return mean, lower, or upper quantile
# Default quantile type is 7, which gives a continuous result
# Choose type
returnstat4matrix <- function(matrix, stat, quanttype=7) {
    if (!quanttype%in%c(1,7)) stop('Choose quantile type 1 or 7')
    if (!stat%in%c('mean', 'lower', 'upper')) stop('stat not supported')
    switch(stat, 
           mean=colMeans(matrix),
           lower=apply(matrix, 2, quantile, na.rm=TRUE, probs=0.025, type=quanttype),
           upper=apply(matrix, 2, quantile, na.rm=TRUE, probs=0.975, type=quanttype))
}

# Create documentation using the test from test_outcomes.R
# "all" is a list where each element besides Cumulative Incidence is also a list,
# one element for each policy. Within each element, rows for sims and cols for futimes.
compile_outcomes <- function(all, pop_size, futimes, policynames,
                             stats=c('mean'), # can add 'lower' and 'upper'
                             per=100000) {

    # Each set of outcomes has to be summarized across sims and 
    # multipled by the "per" denominator
    standardized <- sapply(stats, function(stat) {
                                lapply(all, function(x, s=stat) { 
                                           if (is.list(x)) {
                                                lapply(x, function(y, ss=s) { 
                                                    returnstat4matrix(y, ss)*per/pop_size 
                                                })
                                           } else return(returnstat4matrix(x, s)*per/pop_size)
                                })
                           }, USE.NAMES=TRUE, simplify=FALSE)
    # Rbind futimes, separately for each stat. Columns are policies.
    # This probably needs to be fiddled with for years of life saved, since
    # the sign makes it confusing?
    condensed <- sapply(stats, function(stat) {
                                lapply(standardized[[stat]], function(x) {
                                        if (is.list(x)) {
                                            r <- t(do.call('rbind', x))
                                        } else {
                                            r <- cbind(matrix(x, nrow=length(x)), 
                                                  matrix(NA, nrow=length(x), 
                                                         ncol=length(policynames)-1))
                                        }
                                        rownames(r) <- futimes
                                        return(r)
                                })
                           }, USE.NAMES=TRUE, simplify=FALSE)
    condenseall <- sapply(stats, 
                           function(stat) {
                                do.call(rbind, condensed[[stat]])
                           }, USE.NAMES=TRUE, simplify=FALSE)

    # Previously just means - can delete this later
    if (1==0) {
        standardized <- lapply(all, function(x) { 
                                   if (is.list(x)) {
                                    lapply(x, function(y) { colMeans(y)*per/pop_size })
                                   } else return(colMeans(x)*per/pop_size)
        })
        condensed <- lapply(standardized, function(x) {
                                if (is.list(x)) {
                                    r <- t(do.call('rbind', x))
                                } else {
                                    r <- cbind(matrix(x, nrow=length(x)), 
                                          matrix(NA, nrow=length(x), 
                                                 ncol=length(standardized[[5]])-1))
                                }
                                rownames(r) <- futimes
                                return(r)
                            })
        condenseall <- do.call(rbind, condense1)
    }

    # Compile into tables separate for each follow-up time
    byfu <- lapply(futimes, function(fu) {
        sapply(stats, function(stat, i=fu) {
                   table <- condenseall[[stat]][
                                which(rownames(condenseall[[stat]])== 
                                      as.character(i)),]
                            rownames(table) <- names(all)
                            colnames(table) <- policynames
                            table <- round(table, 2)
                            return(table)
        }, USE.NAMES=TRUE, simplify=FALSE)
    })
    names(byfu) <- futimes
    return(byfu)

}

# Round a matrix, specifying digits per row
round_matrix <- function(mat, digbyrow) {
    newmat <- t(sapply(1:nrow(mat), function(x) {
        as.character(round(mat[x,], digbyrow[x]))
        }))
    rownames(newmat) <- rownames(mat)
    return(newmat)
}

# Lower and upper are tables of the same dimensions with rownames
format_bounds <- function(lower, upper, digits=NULL) {
    # If digits is not null, round
    lower <- round_matrix(lower, digits)
    upper <- round_matrix(upper, digits)
    bounds <- sapply(1:ncol(lower), function(c) { 
                paste(lower[,c], upper[,c], sep=', ') 
              })
    rownames(bounds) <- rownames(lower)
    return(bounds)
}

# Return bounds for each element of the list
format_bounds_list <- function(list, digits=NULL) {
    lapply(list, function(l) format_bounds(l$lower, l$upper, digits))
}


