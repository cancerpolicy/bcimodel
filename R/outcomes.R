
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

calcmrr <- function(cummortlist, controlindex, reverse=FALSE) {
    return(
           lapply(1:length(cummortlist),
           function(x) {
               r <- cummortlist[[x]]/cummortlist[[controlindex]]
               if (reverse) r <- 1-r
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

compile_outcomes <- function(all, pop_size, futimes, policynames,
                             per=100000) {
    # Each set of outcomes has to be summarized across sims and 
    # multipled by the "per" denominator
    standardized <- lapply(all, function(x) { 
                               if (is.list(x)) {
                                lapply(x, function(y) { colMeans(y)*per/pop_size })
                               } else return(colMeans(x)*per/pop_size)
    })
    condense1 <- lapply(standardized, function(x) {
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
    byfu <- vector('list', length=length(futimes))
    for (i in 1:length(futimes)) {
        byfu[[i]] <- condenseall[which(rownames(condenseall)==
                                       as.character(futimes[i])),]
        rownames(byfu[[i]]) <- names(all)
        colnames(byfu[[i]]) <- policynames
        byfu[[i]] <- round(byfu[[i]],2)
    }
    names(byfu) <- futimes

    return(byfu)

}

