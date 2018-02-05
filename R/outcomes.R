
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
#' @param ratematrix Matrix of exponential rates to simulate from
#' 
#' @examples
#' # Cell-specific rates: two columns of 0.05, two columns of 0.2
#' m <- cbind(matrix(0.05, nrow=100, ncol=2),
#'          matrix(0.2, nrow=100, ncol=2))
#' # Sims reflect cell-specific rates
#' round(1/colMeans(rexp_matrix(m)), 2)
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
#'  strategies with early detection
#' @param A matrix of pop_size by nsim; cells are treatment-specific mortality rates for each individual in each sim
#' @param pop_size Population size
#' @param nsim Number of sims
#'
#' @examples
#' See simpolicies.R
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
#' For stage-shifted cases in an early detection scenario, update time to cancer death
#'
#' Takes treatments from a paired scenario with no early detection and
#' updates times to cancer death only for those advanced-stage cases who were
#' stage-shifted in the early detection scenario. Updated times are quantile-
#' correlated to the previous times
#' 
#' @param policies Data frame specifying policies - see ex1$pol. If 'pairnum' = NA, no treatments are updated. If 'pairnum' is numeric, treatments from that number scenario will be the starting point, with early-detected cases shifted according to the next input
#' @param shifts List of treatment-shift indicators (see shifttreatment_indicator)
#' @param rates Matrix of treatment-specific rates, size popsize by nsim
#' @param rates Matrix of times to cancer death under no stage-shifts, size popsize by nsim
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
# Summary functions
#-------------------------------------------------------------------------------
#' Tally cumulative incidence of cancer
#'
#' Count incident cancer cases within a follow-up period
#' 
#' @param futimes Vector of follow-up times to assess
#' @param times Matrix of times to cancer incidence
#'
#' @return Matrix where rows tally incident cases for each sim; columns represent the various follow-up times
#' @examples
#' # Assess 5- and 10-year follow-ups
#' fu <- c(5, 10)
#' times <- matrix(rexp(n=25, rate=1/5), nrow=5, ncol=5)
#' # Rows are simulations, columns are 5- and 10-year follow-up times, 
#' # and cells are incidence counts within the periods
#' cuminc(fu, times)
#' 
#' @return Matrix of incidence counts for sims (rows) and follow-up times (columns)
#' @export

cuminc <- function(futimes, times) {
    sapply(futimes, function(x) colSums(times<=x))
}

#' Tally cumulative incidence of cancer mortality
#'
#' Count cancer mortality within a follow-up period
#' 
#' @param futimes Vector of follow-up times to assess
#' @param times Matrix of times to cancer mortality
#' @param cancerdeath Logical matrix indicating death from cancer
#'
#' @examples
#' 
#' @return Matrix of mortality counts for sims (rows) and follow-up times (columns)
#' @export
cummort <- function(futimes, times, cancerdeath) {
    sapply(futimes, function(x) colSums(times<=x & cancerdeath))
}

#' Tally years lived within a follow-up period
#'
#' Count years alive within a follow-up period, summed over all individuals per sim
#' 
#' @param futimes Vector of follow-up times to assess
#' @param times Matrix of times to cancer mortality
#'
#' @examples
#' # Assess 5- and 10-year follow-ups
#' fu <- c(5, 10)
#' times <- matrix(rexp(n=25, rate=1/5), nrow=5, ncol=5)
#' cumyears(fu, times)
#' 
#' @return Matrix of cumulative years lived within the follow-up period for sims (rows) and follow-up times (columns)
#' @export
cumyears <- function(futimes, times) {
    sapply(futimes, function(x) colSums(ifelse(times<=x, times, x)))
}

#' Calculate mortality rate ratios
#'
#' For a list of scenarios, calculate mortality rate ratios (MRRs) compared to the control scenario
#' 
#' @param cummortlist List of cummort() outputs
#' @param controlindex Number indicating which cummortlist element is the control scenario
#' @param reverse Set to TRUE to return 1-MRR
#' @param perc Set to TRUE to report a percent, instead of a ratio
#'
#' @examples
#' 
#' @return Matrix of MRRs for each scenario/cummortlist element
#' @export
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

#' Calculate absolute risk reduction
#'
#' For a list of scenarios, calculate absolute risk reductions (ARRs) compared to the control scenario
#' 
#' @param cummortlist List of cummort() outputs
#' @param controlindex Number indicating which cummortlist element is the control scenario
#' @param reverse Set to TRUE to return -1*ARR
#'
#' @examples
#' 
#' @return Matrix of ARRs for each scenario/cummortlist element
#' @export
calcarr <- function(cummortlist, controlindex, reverse=FALSE) {
    return(
           lapply(1:length(cummortlist),
           function(x) {
               diff <- cummortlist[[controlindex]]-cummortlist[[x]]
               if (reverse) return(-1*diff) else return(diff)
           })
    )
}

#' Compute mean and/or quantiles
#' Return mean, lower, or upper quantile. Default quantile type is 7, which gives a continuous result
#' @param matrix Matrix of values to summarize
#' @param stat Statistic to return: 'mean', 'lower' or 'upper'
#' @param quanttype Quantile type - see R help page for quantile()
#' 
#' @export
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
#' @export
round_matrix <- function(mat, digbyrow) {
    newmat <- t(sapply(1:nrow(mat), function(x) {
        as.character(round(mat[x,], digbyrow[x]))
        }))
    rownames(newmat) <- rownames(mat)
    colnames(newmat) <- colnames(mat)
    return(newmat)
}

#' Format bounds for reporting
#'
#' Lower and upper are tables of the same dimensions with rownames. For including a meanmat, the line break "\\n" between mean and uncertainty works for exporting as csv (it will get stripped) or for printing with pander but NOT printing with kable.
#' @export
format_bounds <- function(lower, upper, digits=NULL, paren=FALSE,
                          meanmat=NULL) {
    # If digits is not null, round
    if (!is.null(digits)) {
        if (length(digits)!=nrow(lower)) stop('In format_bounds, digits must
                                      have length of nrow(lower)')
        lower <- round_matrix(lower, digits)
        upper <- round_matrix(upper, digits)
        if (!is.null(meanmat)) meanmat <- round_matrix(meanmat, digits)
    }
    bounds <- sapply(1:ncol(lower), function(c) {
                beginning <- ifelse(paren, '(', '')
                ending <- ifelse(paren, ')', '')
                if (is.null(meanmat)) {
                    paste0(beginning, lower[,c], ', ', upper[,c], ending)
                } else {
                    paste0(meanmat[,c], ' ',
                           beginning, lower[,c], ', ', upper[,c], ending)
                }
              })
    rownames(bounds) <- rownames(lower)
    colnames(bounds) <- colnames(lower)
    # An aside - replace NA's with blanks
    replacewithblank <- bounds=='NA' | bounds=='NA (NA, NA)'
    bounds[replacewithblank] <- ''
    return(bounds)
}

#' Return bounds for each element of the list
#' Options are: lower, upper; (lower, upper); mean (lower, upper)
#' @export
format_bounds_list <- function(thislist, digits=NULL, paren=FALSE,
                               includemean=FALSE, compileall=FALSE) {
    if (!includemean) {
        flist <- lapply(thislist, function(x) {
                   format_bounds(x$lower, x$upper, digits, paren)
        })
    } else {
        flist <- lapply(thislist, function(x) {
                   format_bounds(x$lower, x$upper, digits, paren,
                                 meanmat=x$mean)
        })
    }
    if (compileall) {
        for (i in 1:length(flist)) {
            flist[[i]] <- as.data.frame(flist[[i]], check.names=FALSE)
            flist[[i]]$Measure <- rownames(flist[[i]])
        }
        ftable <- ldply(flist)
        ftable <- rename(ftable, c('.id'='Follow-Up Year'))
        addblanks <- rep('', nrow(ftable))
        keepthese <- seq(1,nrow(ftable),by=nrow(flist[[1]]))
        addblanks[keepthese] <- ftable[['Follow-Up Year']][keepthese]
        ftable[['Follow-Up Year']] <- addblanks
        statcol <- which(colnames(ftable)=='Measure')
        ftable <- ftable[,c(1,statcol, 2:(statcol-1))]
    } else ftable <- flist

    return(ftable)
}

#' Compile results long
#' @export
compile_long <- function(thislist) {
    for (i in 1:length(thislist)) {
        for (j in 1:length(thislist[[i]])) {
        thislist[[i]][[j]] <- 
            as.data.frame(thislist[[i]][[j]], check.names=FALSE)
        thislist[[i]][[j]]$Measure <- rownames(thislist[[i]][[j]])
        }
    }
    thislist2 <- lapply(thislist, ldply, .id='Statistic')
    df <- ldply(thislist2, .id='Year')
    df <- melt(df, id.vars=c('Year', 'Statistic', 'Measure'), 
               variable_name='Scenario')
    df <- cast(df, Year+Scenario+Measure~Statistic)
    return(df)
}

#' Plot results
#' Select on measures and panel them. Group on f-u year?
#' @export
plot_results <- function(rlong, measures=NULL, type='bar',
                         rotate_xlab=TRUE) {
    if (is.null(measures)) measures <- unique(rlong$Measure)
    rlong <- subset(rlong, Measure %in% measures)
    switch(type,

    'bar'= {
        # I don't love the bar
        g <- ggplot(rlong, aes(x=Year, y=mean, fill=Scenario)) +
            geom_bar(position=position_dodge(), stat='identity',
                     colour='black', size=0.3) +
            geom_errorbar(aes(ymin=lower, ymax=upper),
                          size=0.3, width=0.2, 
                          position=position_dodge(0.9)) + 
            facet_grid(Measure~., scales='free_y') +
            xlab('Follow-Up Year') + ylab('') +
            scale_fill_hue(name='Scenario') + 
            theme_bw()

    },

    'line'={

        # Line
        pd <- position_dodge(0.1)
        g <- ggplot(rlong, aes(x=Scenario, y=mean, color=Year, group=Year)) +
            geom_errorbar(aes(ymin=lower, ymax=upper), colour='black',
                          position=pd, width=0.1) +
            geom_line(position=pd) +
            geom_point(position=pd, size=2, shape=21, fill='white') +
            facet_grid(Measure~., scales='free_y') +
            xlab('Scenario') + ylab('') +
            scale_fill_hue(name='Follow-Up Year') + 
            expand_limits(y=0) +
            theme_bw()
        if (rotate_xlab) g <- g+ theme(axis.text.x=element_text(angle=90, 
                                                                hjust=1))
        
    })
return(g)
}

