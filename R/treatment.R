
################################################################################
# Functions to simulate treatments received
################################################################################

return_value_from_id = function# Use an ID no to look up a value

##description<< Using a matrix of IDs, look up a value in a dataframe

(ids,
    ### Matrix of ids
df,
    ### Data frame containing the values of interest, sorted
    ### by ID (i.e. row #s == ID #s)
value
    ### Column name of value to return from the df
) {

    return(apply(ids, 2, function(x) df[,value][x]))
    
}

#-------------------------------------------------------------------------------
# shifttreatment_indicator
#-------------------------------------------------------------------------------
#' Create indicator of whether treatment should be re-simulated
#' 
#' Determines TRUE or FALSE indicator based on scenario settings.
#' 
#' @param x Numeric ID
#' @param type Vector of either NA or numeric ID referring to scenario number that is
#  the scenario with the same treatment without any early detection
#' @param shifts List of shift matrices
#' @param basecase Base case matrix of stage-subgroup indicators
#' @param map Map indicating stage-subgroup pairs
#'
#' @examples
# # Port over from test
#'
#' @export

shifttreatment_indicator <- function(x, type, shifts, basecase, map) {

    t <- type[x]
    # NA means there's no early detection so we don't do stage-shifting
    # and subsequent treatment-shifting
    if (is.na(t)) {
        ind <- matrix(NA, nrow=nrow(basecase), ncol=ncol(basecase))
    } else {
        shift <- shifts[[x]]
        ind <- (shift==1 & basecase%in%map['Advanced',])
    }
    return(ind)
}

#-------------------------------------------------------------------------------
# sim_treatment_by_subgroup
#-------------------------------------------------------------------------------
#' Simulate treatment by subgroups

#' Given the proportions of treatments by subgroups, simulate treatment
#' @param treat_chars Data frame indicating stage-subgroup #s in the SSno 
#' variable and stage-subgroup-treatment #s in txSSno, and 
#' prop_* columns where * is the trial name. 
#' @param stage_subgroup_rows Matrix indicating, for each person-sim, 
#' which stage-subgroup
#' they are in (using the SSid from treat_chars)
#' @param thisprop Character indicating the prop_* column name from which
#' to take treatment proportions
#' @param pop_size Population size
#' @param nsim Number of sims
#' 
#' @export

sim_treatment_by_subgroup = function(treat_chars, stage_subgroup_rows, thisprop, 
                                     pop_size, nsim) {
    
    # Matrix for results
    results <- matrix(NA, ncol=nsim, nrow=pop_size)

    if (is.factor(thisprop)) stop('Need character for thisprop, not factor')
    # Loop through stage-subgroups
    for (g in sort(unique(treat_chars$SSno))) { 
        # Subset to this stage-subgroup (SS)
        tempid <- subset(treat_chars, SSno==g)        
        # Simulate treatment for everyone
        # using this group's probs
        treat_choice <- 
            sim_multinom(nsims=pop_size, 
                         nreps=nsim,
                         probs=tempid[,thisprop],
                         names=tempid[,'txSSno'])
        # Assign values to those in group g
        results[stage_subgroup_rows==g] <-
            treat_choice[stage_subgroup_rows==g]
    }

    return(results)
}

#-------------------------------------------------------------------------------
# treatments_by_policy
#-------------------------------------------------------------------------------
#' Use sim_treatment_by_subgroup to simulate treatments
#' 
#' Simulate treatments according to specified policy rules
#' 
#' @param policies A "scenarios" data frame containing an 'id' for the policies
#' and a 'pairnum' column indicating either NA or the paired policy, for
#  strategies with early detection
#'
#' @export

treatments_by_policy <- function(policies, treat_chars, stagegroups,
                                 map, pop_size, nsim) {

    # Two possible policy types: if policies$pairnum is NA, then 
    # sim treatments for everyone. If it's not, that indicates that
    # we need to sim only early-stage treatments to prepare for
    # stage-shifting.
    treats <- lapply(policies$id, 
                     function(x) {
                         policy <- policies$num[which(policies$id==x)]
                         pairnum <- policies$pairnum[which(policies$id==x)]
                         if (is.na(policies$pairnum[which(policies$id==x)])) {
                             t <- sim_treatment_by_subgroup(treat_chars, 
                                                            stagegroups[[1]], 
                                                            x, pop_size, nsim)
                         } else {
                             # This sims treatments only for the early
                             # stages - leaves NA's for the advanced stages
                             # Use the stagegroup matrix that reflects the
                             # stages AFTER stage-shifting
                             earlytreat <- subset(treat_chars, 
                                                  SSno%in%map['Early',])
                             t <- sim_treatment_by_subgroup(earlytreat, 
                                                            stagegroups[[policy]], 
                                                            x, pop_size, nsim)
                         }
                         return(t)
                     })
    return(treats)

}

#-------------------------------------------------------------------------------
# update_treat_stageshift 
#-------------------------------------------------------------------------------
#' For stage-shifted cases in an early detection scenario, update treatment
#' 
#' Takes treatments from a paired scenario with no early detection and 
#' updates treatments only for those advanced-stage cases who were
#' stage-shifted in the early detection scenario
#'
#' @export

update_treat_stageshift <- function(policies, shifts, treats) {
    for (i in policies$num) {
        pairnum <- policies$pairnum[i]
        if (!is.na(pairnum)) { 
            # If there's a paired scenario, start with treatments from there
            temp <- treats[[pairnum]]
            # Now, for stage-shifted cases, insert new early-stage treatments
            # Remember, these were simulated keeping the subgroup constant
            temp[shifts[[i]]] <- treats[[i]][shifts[[i]]]
            # Update treats
            treats[[i]] <- temp
        }
    }
    return(treats)
}
                                    




