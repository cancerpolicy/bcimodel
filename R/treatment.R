
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
#' After treatment has been simulated for a base scenario, this function determines whether treatment should be re-simulated due to early detection. Only stage-shifted cases need a new treatment simulated. 
#' 
#' @param x Numeric ID indicating which scenario
#' @param type Vector, one element for each scenario. NA indicates no early detection. Numeric ID implies early detection, and the number refers to the scenario number that has the same treatment but NO early detection
#' @param shifts List of shift matrices (see shifttreatment_indicator)
#' @param basecase Base case matrix of stage-subgroup indicators
#' @param map Map indicating stage-subgroup pairs
#'
#' @examples
#' # Use ex1. ex1$pol shows two scenarios with no early detection (1 and 2), but the 3rd has a 30% stageshift (see earlydetHR). The "pairnum" variable indicators that #3 has the same treatments as #2.
#' library(bcimodel)
#' data(ex1)
#' # Create stageshift indicator matrices for all 3 scenarios: no stage shifts for #1 and #2, but 30% stageshift for #3. Use a small population of size 20, and 2 sims
#' stageshifts <- list(base=matrix(0, nrow=20, ncol=2),
#'                     tam=matrix(0, nrow=20, ncol=2), 
#'                     tamshift=stageshift_indicator(0.85, 20, 2))
#' # ex1$nh shows that there are 4 stage-subgroups. Use a fake random distribution of groups 1:4 for the population before stage-shifting.
#' popdistr <- matrix(sample.int(4, size=40, replace=TRUE), nrow=20, ncol=2)
#' # Now create shifttreatment_indicators for each scenario. We expect no treatment 
#' # shifting for 1-2, since there's no early detection.
#' shifttreatment_indicator(1, ex1$pol$pairnum, stageshifts, popdistr, ex1$map)
#' shifttreatment_indicator(2, ex1$pol$pairnum, stageshifts, popdistr, ex1$map)
#' shifttreatment_indicator(3, ex1$pol$pairnum, stageshifts, popdistr, ex1$map)
#'
#' @return Matrix of either NAs if no early detection-induced treatment-shifting, or matrix of TRUE/FALSE where TRUE indicates a need for treatment shifting
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
#' prop_* columns where * is the trial name. See ex1$tx
#' @param stage_subgroup_rows Matrix indicating, for each person-sim, 
#' which stage-subgroup they are in (corresponding to the SSid from treat_chars)
#' @param thisprop Character indicating the column name from which
#' to take treatment proportions
#' @param pop_size Population size
#' @param nsim Number of sims
#' 
#' @examples
#' # ex1$nh shows that there are 4 stage-subgroups. Use a fake random distribution of groups 1:4 for the population before stage-shifting. Use population of size 20 and 2 sims
#' popdistr <- matrix(sample.int(4, size=40, replace=TRUE), nrow=20, ncol=2)
#' # Simulate treatments for the 'base' scenario of ex1
#' treats <- sim_treatment_by_subgroup(ex1$tx, popdistr, 'base', 20, 2)
#' # Proportions reflect the treatment distribution and the stage-subgroup distribution
#' prop.table(table(treats))
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
#' and a 'pairnum' column indicating either NA or the paired policy ID, for strategies with early detection. See ex1$pol
#' @param treat_chars Data frame with "txSSno" column indicating treatment numbers and subsequent columns with treatment proportions WITHIN stage-subgroups. Each of these columns should correspond to a row in the "policies" data frame, with their names taken fro policies$id. See ex1$tx
#' @param stagegroups List of stage-subgroup matrices, one for each policy/row in the "scenarios" data frame
#' @param map Stage-subgroup map indicating allowed stage-shifts. See ex1$map. 
#' @param pop_size Population size (number of rows)
#' @param nsim Number of sims (number of columns)
#'
#' @return List of treatment matrices, one for each policy in the "scenarios" data frame. Each matrix contains treatment IDs corresponding to treat_chars$txSSno. Early detection scenarios will have NAs for advanced-stage cases who aren't stage-shifted.
#' 
#' @examples
#' library(bcimodel)
#' data(ex1) 
#' # ex1$nh shows that there are 4 stage-subgroups. Use a fake random distribution of groups 1:4 for the population before stage-shifting.
#' popdistr <- matrix(sample.int(4, size=40, replace=TRUE), nrow=20, ncol=2)
#' # Create stageshift indicator matrices for all 3 scenarios: no stage shifts for #1 and #2, but 30% stageshift for #3. Use a small population of size 20, and 2 sims
#' stageshifts <- list(base=matrix(0, nrow=20, ncol=2),
#'                     tam=matrix(0, nrow=20, ncol=2), 
#'                     tamshift=stageshift_indicator(0.85, 20, 2))
#' # Get the actual stages - only policy #3 has stage-shifting
#' stages <- lapply(stageshifts, shift_stages, original=popdistr, map=ex1$map)
#' lapply(stages, table)
#' t <- treatments_by_policy(policies=ex1$pol, 
#'                           treat_chars=ex1$tx, 
#'                           stagegroups=stages, 
#'                           map=ex1$map, 
#'                           pop_size=20, nsim=2)
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
#' @param policies Data frame specifying policies - see ex1$pol. If 'pairnum' = NA, no treatments are updated. If 'pairnum' is numeric, treatments from that number scenario will be the starting point, with early-detected cases shifted according to the next input
#' @param shifts List of treatment-shift indicators (see shifttreatment_indicator)
#' @param treats List of original treatment matrices (see treatments_by_policy)
#'
#' @examples
#' library(bcimodel)
#' data(ex1) 
#' # ex1$nh shows that there are 4 stage-subgroups. Use a fake random distribution of groups 1:4 for the population before stage-shifting.
#' popdistr <- matrix(sample.int(4, size=40, replace=TRUE), nrow=20, ncol=2)
#' 
#' # Create stageshift indicator matrices for all 3 scenarios: no stage shifts for #1 and #2, but 30% stageshift for #3. Use a small population of size 20, and 2 sims
#' stageshifts <- list(base=matrix(0, nrow=20, ncol=2),
#'                     tam=matrix(0, nrow=20, ncol=2), 
#'                     tamshift=stageshift_indicator(0.85, 20, 2))
#' 
#' # Get the actual stages - only policy #3 has stage-shifting
#' stages <- lapply(stageshifts, shift_stages, original=popdistr, map=ex1$map)
#' 
#' # First 2 scenarios have no shifting, same treatment. Randomly distribute treatment
#' # Scenario 3 has shifting. NAs indicate that treatments will be same as the scenario
#' # with no early detection but same treatments
#' treats <- treatments_by_policy(policies=ex1$pol, 
#'                           treat_chars=ex1$tx, 
#'                           stagegroups=stages, 
#'                           map=ex1$map, 
#'                           pop_size=20, nsim=2)
#' 
#' # Indicators of treatment shifting
#' treatshifts <- lapply(ex1$pol$num,
#'                           shifttreatment_indicator,
#'                           type=ex1$pol$pairnum,
#'                           shifts=stageshifts, 
#'                           basecase=popdistr, 
#'                           map=ex1$map)
#' 
#' newtreat <- update_treat_stageshift(ex1$pol, treatshifts, treats) 
#' # Two cases have been shifted in #3 compared to #2 
#' lapply(treatshifts, table)
#' lapply(newtreat, table)
#' 
#' @return List of updated treatment matrices
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
                                    




