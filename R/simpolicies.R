
################################################################################
# Function to simulate a series of policies on the same population
# of incident cases
################################################################################

#-------------------------------------------------------------------------------
# check_scenarios
#-------------------------------------------------------------------------------
#' Check accurate specification of 'scenarios' parameter to simpolicies()
#' 
#' The 'scenarios' parameter must contain paired policies, i.e. if one scenario has early detection, a scenario must exist that has the same treatments but no early detection. This function checks for the correct pairing
#' 
#' @param scenarios Data frame of scenarios to simulate (see define_scenarios). First scenario should be the base case.
#' @param treatinfo Data frame with treatment hazard ratios, and for each scenario, treatment proportions, by stage-subgroups. See the example and the vignette
#' @return NULL if scenarios are paired properly, or vector listing the problematic pairs
#' @examples
#' data(ex1)
#' # See scenarios
#' ex1$pol
#' # See treatment
#' ex1$tx
#' # Create version that have the "tam" scenario deleted and specify
#' # "base" and "tamandshift" scenarios as paired. This should fail.
#' fail <- ex1
#' fail$pol <- fail$pol[-2,]
#' fail$pol$num[2] <- 2
#' fail$pol$pairnum[2] <- 1
#' fail$tx$tam <- NULL
#' 
#' check_scenarios(ex1$pol, ex1$tx)
#' check_scenarios(fail$pol, fail$tx)

check_scenarios <- function(scenarios, treatinfo) {
    # Check for earlydetection, i.e. any earlydetHR!=1
    warn <- c()
    check_these <- scenarios$id[scenarios$earlydetHR!=1]
    # For each scenario to check, the pairnum
    if (!is.null(check_these)) {
        for (s in check_these) {
            pairnum <- 
                  scenarios$id[scenarios$pairnum[which(scenarios$id==s)]]
            # Check if treatments are the same
            check_tx  <- treatinfo[,s]==treatinfo[,pairnum]
            if (sum(!check_tx)!=0) {
                warn <- c(warn,
                          paste(s, "and", pairnum, 
                                "DO NOT have the same treatment distribution"))
            }
        }
    }
    return(warn)
}

#-------------------------------------------------------------------------------
# simpolicies
#-------------------------------------------------------------------------------
#' Run the model for a series of policies and compare outcomes to the base case
#' Wrapper function for the whole model
#' 
#' @param scenarios Data frame of scenarios to simulate (see define_scenarios). First scenario should be the base case.
#' @param naturalhist Data frame with natural history parameters (see compile_naturalhist)
#' @param treatinfo Data frame with treatment hazard ratios, and for each scenario, treatment proportions, by stage-subgroups. See the example and the vignette
#' @param agesource Country to use for age structure (see data(agestructure) )
#' @param minage Lower age limit for population at sim start
#' @param maxage Upper age limit for population at sim start
#' @param incsource Country to use for incidence rates (see data(incratesf) )
#' @param mortsource Country to use for life table (see data(allmortratesf) )
#' @param popsize Size of population to simulate 
#' @param sims Number of simulations
#' @param futimes Follow-up times at which to tally outcomes, in years
#' @param returnstats Defaults to simulation means, but you can also ask for c("mean", "lower", "upper"). Lower and upper report the 2.5\% and 97.5\% quantiles respectively.
#' @param Denominator/population size by which to report outcomes; default is 100,000
#'
#' @examples
#' # Use example input data loaded in package, object name "ex1"
#' data(ex1) 
#' ex1
#' 
#' # Element "pol" shows the 3 scenarios being modeled
#' "nh" contains the natural history information
#' "map" contains numeric indicators for the stage-subgroups
#' "tx" shows treatment hazard ratios and, for each scenario, the treatment distribution. See vignette for more details
#' 
#' # Run a model using ex1 inputs and the defaults
#' uganda_stdpop <- simpolicies(ex1$pol, ex1$nh, ex1$tx)
#'
#' @return Data frame if returnstats="mean", or list of data frames of length(returnstats)
#' 
#' @export
simpolicies <- function(scenarios, naturalhist, treatinfo, 
                        agesource='Standard', minage=0, maxage=100,
                        incsource='Uganda', mortsource='Uganda',
                        popsize=100000, sims=5, futimes=c(5,10), 
                        returnstats=c('mean'),
                        denom=100000) {

set.seed(98103)
#-------------------------------------------------------------------------------
# Check that scenarios have been specified correctly
#-------------------------------------------------------------------------------
warn <- check_scenarios(scenarios, treatinfo)
if (!is.null(warn)) stop('Error in simpolicies: ', warn)

#-------------------------------------------------------------------------------
# Initialize population
#-------------------------------------------------------------------------------
cat('\nInitializing population...')

pop <- initialize_pop(pop_size=popsize,
                      nsim=sims, 
                      agesource,
                      minage, maxage,
                      incsource, mortsource)


#-------------------------------------------------------------------------------
# Assign base case stage and subgroup
#-------------------------------------------------------------------------------
cat('\nSimulating base case stage and subgroup...')

# Used to be called control_notreat_rows
pop$stagetumor <- add_features(pop$ageclin, probs=naturalhist$prop)

#-------------------------------------------------------------------------------
# Determine stage shifts, if any
#-------------------------------------------------------------------------------
cat('\nDetermining stage shifts...')

# Map of how to shift stages, keeping subgroups the same
stagepairs <- create_stageshift_map(naturalhist)

# Create shift indicators (1=yes, 0=no)
# Generate for all, for simplicity
stageshifts <- lapply(scenarios$earlydetHR, 
                      stageshift_indicator, pop_size=popsize, nsim=sims)

# Shift stage, for advanced cases only
newstages <- lapply(stageshifts, shift_stages, original=pop$stagetumor, 
                    map=stagepairs)
# lapply(newstages, function(x) table(x)/table(pop$stagetumor))

#-------------------------------------------------------------------------------
# Simulate treatments
#-------------------------------------------------------------------------------
cat('\nSimulating treatment received...')

# First create logical indicator of needing to change treatment
shift_treatment <- lapply(scenarios$num,
                          shifttreatment_indicator,
                          type=scenarios$pairnum,
                          shifts=stageshifts, 
                          basecase=pop$stagetumor, 
                          map=stagepairs)

# Sim treatments - THIS IS SLOW
# For early-detection scenarios, we'll sim only early stage treatments
# In the next step, we'll insert those for the stage-shifted cases only
treatments <- treatments_by_policy(policies=scenarios, 
                                   treat_chars=treatinfo, 
                                   stagegroups=newstages, 
                                   map=stagepairs,
                                   popsize,
                                   sims)

# Replace screen_treatments with control_treatments for non-shifted 
# early-stage cases
treatments <- update_treat_stageshift(policies=scenarios,
                                      shifts=shift_treatment,
                                      treats=treatments)

#-------------------------------------------------------------------------------
# Simulate mortality
#-------------------------------------------------------------------------------
cat('\nSimulating cancer mortality...')

# Baseline mortality
mortrates <- lapply(newstages, return_value_from_id, df=naturalhist, 
                    value='mortrate')

# Hazard ratios
HRs <- lapply(treatments, return_value_from_id, df=treatinfo, value='txHR')

# Final rate: hazard ratio times baseline mortality rate
finalmortrates <- lapply(scenarios$num, 
                       FUN=function(x, HR, rate) {
                           HR[[x]]*rate[[x]]
                       },
                       HR=HRs, rate=mortrates)

# Time to cancer death: for scenarios$pairnum=NA, just sim from mortrate.
# For early detection scenarios, keep the same time for the non-shifted 
# cases; for shifted cases, use the quantile from non-early-detection 
# scenario to sim the new time. Two step process, similar to treatment sim
clin2cd <- timetocancerdeath_by_policy(policies=scenarios,
                                       finalmortrates, popsize, sims)

# Replace with new times to cancer death for shifted cases in early detection
# scenario
clin2cd <- update_time_stageshift(policies=scenarios,
                                      shifts=shift_treatment,
                                      rates=finalmortrates,
                                      times=clin2cd)

#-------------------------------------------------------------------------------
# Tabulate time to and cause of death
#-------------------------------------------------------------------------------
cat('\nTabulating time to and cause of death...')

# Compute age at cancer death
ageCD <- lapply(scenarios$num,
            FUN=function(x, ageinc, timetocd) {
                ageinc+timetocd[[x]]
            }, ageinc=pop$ageclin, timetocd=clin2cd)

# Compute time from study start to cancer incidence
timetoInc <- pop$ageclin-pop$ageentry

# Compute time from study start to cancer death
timetoCD <- lapply(scenarios$num,
            FUN=function(x, cancer, entry) {
                cancer[[x]]-entry
            }, cancer=ageCD, entry=pop$ageentry)

# Cause of death
cancerD <- lapply(scenarios$num,
            FUN=function(x, ageCD, ageOC) {
                ifelse(ageCD[[x]]<ageOC,1,0)
            }, ageCD, pop$ageOC)

# Time from trial start to all-cause death
timetoD <- lapply(scenarios$num,
            FUN=function(x, ageCD, ageOC, ageentry) {
                ifelse(ageCD[[x]]<ageOC,
                       ageCD[[x]]-ageentry,
                       ageOC-ageentry)
            }, ageCD, pop$ageOC, pop$ageentry)

#-------------------------------------------------------------------------------
# Summarize outcomes
#-------------------------------------------------------------------------------
cat('\nSummarizing outcomes...')

cumulative_incidence <- cuminc(futimes, timetoInc)
cumulative_mortality <- lapply(scenarios$num,
                               function(x) {
                                   cummort(futimes, timetoD[[x]],
                                           cancerD[[x]])
                               })
cummortandinc <- cumulative_mortality
cummortandinc[[length(cumulative_mortality)+1]] <- 
    cumulative_incidence
cumulative_yearslived <- lapply(scenarios$num,
                               function(x) {
                                   cumyears(futimes, timetoD[[x]])
                               })
mrr <- calcmrr(cumulative_mortality, 1)
arr <- calcarr(cumulative_mortality, 1)
percsurvival <- calcmrr(cummortandinc, 
                    length(cumulative_mortality)+1,
                    reverse=TRUE, perc=TRUE)[1:length(cumulative_mortality)]
yearssaved <- calcarr(cumulative_yearslived, 1, reverse=TRUE)

# Have a scaled version of survival and mrr, so it comes out correctly 
# after being run through compile_outcomes
percsurvivalscaled <- lapply(percsurvival, 
                             function(x) x*popsize/100000)
mrrscaled <- lapply(mrr, function(x) x*popsize/100000)

# Summarized!
table <- compile_outcomes(list(
                            `Cumulative BC Incidence`=cumulative_incidence,
                            `Cumulative BC Mortality`=cumulative_mortality,
                            `% Incident Surviving`=percsurvivalscaled,
                            `MRR`=mrrscaled, `ARR`=arr, 
                            `Years of Life Saved`=yearssaved),
                          futimes,
                          policynames=scenarios$name,
                          pop_size=popsize,
                          stats=returnstats)



return(table)

} # end simpolicy








#-------------------------------------------------------------------------------
# parsimpolicies
#-------------------------------------------------------------------------------
#' Same as simpolicies but calling parallelized sub-functions
#' 
#' @param ncores Number of cores to use
#' @note See help file for simpolicies
#' 
#' 
#' @export

parsimpolicies <- function(scenarios, naturalhist, treatinfo, 
                        agesource='Standard', minage=0, maxage=100,
                        incsource='Uganda', mortsource='Uganda',
                        popsize=100000, sims=5, futimes=c(5,10), 
                        returnstats=c('mean'),
                        denom=100000, ncores=4) {

set.seed(98103)
#-------------------------------------------------------------------------------
# Check that scenarios have been specified correctly
#-------------------------------------------------------------------------------
warn <- check_scenarios(scenarios, treatinfo)
if (!is.null(warn)) stop('Error in parsimpolicies: ', warn)

#-------------------------------------------------------------------------------
# Initialize population
#-------------------------------------------------------------------------------
cat('\nInitializing population...')

pop <- parinitialize_pop(pop_size=popsize,
                      nsim=sims, 
                      agesource,
                      minage, maxage,
                      incsource, mortsource, ncores=ncores)


#-------------------------------------------------------------------------------
# Assign base case stage and subgroup
#-------------------------------------------------------------------------------
cat('\nSimulating base case stage and subgroup...')

# Used to be called control_notreat_rows
pop$stagetumor <- paradd_features(pop$ageclin, probs=naturalhist$prop)

#-------------------------------------------------------------------------------
# Determine stage shifts, if any
#-------------------------------------------------------------------------------
cat('\nDetermining stage shifts...')

# Map of how to shift stages, keeping subgroups the same
stagepairs <- create_stageshift_map(naturalhist)

# Create shift indicators (1=yes, 0=no)
# Generate for all, for simplicity
stageshifts <- lapply(scenarios$earlydetHR, 
                      stageshift_indicator, pop_size=popsize, nsim=sims)

# Shift stage, for advanced cases only
newstages <- lapply(stageshifts, shift_stages, original=pop$stagetumor, 
                    map=stagepairs)
# lapply(newstages, function(x) table(x)/table(pop$stagetumor))

#-------------------------------------------------------------------------------
# Simulate treatments
#-------------------------------------------------------------------------------
cat('\nSimulating treatment received...')

# First create logical indicator of needing to change treatment
shift_treatment <- lapply(scenarios$num,
                          shifttreatment_indicator,
                          type=scenarios$pairnum,
                          shifts=stageshifts, 
                          basecase=pop$stagetumor, 
                          map=stagepairs)

# Sim treatments - THIS IS SLOW
# For early-detection scenarios, we'll sim only early stage treatments
# In the next step, we'll insert those for the stage-shifted cases only
treatments <- partreatments_by_policy(policies=scenarios, 
                                   treat_chars=treatinfo, 
                                   stagegroups=newstages, 
                                   map=stagepairs,
                                   popsize,
                                   sims)

# Replace screen_treatments with control_treatments for non-shifted 
# early-stage cases
treatments <- update_treat_stageshift(policies=scenarios,
                                      shifts=shift_treatment,
                                      treats=treatments)

#-------------------------------------------------------------------------------
# Simulate mortality
#-------------------------------------------------------------------------------
cat('\nSimulating cancer mortality...')

# Baseline mortality
mortrates <- lapply(newstages, return_value_from_id, df=naturalhist, 
                    value='mortrate')

# Hazard ratios
HRs <- lapply(treatments, return_value_from_id, df=treatinfo, value='txHR')

# Final rate: hazard ratio times baseline mortality rate
finalmortrates <- lapply(scenarios$num, 
                       FUN=function(x, HR, rate) {
                           HR[[x]]*rate[[x]]
                       },
                       HR=HRs, rate=mortrates)

# Time to cancer death: for scenarios$pairnum=NA, just sim from mortrate.
# For early detection scenarios, keep the same time for the non-shifted 
# cases; for shifted cases, use the quantile from non-early-detection 
# scenario to sim the new time. Two step process, similar to treatment sim
clin2cd <- timetocancerdeath_by_policy(policies=scenarios,
                                       finalmortrates, popsize, sims)

# Replace with new times to cancer death for shifted cases in early detection
# scenario
clin2cd <- update_time_stageshift(policies=scenarios,
                                      shifts=shift_treatment,
                                      rates=finalmortrates,
                                      times=clin2cd)

#-------------------------------------------------------------------------------
# Tabulate time to and cause of death
#-------------------------------------------------------------------------------
cat('\nTabulating time to and cause of death...')

# Compute age at cancer death
ageCD <- lapply(scenarios$num,
            FUN=function(x, ageinc, timetocd) {
                ageinc+timetocd[[x]]
            }, ageinc=pop$ageclin, timetocd=clin2cd)

# Compute time from study start to cancer incidence
timetoInc <- pop$ageclin-pop$ageentry

# Compute time from study start to cancer death
timetoCD <- lapply(scenarios$num,
            FUN=function(x, cancer, entry) {
                cancer[[x]]-entry
            }, cancer=ageCD, entry=pop$ageentry)

# Cause of death
cancerD <- lapply(scenarios$num,
            FUN=function(x, ageCD, ageOC) {
                ifelse(ageCD[[x]]<ageOC,1,0)
            }, ageCD, pop$ageOC)

# Time from trial start to all-cause death
timetoD <- lapply(scenarios$num,
            FUN=function(x, ageCD, ageOC, ageentry) {
                ifelse(ageCD[[x]]<ageOC,
                       ageCD[[x]]-ageentry,
                       ageOC-ageentry)
            }, ageCD, pop$ageOC, pop$ageentry)

#-------------------------------------------------------------------------------
# Summarize outcomes
#-------------------------------------------------------------------------------
cat('\nSummarizing outcomes...')

cumulative_incidence <- cuminc(futimes, timetoInc)
cumulative_mortality <- lapply(scenarios$num,
                               function(x) {
                                   cummort(futimes, timetoD[[x]],
                                           cancerD[[x]])
                               })
cummortandinc <- cumulative_mortality
cummortandinc[[length(cumulative_mortality)+1]] <- 
    cumulative_incidence
cumulative_yearslived <- lapply(scenarios$num,
                               function(x) {
                                   cumyears(futimes, timetoD[[x]])
                               })
mrr <- calcmrr(cumulative_mortality, 1)
arr <- calcarr(cumulative_mortality, 1)
percsurvival <- calcmrr(cummortandinc, 
                    length(cumulative_mortality)+1,
                    reverse=TRUE, perc=TRUE)[1:length(cumulative_mortality)]
yearssaved <- calcarr(cumulative_yearslived, 1, reverse=TRUE)

# Have a scaled version of survival and mrr, so it comes out correctly 
# after being run through compile_outcomes
percsurvivalscaled <- lapply(percsurvival, 
                             function(x) x*popsize/100000)
mrrscaled <- lapply(mrr, function(x) x*popsize/100000)

# Summarized!
table <- compile_outcomes(list(
                            `Cumulative BC Incidence`=cumulative_incidence,
                            `Cumulative BC Mortality`=cumulative_mortality,
                            `% Incident Surviving`=percsurvivalscaled,
                            `MRR`=mrrscaled, `ARR`=arr, 
                            `Years of Life Saved`=yearssaved),
                          futimes,
                          policynames=scenarios$name,
                          pop_size=popsize,
                          stats=returnstats)



return(table)

} # end simpolicy


