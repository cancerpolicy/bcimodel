
################################################################################
# Function to simulate a series of policies on the same population
# of incident cases
################################################################################


simpolicies <- function(scenarios, naturalhist, treatinfo, 
                        agesource='Standard', minage=0, maxage=100,
                        incsource='Uganda', mortsource='Uganda',
                        popsize=100000, sims=5, futimes=c(5,10), 
                        denom=100000) {

set.seed(98103)
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
survival <- calcmrr(cummortandinc, 
                    length(cumulative_mortality)+1,
                    reverse=TRUE)[1:length(cumulative_mortality)]
yearssaved <- calcarr(cumulative_yearslived, 1, reverse=TRUE)

# Summarized!
table <- compile_outcomes(list(
                    `Cumulative Incidence`=cumulative_incidence,
                    `Cumulative Mortality`=cumulative_mortality,
                    `% Incident Surviving`=survival,
                    `MRR`=mrr, `ARR`=arr, 
                    `Years of Life Saved`=yearssaved),
                          futimes,
                          policynames=ex1$pol$name,
                          pop_size=popsize)



return(table)

} # end simpolicy
