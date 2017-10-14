
#-------------------------------------------------------------------------------
# partreatments_by_policy
#-------------------------------------------------------------------------------
#' Use sim_treatment_by_subgroup to simulate treatments, parallelized
#' 
#' Simulate treatments according to specified policy rules
#' 
#' @param policies A "scenarios" data frame containing an 'id' for the policies
#' and a 'pairnum' column indicating either NA or the paired policy, for
#  strategies with early detection
#' @param method Either 'outer' or 'inner'. If outer, will parallelize the outer
#' lapply. If inner, will parallelize repeated calls to sim_multinom

#' @export

partreatments_by_policy <- function(policies, treat_chars, stagegroups,
                                 map, pop_size, nsim, ncores=4, method='inner') {

    if (method=='inner') {
        cl <- makeCluster(getOption("cl.cores", ncores)) 
        clusterSetRNGStream(cl, iseed = 98103)
        clusterExport(cl=cl, varlist=c("treat_chars", "stagegroups", "map",
                                       "policies", "pop_size", "nsim"),
                      envir=environment())

        # Should really make this function separate
        treats <- parLapply(cl, policies$id, 
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
        stopCluster(cl)

    } else if (method=='inner') {
        treats <- parLapply(cl, policies$id, 
                         function(x) {
                             policy <- policies$num[which(policies$id==x)]
                             pairnum <- policies$pairnum[which(policies$id==x)]
                             if (is.na(policies$pairnum[which(policies$id==x)])) {
                                 t <- parsim_treatment_by_subgroup(treat_chars, 
                                                                stagegroups[[1]], 
                                                                x, pop_size, nsim,
                                                                ncores)
                             } else {
                                 # This sims treatments only for the early
                                 # stages - leaves NA's for the advanced stages
                                 # Use the stagegroup matrix that reflects the
                                 # stages AFTER stage-shifting
                                 earlytreat <- subset(treat_chars, 
                                                      SSno%in%map['Early',])
                                 t <- parsim_treatment_by_subgroup(earlytreat, 
                                                                stagegroups[[policy]], 
                                                                x, pop_size, nsim,
                                                                ncores)
                             }
                             return(t)
                         })
    }
    return(treats)

}


#-------------------------------------------------------------------------------
# parsim_multinom
#-------------------------------------------------------------------------------

#' Simulate from mutinomial, parallelized
#'
#' Simulate matrix of draws from multinomial
#' 
#' @param nsims Number of size=1 sims from the multinomials
#' @param nreps Number of times to replicate the nsims. Put the smaller number here, the bigger one as nsims
#' @param probs Vector of probabilities for each of the K categories
#' @param names Vector of names associated with each of the categories. Must be of same length as probs. Can be character or numeric. Use numeric row IDs if you're referring to a group defined by multiple variables whose groupings are defined in a separate data frame

#' @export 

parsim_multinom <- function(nsims, nreps, probs, names, ncores=4) {

    # Helper function to do one rep
    one_rep <- function(id, nsims, probs, names) {

        # Do the draw from rmultinom
        draws <- rmultinom(nsims, size=1, probs)

        # Match to names
        indices <- apply(draws, 2, function(x) which(x==1))
        return(t(names[indices]))
    }

    cl <- makeCluster(getOption("cl.cores", ncores)) 
    clusterSetRNGStream(cl, iseed = 98103)
    clusterExport(cl=cl, 
                  varlist=c("nreps", "one_rep", "nsims", "probs", "names"), 
                  envir=environment())

    # Now do the replicates
    all_reps <- parSapply(cl, 1:nreps, FUN=one_rep, nsims, probs, names)

    stopCluster(cl)
    return(all_reps)
}


#-------------------------------------------------------------------------------
# parinitialize_pop
#-------------------------------------------------------------------------------
#' Initialize a FEMALE population with an age structure, dates of cancer 
#' incidence, and dates of all-cause mortality - parallelized version
#' 
#' @param agesource Country to use for age structure (see data(agestructure) )
#' @param minage Lower age limit for population at sim start
#' @param maxage Upper age limit for population at sim start
#' @param incsource Country to use for incidence rates (see data(incratesf) )
#' @param mortsource Country to use for life table (see data(allmortratesf) )
#' @param pop_chars A list of data frames that specify additional features
#' to simulate in the population. Defaults to giving the whole population
#' male=0, i.e. all female sex. 
#' param ncores Number of cores
#' @examples
#'pop <- initialize_pop(pop_size=100000,
#'                      nsim=2, 
#'                      agesource='Standard', 
#'                      minage=0, maxage=100, 
#'                      incsource='Uganda', 
#'                      mortsource='Uganda')

#' @export

parinitialize_pop <- function(pop_size, nsim,
                           agesource, minage, maxage, 
                           incsource, mortsource,
                           pop_chars=list(male=data.frame(male=c(0), 
                                                          prop=c(1))),
                           ncores=4) 
{
    # Load databases
    data(incratesf)
    data(allmortratesf)
    data(agestructure)

    # Compute survival from incidence/mortality databases
    # Edit 8/21/17: removed "maxage=100" argument from interpolate_cumsurv
    # Edit 10/10/17: returned the maxage arg to interpolate_cumsurv for incidence,
    # since otherwise it breaks down for maxage > 87 (the incidence data limit). 
    # Use only if maxage > 87.
    if (maxage>87) maxIncAge=100 else maxIncAge=NULL
    inc <- interpolate_cumsurv(incratesf, 
                              ratevar='Female.Rate.Per.100K',
                              country=incsource,
                              maxage=maxIncAge)
    # Edit 8/21/17: allow use of BMD cohort life tables
    if (!grepl('Birth Cohort 1950', mortsource)) {
        mort <- interpolate_cumsurv(allmortratesf, 
                                  ratevar='Rate.Per.100K',
                                  country=mortsource)
        # Code compatibility tweaks
        mort <- transform(mort, Age=age, Survival=cumsurv, Male=0)
        interpolate_to_100=TRUE
    } else {
        data(cohortltf)
        mort <- subset(cohortltf, BirthCohort==1950)
        interpolate_to_100=FALSE
    }
        

    # Add age to pop_chars using parameter choices
    pop_chars[['age']] <- format_age(subset(agestructure, Country==agesource,
                                            select=c('age', 'pop', 'prop')), 
                                     minAge=minage, maxAge=maxage)

    # First, simulate the indepenent characteristics given in
    # pop_chars. Right now it assumes that each element in the list
    # refers to a single variable rather than a joint distribution
    # of variables. Very easy to generalize to return the row #
    # of the original dataset rather than the value of a single
    # variable. It would be easy to instead use Leslie's 
    # create_pop_list() function, and that would work with
    # a more complex age pattern, too
    pop_chars_rows <- lapply(pop_chars, function(x, Npop, Nsim) {
                            parsim_multinom(nsims=Npop, 
                                         nreps=Nsim,
                                         probs=x$prop,
                                         names=1:nrow(x),
                                         ncores)
                            }, pop_size, nsim)

    # Ages at entry 
    ageentry <- return_value_from_id(ids=pop_chars_rows[['age']],
                                     df=pop_chars[['age']],
                                     value='age')
    
    # Ages at other-cause death (load alternative life table if desired)
    ageOC <- calc_ac_lifespan_pop(popdata=pop_chars,
                                    bootrows=pop_chars_rows,
                                    life_table=mort,
                                    results_as_matrix=TRUE,
                                    max100_topass=interpolate_to_100)

    # Ages at cancer incidence
    ageclin <- sim_clinical_incidence(popdata=pop_chars,
                               bootrows=pop_chars_rows,
                               incidence=inc,
                               results_as_matrix=TRUE)

    return(list(ageentry=ageentry, ageOC=ageOC, ageclin=ageclin))
}

#-------------------------------------------------------------------------------
# paradd_features
#-------------------------------------------------------------------------------
#' Add features like stage and tumor type to an initialized population - parallelized
#' 
#' @param popbysim A matrix of dimensions (population size)x(number of sims)
#' @param probs A vector of probabilities summing to 1. Each probability
#' represents the probability of being in each feature group
#' stage and tumor-type subgroup
#' @param stagetumor Optional vector of names corresponding to the 
#' stage-tumor subgroups in probs. Defaults to 1:length(probs)
#' param ncores Number of cores
#' @return Matrix of dimensions (pop size)x(nsim) indicating stage-tumor 
#' status according to the 'names' parameter
#' Defaults to 1:length(probs)
#'
#' @examples
#' Say there are four stage-ER receptor subgroups that have equal probability,
#' in a population of size 20 and 2 simulations
#' add_features(matrix(1, nrow=20, ncol=2), probs=c(0.25, 0.25, 0.25, 0.25))
#' Could ID the groups by name
#' add_features(matrix(1, nrow=20, ncol=2), probs=c(0.25, 0.25, 0.25, 0.25),
#'              names=c('EarlyER+', 'EarlyER-', 'AdvancedER+', 'AdvancedER-'))

#' @export

paradd_features <- function(popbysim, probs, names=NULL, ncores=4) {

    if (is.null(names)) names <- 1:length(probs)
    if (!is.null(names)) {
        if (length(names)!=length(probs)) stop('No 1:1 match between names 
                                              and probs')
    }
    if (round(sum(probs),1)!=1) stop('Probabilities do not sum to 1')
    
    return(parsim_multinom(nsims=nrow(popbysim), 
                        nreps=ncol(popbysim), probs, names, ncores))
}

#-------------------------------------------------------------------------------
# parsim_treatment_by_subgroup
#-------------------------------------------------------------------------------
#' Simulate treatment by subgroups - parallelized

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
#' @param ncores Number of cores
#' 
#' @export

parsim_treatment_by_subgroup = function(treat_chars, stage_subgroup_rows, thisprop, 
                                     pop_size, nsim, ncores=4) {
    
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
            parsim_multinom(nsims=pop_size, 
                         nreps=nsim,
                         probs=tempid[,thisprop],
                         names=tempid[,'txSSno'],
                         ncores)
        # Assign values to those in group g
        results[stage_subgroup_rows==g] <-
            treat_choice[stage_subgroup_rows==g]
    }

    return(results)
}
