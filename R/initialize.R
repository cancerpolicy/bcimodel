############################################################
# Functions to initialize a population with specified 
# age structure, cancer incidence, all-cause mortality,
# and stage-receptor status
############################################################

#-------------------------------------------------------------------------------
# initialize_pop
#-------------------------------------------------------------------------------
#' Initialize a FEMALE population with an age structure, dates of cancer 
#' incidence, and dates of all-cause mortality
#' 
#' @param agesource Country to use for age structure (see data(agestructure) )
#' @param minage Lower age limit for population at sim start
#' @param maxage Upper age limit for population at sim start
#' @param incsource Country to use for incidence rates (see data(incratesf) )
#' @param mortsource Country to use for life table (see data(allmortratesf) )
#' @param pop_chars A list of data frames that specify additional features
#' to simulate in the population. Defaults to giving the whole population
#' male=0, i.e. all female sex. 
#' @examples
#'pop <- initialize_pop(pop_size=100000,
#'                      nsim=2, 
#'                      agesource='Standard', 
#'                      minage=0, maxage=100, 
#'                      incsource='Uganda', 
#'                      mortsource='Uganda')
#'
#' @export

initialize_pop <- function(pop_size, nsim,
                           agesource, minage, maxage, 
                           incsource, mortsource,
                           pop_chars=list(male=data.frame(male=c(0), 
                                                          prop=c(1)))) 
{
    # Load databases
    data(incratesf)
    data(allmortratesf)
    data(agestructure)

    # Compute survival from incidence/mortality databases
    inc <- interpolate_cumsurv(incratesf, 
                              ratevar='Female.Rate.Per.100K',
                              country=incsource,
                              maxage=100)
    mort <- interpolate_cumsurv(allmortratesf, 
                              ratevar='Rate.Per.100K',
                              country=mortsource)
    # Code compatibility tweaks
    mort <- transform(mort, Age=age, Survival=cumsurv, Male=0)

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
                            sim_multinom(nsims=Npop, 
                                         nreps=Nsim,
                                         probs=x$prop,
                                         names=1:nrow(x))
                            }, pop_size, nsim)

    # Ages at entry 
    ageentry <- return_value_from_id(ids=pop_chars_rows[['age']],
                                     df=pop_chars[['age']],
                                     value='age')
    
    # Ages at other-cause death (load alternative life table if desired)
    ageOC <- calc_ac_lifespan_pop(popdata=pop_chars,
                                    bootrows=pop_chars_rows,
                                    life_table=mort,
                                    results_as_matrix=TRUE)

    # Ages at cancer incidence
    ageclin <- sim_clinical_incidence(popdata=pop_chars,
                               bootrows=pop_chars_rows,
                               incidence=inc,
                               results_as_matrix=TRUE)

    return(list(ageentry=ageentry, ageOC=ageOC, ageclin=ageclin))
}

#-------------------------------------------------------------------------------
# compile_naturalhist
#-------------------------------------------------------------------------------
#' Creates a 'naturalhist' object that describes key natural history parameters
#' for incident cases
#' 
#' Takes in stats on cancer stage, subgroup (e.g. tumor type) and survival
#' and compiles into a data frame of class 'naturalhist'. Note that this function
#' could be modified to allow for more complex situations. Right now it presumes
#' that mortality rates vary by stage only, and that stage and subgroup are
#' uncorrelated.
#'
#' @param prop_adv Advanced Proportion of cancers presenting as Advanced stage
#' @param mortrates Named vector of cancer mortality rates by stage,
#' e.g. c(Early=.05, Advanced=0.21)
#' @param subgroup_probs Named vector of subgroup probabilities that occur within stage groups, e.g. c(`ER+`=0.5, `ER-`=0.5)
#' @return Data frame of class 'naturalhist' with columns stage, subgroup,
#' mortrate and prop
#' 
#' @examples
#' compile_naturalhist(prop_adv=0.85, mortrates=c(Early=0.05, Advanced=0.21), 
#'                    subgroup_probs=c(`ER+`=0.5, `ER-`=0.5))
#'
#' @export

compile_naturalhist <- function(prop_adv, mortrates, subgroup_probs) {
    
    stage_probs <- c(Early=1-prop_adv, Advanced=prop_adv)
    df <- lapply(names(mortrates),
                       function(x) {
                           return(data.frame(prop=subgroup_probs*stage_probs[x],
                                              stage=x,
                                              subgroup=names(subgroup_probs),
                                              mortrate=mortrates[x],
                                              stringsAsFactors=FALSE))
                        })
    
    df <- ldply(df)
    if (round(sum(df$prop),1)!=1) stop('Check that subgroup_probs sum to 1')
    class(df) <- append(class(df), 'naturalhist')
    return(df)
}

#-------------------------------------------------------------------------------
# add_features
#-------------------------------------------------------------------------------
#' Add features like stage and tumor type to an initialized population
#' 
#' @param popbysim A matrix of dimensions (population size)x(number of sims)
#' @param probs A vector of probabilities summing to 1. Each probability
#' represents the probability of being in each feature group
#' stage and tumor-type subgroup
#' @param stagetumor Optional vector of names corresponding to the 
#' stage-tumor subgroups in probs. Defaults to 1:length(probs)
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
#'
#' @export

add_features <- function(popbysim, probs, names=NULL) {

    if (is.null(names)) names <- 1:length(probs)
    if (!is.null(names)) {
        if (length(names)!=length(probs)) stop('No 1:1 match between names 
                                              and probs')
    }
    if (round(sum(probs),1)!=1) stop('Probabilities do not sum to 1')
    
    return(sim_multinom(nsims=nrow(popbysim), 
                        nreps=ncol(popbysim), probs, names))
}

