
#############################################################################
# Functions from the cantrance package that need to be improved ultimately
#############################################################################



############################################################
## build_simpop
############################################################
#' Build a simulated population by bootstrapping from a list
#'
#' Cantrance remnant; use row IDs to reconstruct a single bootstrap of data
#'
#' @param datalist List of data frames containing the data to be bootstrapped then combined into a single data frame. Each data frame should refer to a different characterisic of individuals in the population, e.g. 1st data frame is sex, 2nd data frame is age
#' @param rowlist List of data frames containing row numbers for bootstrapping the data in datalist. Row numbers of rowlist[[i]] refer to datalist[[i]]
#' @param sim Simulation number, i.e. which column of rowlist elements to use
#'
#' @return Data frame with characteristics from datalist in proportion to their appearance in rowlist[[sim]]
#'
#' @examples
#' # Possible characteristics
#' pop <- list(data.frame(male=0),
#'             data.frame(age=c(20,30,40)))
#' # Actual characteristics for each of 2 sims, specified by rows, for pop size = 10
#' rows <- list(matrix(1, nrow=10, ncol=2),
#'              matrix(sample.int(3, 20, replace=TRUE), nrow=10, ncol=2))
#' # Population using sim 1
#' build_simpop(pop, rows, 1)
#' build_simpop(pop, rows, 2)
#' 
#' @export

build_simpop = function(datalist, rowlist, sim) {
    do.call('cbind',
            lapply(1:length(datalist),
                   function(i) {
                       datalist[[i]][rowlist[[i]][, sim], , drop=FALSE]
                   }))
}


############################################################
## sim_KM
############################################################
#'Returns random deviates from a survival curve
#'
#'Takes in event times and corresponding survival probs (as from a KM curve) and returns n random deviates
#'
#'@param survival Vector of survival probabilities, as from a KM curve
#'@param time Vector of event times corresponding to the survival probabilities
#'@param smalltimes Value to be returned if survival is higher/time is smaller than any observed. Suggestions: 0, min(times), NA
#'@param bigtimes Value to be returned if survival is lower/time is higher than any observed. Suggestions: Inf, max(times), NA
#'@param nsims Number of deviates to return
#'@param mindraw Minimum draw allowed on the survival scale
#'@param maxdraw Maximum draw allowed on the survival scale
#'@param draws If draws on the survival scale have already been made, specify as a vector here
#'
#'@return Times to event representative of the input curve
#'@export

sim_KM = function(survival, time, smalltimes, bigtimes, nsims, mindraw=0, 
                  maxdraw=1, draws=NULL) {

    # Take n_sim random draws between mindraw and maxdraw
    if (is.null(draws))
        draws = runif(nsims, min=mindraw, max=maxdraw)
    
    # Now get the corresponding times 
    ##note<< Uses a linear approximation to interpolate between observed events
    ttrs = approx(x=survival,
                  y=time,
                  method="linear",
                  xout=draws,
                  yleft=bigtimes,
                  yright=smalltimes)[["y"]]
    return(ttrs)
}


############################################################
## sim_same_qexp
############################################################

sim_same_qexp = function# Simulate new time to event using quantile from old time to event

(oldtime,
    ### Vector of time to event using old estimation method
 oldrate,
    ### Vector of rates for exponential distribution used to
    ### estimate oldtime
 newrate,
    ### Vector of rates for exponential distribution to be
    ### used to estimate new times to event
 prefix
    ### Prefix for column names in results vector/matrix
) {

    # Estimate new time to event using exponential
    # distribution quantile from old time to event
    newtime = qexp(p=pexp(oldtime, rate=oldrate),
                    rate=newrate)
    
    # Name columns with prefix
    colnames(newtime) = paste0(prefix, 1:ncol(newtime))
    
    # Correct nuances due to rounding
    newtime[abs(oldtime-newtime)<0.001] = 
        oldtime[abs(oldtime-newtime)<0.001]
    
    # Return
    return(newtime)
    
### Vector of time to event using same quantile as old time
### to event
}


############################################################
## calc_ac_lifespan
############################################################

#' Generate an age at other-cause death from a lifetable
#' 
#' Given a year of birth and current age, draws a random age at other-cause death from a lifetable
#' 
#' @param ageentry An integer age
#' @param male Sex indicator: 1 for male, 0 for female
#' @param lifetable Data frame with columns Survival, BirthCohort, Age, and Male
#' @param n_sim Number of random draws to return
#' @param time Set to TRUE to return time between ageentry and death, instead of age at death
#' @param haz Hazard to apply to the life table survival
#' @param max100 If TRUE, interpolate out to 100 as the max survival age
#' 
#' @return A vector of length n_sim containing ages at other-cause death (or times to other-cause death, if time=TRUE)
#' 
#' @examples
# Use a Ugandan life table
#' library(bcimodel)
#' data(allmortratesf)
#' surv <- interpolate_cumsurv(allmortratesf, 
#'                           ratevar='Rate.Per.100K',
#'                           country='Uganda')
#' # Code compatibility tweaks for life table
#' surv <- transform(surv, Age=age, Survival=cumsurv, Male=0)
#' # Generate 5 death times for a female of age 30
#' calc <- ac <- lifespan(30, 0, surv, 5)
#'
#' @export

calc_ac_lifespan = function(ageentry, male, lifetable, n_sim=100, time=FALSE, 
                            haz=1, max100=TRUE) {

    # Check that sort is as expected
    if (lifetable$Survival[nrow(lifetable)]>
        lifetable$Survival[1]) stop('Sort by descending survival')

    if (max100) {
        # Interpolate survival out to age 100
        itime <- (max(lifetable$Age)+1):100
        isurv <- approx(c(max(lifetable$Age), 100),
                        c(min(lifetable$Survival), 0),
                        xout=itime)
        times <- c(lifetable$Age, isurv$x)
        survival <- c(lifetable$Survival, isurv$y)
        maxage=100
    } else {
        times <- lifetable$Age
        survival <- lifetable$Survival
        maxage=max(lifetable$Age)
    }
    
    # Apply hazard modifier
    if (haz!=1) survival <- survival^haz

    # Get the max cumulative survival for specified age 
    if (!round(ageentry)==ageentry) stop('Age at entry must be integer')
    maxu = survival[which(times==ageentry)]
    
    # Take random draws between 0 and the survival
    # probability. Because this is cumulative survival, 0
    # survival corresponds to the upper bound for age at
    # death. Next,linearly interpolate between age and
    # survival estimates and determine the age that
    # corresponds to the draw. Use maxage as the maximum age
    # possible (survival=0)
    death <- sim_KM(survival, times, 
                    smalltimes=ageentry, 
                    bigtimes=maxage, 
                    nsims=n_sim, 
                    mindraw=0, 
                    maxdraw=maxu)
    
    # Convert the age into a time from entry age until
    # death?
    if (time) death = death-rep(ageentry, n_sim)

    # Check for errors
    if (sum(is.na(death))!=0 | sum(death==0)!=0) 
        browser("In calc_ac_lifespan(), death values equal NA or O")
    return(death)

}


############################################################
# calc_ac_lifespan_pop
############################################################

#' Use the individual calc_ac_lifespan function to estimate lifespans of a population
#' 
#' Given a population with birth years and ages, returns random draws of their ages at other-cause death
#'
#' @param popdata Data frame where individuals are rows, with columns birth_year, age, and male, OR a list of data frames that allow build_simpop() to construct this
#' @param bootrows Matrix or data frame of row indicators that can be applied to the data to recover different bootstraps of the data. Each column is a different bootstrap of thedata. OR, list of data frames with these row indicators for use with build_simpop()
#' @param life_table Data frame with columns Survival, BirthCohort, Age, and Male
#' @param results_as_matrix Set to TRUE to convert results from matrix to data frame
#' @param survHR Hazard to apply to the life table survival
#' @param max100_topass If TRUE, interpolate out to 100 as the max survival age
#' 
#' @return A data frame or matrix with nsim columns of randomly drawn ages at other-cause death for individuals (rows)
#' @examples 
# Use a Ugandan life table
#' library(bcimodel)
#' data(allmortratesf)
#' surv <- interpolate_cumsurv(allmortratesf, 
#'                           ratevar='Rate.Per.100K',
#'                           country='Uganda')
#' # Code compatibility tweaks for life table
#' surv <- transform(surv, Age=age, Survival=cumsurv, Male=0)
#' 
#' # Example 1: use the list approach
#' # Simple age distribution
#' ages <- data.frame(age=c(20,30,40), prop=c(0.2, 0.5, 0.3))
#' # Simulate a population of size 10 twice, using these proportions and the row IDs
#' sim.age.rows <- sim_multinom(nsims=10, 2, ages$prop, names=1:nrow(ages))
#' # Now specify they're all female
#' sex <- data.frame(male=0, prop=1)
#' sim.sex.rows <- matrix(1, nrow=10, ncol=2)
#' # Create the lists
#' poplist <- list(sex, ages)
#' rowlist <- list(sim.sex.rows, sim.age.rows)
#' # Get ages at OC death
#' oc <- calc_ac_lifespan_pop(poplist, rowlist, surv)
#' 
#' # Example 2
#' # Alternatively, pass a population in data frame form and just use its normal rows
#' # Use build_simpop and the 2nd simulation's row IDs
#' pop <- build_simpop(poplist, rowlist, sim=2)[,c('male', 'age')]
#' poprows <- replicate(2, 1:nrow(pop))
#' oc <- calc_ac_lifespan_pop(pop, poprows, surv)
#' 
#' @export

calc_ac_lifespan_pop = function(popdata, bootrows, life_table, results_as_matrix=FALSE, 
                                survHR=1, max100_topass=TRUE) {
    # Determine number of simulations
    nsim = ncol(bootrows)
    if (is.null(nsim)) nsim = ncol(bootrows[[1]])

    # Estimate for each simulation separately
    newpop = sapply(1:nsim,
        function(y) {
            # Construct the dataset
            if (!is.data.frame(popdata)) {
                # This was the cantrance way of getting one bootstrap of the data
                thisboot = subset(build_simpop(datalist=popdata,
                                               rowlist=bootrows,
                                               sim=y),
                                  select=c(age, male))
            } else {
                thisboot = popdata[bootrows[, y], 
                                   c("age",
                                     "male")]
            }
            thisboot$tempid = 1:nrow(thisboot)
            
            # Simulate ages at other-cause death
            ocs = ddply(thisboot,
                        .(age, male),
                        function(x) {
                            # Extract age, birth year, and 
                            # male specifications
                            age = min(x$age)
                            male = min(x$male)

                            # Calculate all-cause lifespan
                            if (is.data.frame(x)) nrow=nrow(x) else nrow=1
                            ageOCs = calc_ac_lifespan(ageentry=age, 
                                                      male=male, 
                                                      n_sim=nrow, 
                                                      haz=survHR, 
                                                      lifetable=life_table,
                                                      max100=max100_topass) 
                            ageOCs = data.frame(x, ageOC=matrix(ageOCs, 
                                                                nrow=nrow, ncol=1))
                            return(ageOCs)
                        })
            
            # Return to original sort order, and return 
            # ageOC
            ocs = ocs[order(ocs$tempid), ]
            return(ocs$ageOC)
        })
    
    # Return results as a matrix?
    if (results_as_matrix) {
        colnames(newpop) = 1:nsim
        newpop = as.matrix(newpop)
    }

    return(newpop)

}
