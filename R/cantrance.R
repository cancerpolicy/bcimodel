
#############################################################################
# Functions from the cantrance package that need to be improved ultimately
#############################################################################



############################################################
## build_simpop
############################################################

build_simpop = function# Build a simulated population by bootstrapping from a list

(datalist, 
    ### List of data frames containing the data to be 
    ### bootstrapped then combined into a single data frame
 rowlist,
    ### List of data frames containing row numbers for 
    ### bootstrapping the data in datalist
 sim
    ### Simulation number. Indicates which column of rowlist
    ### to use
 ) {
    do.call('cbind',
            lapply(1:length(datalist),
                   function(i) {
                       datalist[[i]][rowlist[[i]][, sim], , drop=FALSE]
                   }))
}


############################################################
## sim_KM
############################################################

sim_KM = function# Returns random deviates from a survival curve

##description<< Takes in event times and corresponding
## survival probs (as from a KM curve) and returns n random
## deviates

(survival, 
    ### Vector of survival probabilities, as from a KM curve
time, 
    ### Vector of event times corresponding to the survival
    ### probabilities
smalltimes, 
    ### Value to be returned if survival is higher/time is
    ### smaller than any observed. Suggestions: o,
    ### min(times), NA
bigtimes, 
    ### Value to be returned if survival is lower/time is
    ### higher than any observed. Suggestions: Inf,
    ### max(times), NA
nsims,
    ### Number of deviates to return
mindraw=0, 
    ### Minimum draw allowed on the survival scale
maxdraw=1, 
    ### Maximum draw allowed on the survival scale
draws=NULL
    ### If draws on the survival scale have already been
    ### made, specify as a vector here
) {

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
### Times to event representative of the input curve
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

calc_ac_lifespan = function# Generates an age at other-cause death from US lifetables

##description<< Given a year of birth and current age, draws 
## a random age at other-cause death from US cohort
## lifetables. 

(ageentry, 
    ### Current age
male, 
    ### Sex: male=1 for male, male=0 for female
lifetable,
    ### Life table dataframe with columns Survival,
    ### BirthCohort, Age, and Male
n_sim=100, 
    ### Number of random draws to return
time=FALSE,
    ### Return age at other-cause death, or time from
    ### current age to other-cause death?
haz=1,
    ### Should the US lifetables be modified by a hazard
    ### ratio before drawing from them? If so, the cohort
    ### survival will be raised to this number. 
max100=TRUE
    ### Interpolate out to 100 as the max
) {

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

### A vector of length n_sim containing ages at other-cause
### death (or times to other-cause death, if time=TRUE)
}


############################################################
# calc_ac_lifespan_pop
############################################################

calc_ac_lifespan_pop = function# Use the individual calc_ac_lifespan function to estimate lifespans of a population

##description<< Given a population with birth years and
## ages, returns random draws of their ages at other-cause
## death

(popdata, 
    ### Data frame where individuals are rows, with columns 
    ### birth_year, age, and male
bootrows,
    ### Matrix/df of row indicators that can be applied to
    ### the data to recover different bootstraps of the
    ### data. Each column is a different bootstrap of the
    ### data
life_table,
results_as_matrix=FALSE,
    ### Convert the results from a data frame to a matrix?
survHR=1,
    ### Hazard ratio for survival as a modification of the
    ### life table in use
max100_topass=TRUE
    ### Interpolate out to max age of 100? TRUE/FALSE pass on to 
    ### calc_ac_lifespan
) {
    # Determine number of simulations
    nsim = ncol(bootrows)
    if (is.null(nsim)) nsim = ncol(bootrows[[1]])

    # Estimate for each simulation separately
    newpop = sapply(1:nsim,
        function(y) {
            # Construct the dataset
            if (!is.data.frame(popdata)) {
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

### A data frame or matrix with nsim columns of randomly
### drawn ages at other-cause death for individuals (rows)
}
