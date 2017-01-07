############################################################
# System time tests
############################################################


#----------------------------------------------------------
# One long vector versus a matrix
#----------------------------------------------------------
popsize <- 100000
nsim <- 100
probs <- c(0.1,0.2,0.6,0.1)
names <- c('a', 'b', 'c', 'd')

# Helper function to do one rep
one_rep <- function(id, popsize, probs, names) {
    # Do the draw from rmultinom
    draws <- rmultinom(popsize, size=1, probs)
    # Match to names
    indices <- apply(draws, 2, function(x) which(x==1))
    return(t(names[indices]))
}

onelong_sim_multinom <- function(popsize,nsim,probs,names) {
    return(one_rep(1,popsize*nsim, probs, names))
}

systest_sim_multinom <- function(popsize,nsim,probs, names) {
    # Now do the replicates
    all_reps <- sapply(1:nsim, FUN=one_rep, popsize, probs, names)
    return(all_reps)
}

parallel_sim_multinom <- function(popsize,nsim,probs,names,ncores) {
    cl <- makeCluster(getOption("cl.cores", ncores)) 
    v <- parSapply(cl,1:nsim, FUN = one_rep, popsize, probs, names)
    stopCluster(cl)
    return(v)
}

# Should move this to a vignette
if (1==0) {
    # One long vector: elapsed is 44.339
    system.time(onelong_sim_multinom(popsize,nsim,probs,names))

    # Matrix: elapsed is 32.505
    system.time(systest_sim_multinom(popsize,nsim,probs,names))

    # Matrix with rows and cols reversed: elapsed is 30.04
    system.time(systest_sim_multinom(nsim,popsize,probs,names))

    # Parallel processing with 10 cores: elapsed is 19.23
    library(parallel)
    system.time(parallel_sim_multinom(popsize,nsim,probs,names,10))
    # Parallel processing with 4 cores: elapsed is 15.360
    system.time(parallel_sim_multinom(popsize,nsim,probs,names,4))
    # Parallel processing with 1 cores: elapsed is 29.27
    system.time(parallel_sim_multinom(popsize,nsim,probs,names,1))
}


