
#-------------------------------------------------------------------------------
# cumsurv_to_exprate
#-------------------------------------------------------------------------------
#' Convert survivals to rates
#' @param cum.surv  Proportion surviving
#' @param year Time corresponding to cum.surv, i.e. time since time=0, cum.surv=1
#' @param haz Hazard ratio, of the cum.surv reflects an original curve raised to a hazard
#' 
#' @export
cumsurv_to_exprate = function(cum.surv, year=10, haz=1) {log(cum.surv)/(-year*haz)}

#-------------------------------------------------------------------------------
# sim_multinom
#-------------------------------------------------------------------------------
#' Simulate a category assignment using a multinomial distribution
#' 
#' @param nsims Number of size=1 sims from the multinomial
#' @param nreps Number of times to replicate the nsims. Put the smaller number here, the bigger one as nsims
#' @param probs Vector of probabilities for each of the K categories
#' @param names Vector of names associated with each of the categories. Must be of same length as probs. Can be character or numeric. Use numeric row IDs if you're referring to a group defined by multiple variables whose groupings are defined in a separate data frame
#' @examples
#' sim_multinom(10, 5, c(0.1, 0.3, 0.6), names=c('a', 'b', 'c'))
#' sim_multinom(5, 10, c(0.1, 0.3, 0.6), names=c('a', 'b', 'c'))
#'
#' @export

sim_multinom <- function(nsims, nreps, probs, names) {

    # Helper function to do one rep
    one_rep <- function(id, nsims, probs, names) {

        # Do the draw from rmultinom
        draws <- rmultinom(nsims, size=1, probs)

        # Match to names
        indices <- apply(draws, 2, function(x) which(x==1))
        return(t(names[indices]))
    }

    # Now do the replicates
    all_reps <- sapply(1:nreps, FUN=one_rep, nsims, probs, names)

    return(all_reps)
}
