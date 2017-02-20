
############################################################
# Functions to model early detection
############################################################

#-------------------------------------------------------------------------------
# create_stageshift_map
#-------------------------------------------------------------------------------
#' Creates a "map" of how the stage shift will be applied
#'
#' Follows the principle that when an advanced cases is shifted to early stage
#' via early detection, the subgroup (e.g., tumor type) remains constant.
#'
#' @param x Data frame of class naturalhist (see compile_naturalhist)
#' @return Matrix that maps which rows of the naturalhist data frame 
#' correspond to Early vs Advanced stages
#' 
#' @examples 
#' nh <- compile_naturalhist(prop_adv=0.85, mortrates=c(Early=0.05, Advanced=0.21), 
#'                    subgroup_probs=c(`ER+`=0.5, `ER-`=0.5))
#' create_stageshift_map(nh)
#'
#' @export

create_stageshift_map <- function(x) {

    if (!'naturalhist'%in%class(x)) stop('x must have class naturalhist')

    subgroups <- as.character(unique(x$subgroup))
    stage_pairs <- sapply(subgroups, function(group, df) which(x$subgroup==group), x) 
    rownames(stage_pairs) <- x$stage[stage_pairs[,1]]

    return(stage_pairs)
}

#-------------------------------------------------------------------------------
# stageshift_indicator
#-------------------------------------------------------------------------------
#' Sim indicator of stageshift from binomial
#' 
#' @param x Hazard ratio for early detection (i.e., 1-proportion stage shifted)
#' @param pop_size Population size
#' @param nsim Number of sims
#' @return Matrix of 1s and 0s
#' 
#' @examples
#' #
#'
#' @export

stageshift_indicator <- function(x, pop_size, nsim) {
    return(matrix(rbinom(nsim*pop_size, 1, prob=1-x), nrow=pop_size, ncol=nsim))
}

#-------------------------------------------------------------------------------
# shift_stages
#-------------------------------------------------------------------------------
#' Use indicator of stage shift to change Advanced to Early stage
#' 
#' @param indicator Matrix of 1s and 0s, 1=shift stage
#' @param original Matrix of numeric stage-subgroup id's. Rows=population size, 
#' columns=nsim.
#' @param map Matrix with rows Early and Advanced, columns indicdating subgroups,
#' and cells with stage-subgroup id's (see create_stageshift_map)
#' @param Matrix of new stage-subgroup id's after shift
#' 
#' @examples
#' mapn <- matrix(1:4, ncol=2, dimnames=list(c('Early', 'Advanced'), c('ER+', 'ER-')))
#' orig <- sim_multinom(10000, 2, c(0.25, 0.25, 0.25, 0.25), 1:4)
#' shifti <- stageshift_indicator(0.85, 10000, 2)
#' new <- shift_stages(shifti, orig, mapn)
#' table(new)/table(orig)
#'
#' @export

shift_stages <- function(indicator, original, map) {
    new <- original
    for (s in colnames(map)) {
        adv_cases <- original==map['Advanced',s]
        new[adv_cases & indicator==1] <- map['Early',s]
    }
    return(new)
}

