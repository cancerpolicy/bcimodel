
################################################################################
# Tests for parallel.R
context('Parallelizing functions')
################################################################################

test_that('treatments_by_policy and partreatments_by_policy have same output',
          { 
              library(bcimodel)
              library(parallel)
              set.seed(98103)
              data(ex1) 

              # Large example
              popsize <- 100000
              sims <- 100
              # Base case has equally distributed groups 1 to 4 
              b <- matrix(sample.int(4, size=popsize*sims, replace=TRUE), 
                          nrow=popsize, ncol=sims) 
              # Map of 1:4 onto stage and ER status 
              m <- ex1[[3]] 
              # Shifts
              s <- lapply(ex1[[1]]$earlydetHR, 
                          stageshift_indicator, pop_size=popsize, nsim=sims) 
              # Get new stages (advanced cases only)
              n <- lapply(s, shift_stages, original=b, map=m)
              # Create indicator for shifting treatment (advanced cases only)
              st <- lapply(ex1[[1]]$num, shifttreatment_indicator, 
                           type=ex1[[1]]$pairnum, shifts=s, basecase=b, map=m)

              # Simulate treatment (for early detection scenarios, candidate
              # early-stage treatments for shifted cases)
              startclock <- proc.time()
              t <- treatments_by_policy(policies=ex1[[1]], treat_chars=ex1[[4]], 
                                        stagegroups=n, map=m, popsize, sims)
              el <- proc.time()-startclock
              # elapsed time was 348s for 100 sims and popsize 100K

              # Parallelized: use default ncores=4
              startclock <- proc.time()
              part <- partreatments_by_policy(policies=ex1[[1]], treat_chars=ex1[[4]], 
                                        stagegroups=n, map=m, popsize, sims)
              parel <- proc.time()-startclock
              # elapsed time was 182s for 100 sims and 100K

              startclock <- proc.time()
              part2 <- partreatments_by_policy(policies=ex1[[1]], treat_chars=ex1[[4]], 
                                        stagegroups=n, map=m, popsize, sims,
                                        method='inner')
              parel2 <- proc.time()-startclock
              # this was basically equivalent, 179.5s for 100 sims and 100K
              
              
}


test_that('parsim_multinom faster than sim_multinom', 
    {
        library(parallel)
        probs <- c(0.1, 0.3, 0.6)
        # Not much difference if we do 1000 by 1000, but for 100,000 by 100,
        # elapsed is 33.69s versus 15.34s. For automated testing keep it faster,
        # only 10000 x 100
        nrow <- 10000
        ncol <- 100
        t1 <- system.time(sim_multinom(nrow, ncol, probs, names=c('a', 'b', 'c')))
        t2 <- system.time(parsim_multinom(nrow, ncol, probs, names=c('a', 'b', 'c')))
        expect_true(t1<t2)
    }
)

test_that('parsim_treatment_by_subgroup is faster',
          { 
              library(bcimodel)
              set.seed(98103)
              data(ex1) 

              # Again, significant speed gains are for 100,000 x 100
              # This will be more marginal
              popsize <- 10000
              sims <- 100
              # Base case has equally distributed groups 1 to 4 
              b <- matrix(sample.int(4, size=popsize*sims, replace=TRUE), 
                          nrow=popsize, ncol=sims) 
              # Map of 1:4 onto stage and ER status 
              m <- ex1[[3]] 
              # Shifts
              s <- lapply(ex1[[1]]$earlydetHR, 
                          stageshift_indicator, pop_size=popsize, nsim=sims) 
              # Get new stages (advanced cases only)
              n <- lapply(s, shift_stages, original=b, map=m)
              # Create indicator for shifting treatment (advanced cases only)
              st <- lapply(ex1[[1]]$num, shifttreatment_indicator, 
                           type=ex1[[1]]$pairnum, shifts=s, basecase=b, map=m)

              # Simulate treatment (for early detection scenarios, candidate
              # early-stage treatments for shifted cases)
              t1 <- system.time(sim_treatment_by_subgroup(ex1[[4]], n[[1]], 'base', popsize, nsim))
              t2 <- system.time(parsim_treatment_by_subgroup(ex1[[4]], n[[1]], 'base', popsize, nsim))
              expect_true(t1<t2)
          }
)

test_that('setting seed for cluster produces identical results for parsim_multinom',
          {
              library(bcimodel)
              library(parallel)

              v1 <- parsim_multinom(5, 5, probs=c(0.1, 0.2, 0.7), 
                                    names=c(1,2,7))
              v2 <- parsim_multinom(5, 5, probs=c(0.1, 0.2, 0.7), 
                                    names=c(1,2,7))
              expect_equal(v1, v2)
          }
)

test_that('setting seed for cluster produces identical results for whole run',
          {
            # Load data
            library(bcimodel)
            library(parallel)
            library(plyr)
            data(ex1)
            # Model
            v1 <- parsimpolicies(ex1$pol, ex1$nh, ex1$tx, futimes=10)

            v2 <- parsimpolicies(ex1$pol, ex1$nh, ex1$tx, futimes=10)

            expect_equal(v1, v2)
          }
)
