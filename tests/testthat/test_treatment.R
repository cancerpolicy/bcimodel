
################################################################################
# Tests for treatment.R
context('Assigning treatment')
################################################################################

test_that('shifttreatment_indicator works', 
          {
            # First has no shift, second has 15% shift (0.85 HR)
            s <- list(matrix(0, nrow=1000, ncol=2),
                           stageshift_indicator(0.85, 1000, 2))
            # Base case has equally distributed groups 1 to 4
            b <- matrix(sample.int(4, size=2000, replace=TRUE), 
                               nrow=1000, ncol=2)
            # Map of 1:4 onto stage and ER status
            m <- matrix(1:4, nrow=2, dimnames=list(c('Early', 'Advanced'), 
                                                     c('ER+', 'ER-')))
            # If type[x]=NA, should return all NA's
            expect_equal(sum(is.na(shifttreatment_indicator(x=1, type=c(NA, 2), 
                                                            s, b, m))), 2000)
            ind <- shifttreatment_indicator(x=2, type=c(NA, 2), s, b, m)
            expect_true(abs(round(mean(ind)-(0.5*0.15),2))<=0.02)
          }
)
test_that('sim_treatment_by_subgroup works',
          { 
              library(bcimodel)
              set.seed(98103)
              data(ex1) 

              # Small example
              popsize <- 1000
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
              t <- sim_treatment_by_subgroup(ex1[[4]], n[[1]], 'base', popsize, nsim)
          }
)

test_that('treatments_by_policy and update_treat_stageshift work', 
          { 
              library(bcimodel)
              set.seed(98103)
              data(ex1) 

              # Small example
              popsize <- 10
              sims <- 5
              # Base case has equally distributed groups 1 to 4 
              b <- matrix(sample.int(4, size=popsize*sims, replace=TRUE), 
                          nrow=popsize, ncol=sims) 
              # Map of 1:4 onto stage and ER status 
              m <- ex1$map
              # Shifts
              s <- lapply(ex1$pol$earlydetHR, 
                          stageshift_indicator, pop_size=popsize, nsim=sims) 
              # Get new stages (advanced cases only)
              n <- lapply(s, shift_stages, original=b, map=m)
              # Create indicator for shifting treatment (advanced cases only)
              st <- lapply(ex1$pol$num, shifttreatment_indicator, 
                           type=ex1$pol$pairnum, shifts=s, basecase=b, map=m)

              # Simulate treatment (for early detection scenarios, candidate
              # early-stage treatments for shifted cases)
              t <- treatments_by_policy(policies=ex1[[1]], treat_chars=ex1[[4]], 
                                        stagegroups=n, map=m, popsize, sims)

              ####### TEST ONE - TO DO
              # Scenarios with early detection should only have early-stage
              # treatments

              ####### TEST TWO - TO DO
              # Scenarios with early detection should only have early-stage 

              ####### TEST THREE
              # Shift treatment indicator is TRUE only for cases that were
              # advanced-stage in the base case
              expect_equal(sum(!b[st[[3]]]%in%m['Advanced',]),0)

              ####### TEST FOUR
              # In paired non-earlydet scenario, shift treatment indicator
              # is TRUE only for advanced-stage treatments
              expect_equal(sum(!t[[2]][st[[3]]]%in%
                           subset(ex1[[4]], SSno%in%m['Advanced',])$txSSno), 0)
              
              ####### TEST FIVE
              # New stages for the shifted cases are early stages
              expect_equal(
                           sum(!n[[3]][st[[3]]]%in%m['Early',]), 0
                           )

              ####### TEST SIX (similar to TEST FOUR)
              # In paired non-earlydet scenario, shift treatment indicator
              # is TRUE only for advanced-stage treatments
              expect_equal(
                           sum(!t[[3]][st[[3]]]%in%
                           subset(ex1[[4]], SSno%in%m['Early',])$txSSno), 0
                           )

              ####### TEST SEVEN (similar to TEST FOUR)
              # Scenario 3 is where txSSno 1 and 4 have prop=0, so we 
              # should only see 2 and 3
              expect_equal(
                           sum(!t[[3]][st[[3]]]%in%c(2,3)), 0
                           )

              # Finalize treatments by pairing the non-stage shifted scenario
              # with the stage-shifted cases for the early detection scenario
              tfinal <- update_treat_stageshift(ex1$pol, shifts=s, treats=t)

              ####### TEST EIGHT
              # Treatments are the same between scenario 2 and 3 for 
              # non-shifted cases
              expect_equal(tfinal[[2]][!s[[3]]], tfinal[[3]][!s[[3]]])

              ####### TEST NINE
              # Treatments are only early-stage for shifted cases in final
              expect_equal(
                           sum(!tfinal[[3]][s[[3]]]%in%c(2,3)), 0
                           )
            }
)


test_that('sim_treatment_by_subgroup works if only 1 treatment', {

    # Set up 1-treatment scenario
    library(bcimodel)
    data(ex1)
    ex1$tx <- subset(ex1$tx, txSSid=='None')
    ex1$tx <- transform(ex1$tx, txSSno=1:4, base=1, tam=1, tamandshift=1)

    # Small example
    popsize <- 10
    sims <- 5
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
    t <- treatments_by_policy(policies=ex1[[1]], treat_chars=ex1[[4]], 
                            stagegroups=n, map=m, popsize, sims)

    # MANUALLY RERAN TESTS FROM ABOVE...

}

