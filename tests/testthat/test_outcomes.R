
################################################################################
# Tests for treatment.R
context('Modeling mortality outcomes')
################################################################################

test_that('rexp_matrix produces expected means',
          {
              # Set up a large and small rate, each twice
              smallrate <- .01
              bigrate <- .5
              means <- c(c(1,1)/smallrate, c(1,1)/bigrate)

              # Sim 100,000 draws for each rate using matrix input
              rows <- 100000
              m <- cbind(matrix(smallrate, nrow=rows, ncol=2),
                         matrix(bigrate, nrow=rows, ncol=2))
              set.seed(98103)
              rmeans <- colMeans(rexp_matrix(m))
            
              # Check relative error
              error <- abs(rmeans-means)/means
              expect_equal(sum(error>.01), 0)
          }
)

test_that('Mortality sim works', 
          { 
              library(bcimodel)
              set.seed(98103)
              data(ex1) 

              # Small example
              popsize <- 10
              sims <- 5
              # Base case has equally distributed groups 1 to 4 
              b <- matrix(sample.int(4, size=2000, replace=TRUE), 
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

              # Finalize treatments by pairing the non-stage shifted scenario
              # with the stage-shifted cases for the early detection scenario
              tfinal <- update_treat_stageshift(ex1$pol, shifts=s, treats=t)

              # Mortality rate
              mortrates <- lapply(n, return_value_from_id,
                                  df=ex1$nh, value='mortrate')

              #### TEST ONE: mortrates should be the same for scenario 1 and 2
              expect_equal(mortrates[[1]], mortrates[[2]])

              #### TEST TWO: mortrates should be different between 2 and 3 
              #### for shifted cases only
              expect_equal(mortrates[[2]][!st[[3]]], mortrates[[3]][!st[[3]]])

              #### TEST THREE: mortrates for shifted cases should all be 
              #### the early mortrate
              expect_equal(sum(mortrates[[3]][st[[3]]]!=
                              subset(ex1$nh, stage=='Early')$mortrate[1]), 0)

              #  Hazard ratios
              HRs <- lapply(tfinal, return_value_from_id, df=ex1$tx, value='txHR')

              # Final mortality rate
              finalmortrates <- lapply(ex1$pol$num, 
                                       FUN=function(x, HR, rate) {
                                           HR[[x]]*rate[[x]]
                                       },
                                       HR=HRs, rate=mortrates)

              # Times to cancer death, no early detection
              cdt <- timetocancerdeath_by_policy(policies=ex1$pol,
                                                   finalmortrates, popsize, sims)

              # Update with for paired early detection
              cdtfinal <- update_time_stageshift(policies=ex1$pol,
                                                rates=finalmortrates, 
                                                shifts=st, times=cdt)
              #### TEST FOUR: paired times are paired
              ratio <- cdtfinal[[3]]/cdtfinal[[2]]
              # If no shift, ratio is 1
              expect_equal(sum(ratio[!st[[3]]]), sum(!st[[3]]))
              # If shift, ratio is something else that's constant
              # (I'm not sure if that's actually supposed to be,
              # but it seems reasonable?)
              expect_true(length(unique(round(ratio[st[[3]]], 2)))==1)

              # Age at cancer death - create ageentry and ageOC
              ageentry <- matrix(50, nrow=popsize, ncol=sims)
              ageclin <- matrix(60, nrow=popsize, ncol=sims)
              ageOC <- matrix(70, nrow=popsize, ncol=sims)
              ageCD <- lapply(ex1$pol$num,
                            FUN=function(x, ageinc, timetocd) {
                                ageinc[[x]]+timetocd[[x]]
                            }, ageinc=ageclin, timetocd=cdtfinal)
              # Time from study start to cancer death
                timetoCD <- lapply(ex1$pol$num,
                            FUN=function(x, cancer, entry) {
                                cancer[[x]]-entry[[x]]
                            }, cancer=ageCD, entry=ageentry)
              # Time to incidence
                timetoInc <- ageclin-ageentry
                # Cause of death
                cancerD <- lapply(ex1$pol$num,
                            FUN=function(x, ageCD, ageOC) {
                                ifelse(ageCD[[x]]<ageOC[[x]],TRUE,FALSE)
                            }, ageCD, ageOC)

                # Time from trial start to all-cause death
                timetoD <- lapply(ex1$pol$num,
                            FUN=function(x, ageCD, ageOC, ageentry) {
                                ifelse(ageCD[[x]]<ageOC[[x]],
                                       ageCD[[x]]-ageentry[[x]],
                                       ageOC[[x]]-ageentry[[x]])
                            }, ageCD, ageOC, ageentry)

                #### TURN INTO TESTS
                # lapply(ageCD, function(x) x>60)
                # lapply(timetoCD, function(x) x>10)

                # Outcomes - not summarized yet, but returned for each sim
                futimes <- c(5, 15)
                cumulative_incidence <- cuminc(futimes, timetoInc)
                cumulative_mortality <- lapply(ex1$pol$num,
                                               function(x) {
                                                   cummort(futimes, timetoD[[x]],
                                                           cancerD[[x]])
                                               })
                cummortandinc <- cumulative_mortality
                cummortandinc[[length(cumulative_mortality)+1]] <- 
                    cumulative_incidence
                cumulative_yearslived <- lapply(ex1$pol$num,
                                               function(x) {
                                                   cumyears(futimes, timetoD[[x]])
                                               })
                mrr <- calcmrr(cumulative_mortality, 1)
                arr <- calcarr(cumulative_mortality, 1)
                percsurvival <- calcmrr(cummortandinc, 
                                    length(cumulative_mortality)+1,
                                    reverse=TRUE, perc=TRUE)[1:length(cumulative_mortality)]
                yearssaved <- calcarr(cumulative_yearslived, 1, reverse=TRUE)

                # Have a scaled version of survival, so it comes out correctly 
                # after being run through compile_outcomes
                percsurvivalscaled <- lapply(percsurvival, 
                                             function(x) x*popsize/100000)

                # Summarized! - Mean only is the default
                table <- compile_outcomes(list(
                                            `Cumulative Incidence`=cumulative_incidence,
                                            `Cumulative Mortality`=cumulative_mortality,
                                            `% Incident Surviving`=percsurvivalscaled,
                                            `MRR`=mrr, `ARR`=arr, 
                                            `Years of Life Saved`=yearssaved),
                                          futimes,
                                          policynames=ex1$pol$name,
                                          pop_size=popsize)
                # Return lower and upper tables, too
                tables <- compile_outcomes(list(
                                            `Cumulative Incidence`=cumulative_incidence,
                                            `Cumulative Mortality`=cumulative_mortality,
                                            `% Incident Surviving`=survival,
                                            `MRR`=mrr, `ARR`=arr, 
                                            `Years of Life Saved`=yearssaved),
                                          futimes,
                                          policynames=ex1$pol$name,
                                          pop_size=popsize,
                                          stats=c('mean', 'lower', 'upper'))

                # Show a single uncertainty table
                bounds <- format_bounds_list(tables,
                                             digits=c(0,0,1,2,1,0))

        }
)


test_that('compile_outcomes works with only 2 policies',
          { 
            uganda_stdpop <- simpolicies(ex1$pol[1:2,], ex1$nh, ex1$tx,
                                         sims=2, popsize=100)

})

test_that('returnstat4matrix works', {
              mat <- replicate(3, c(0, 0.5, 1))
              # Mean
              expect_equal(unique(returnstat4matrix(mat, 'mean')), 0.5)
              # Continuous quantile type
              expect_equal(unique(returnstat4matrix(mat, 'lower')), 0.025)
              expect_equal(unique(returnstat4matrix(mat, 'upper')), 0.975)
              # Empirical quantile type
              expect_equal(unique(returnstat4matrix(mat, 'lower', quanttype=1)), 0)
              expect_equal(unique(returnstat4matrix(mat, 'upper', quanttype=1)), 1)

})

