
################################################################################
# Tests for initialize.R
context('Initializing population')
################################################################################

test_that('calc_ac_lifespan requirements', 
          {
            # Make some life table data
            data(allmortratesf)
            lt <- transform(interpolate_cumsurv(allmortratesf, 
                                      ratevar='Rate.Per.100K',
                                      country='Tanzania'), 
                            birth_year=1967, BirthCohort=1967, 
                            Age=age, Survival=cumsurv, Male=0)

            # Age at OC death
            ages <- calc_ac_lifespan(age=50, male=0, lifetable=lt)

            # No NA's
            expect_true(sum(is.na(ages))==0)

            # Ages at OC death greater than ages at entry
            expect_true(sum(ages<50)==0)
          }
)
test_that('calc_ac_lifespan_pop works',
          {
            library(plyr)
            # Make some life table data
            data(allmortratesf)
            lt <- transform(interpolate_cumsurv(allmortratesf, 
                                      ratevar='Rate.Per.100K',
                                      country='Tanzania'), 
                            birth_year=1967, BirthCohort=1967, 
                            Age=age, Survival=cumsurv, Male=0)

            # Fake data
            pop <- data.frame(age=seq(10,50,by=10),
                              male=0)
            rows <- replicate(100, 1:nrow(pop))
            ages <- calc_ac_lifespan_pop(pop, rows, lt)

            # No NA's
            expect_true(sum(is.na(ages))==0)

            # Ages at OC death greater than ages at entry
            compare <- ages>=replicate(100, pop$age)
            expect_true(sum(!compare)==0)
          }
)
test_that('Initialization works', 
          {
            pop <- initialize_pop(pop_size=100000,
                                  nsim=2, 
                                  agesource='Standard', 
                                  minage=0, maxage=100, 
                                  incsource='Uganda', 
                                  mortsource='Uganda')
          }
)

test_that('Adding stage-tumor groups returns correct proportions',
          {
              # Proportions sum to 1
              prop_in <- c(0.1, 0.9, 0.5)
              expect_error(add_features(matrix(1, nrow=1, ncol=1),
                                                     probs=prop_in),
                           'Probabilities do not sum to 1')

              # Proportions in yield proportions out, within a 10%
              # margin for 1000*2=2000 random number generations
              prop_in <- c(0.25, 0.25, 0.5)
              distribution <- add_features(matrix(1, nrow=1000, ncol=2), 
                                                 probs=prop_in)
              prop_out <- round(prop.table(table(distribution)),2)
              expect_true(sum(abs(prop_in-prop_out))<=0.1)

              # Works with character names
              distribution <- add_features(matrix(1, nrow=1000, ncol=2), 
                                                 probs=prop_in,
                                                 names=c('G1', 'G2', 'G3'))
              prop_out <- round(prop.table(table(distribution)),2)
              expect_true(sum(abs(prop_in-prop_out))<=0.1)

          }
)

test_that('Natural history parameters are regulated by naturalhist class',
          {
             nh <- compile_naturalhist(prop_adv=0.85, 
                                 mortrates=c(Early=0.05, Advanced=0.21),
                                 subgroup_probs=c(`ER+`=0.5, `ER-`=0.5))
             expect_true('naturalhist'%in%class(nh))

             expect_error(compile_naturalhist(prop_adv=0.85, 
                                 mortrates=c(Early=0.05, Advanced=0.21),
                                 subgroup_probs=c(`ER+`=0.25, `ER-`=0.5)),
                          'Check that subgroup_probs sum to 1')
          })
