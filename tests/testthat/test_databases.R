################################################################################
# Tests for incidence and lifetable
context('Incidence and lifetable database construction')
################################################################################

test_that('Rates transform to cumulative survival', 
          {
              data(incratesf)
              data(allmortratesf)

              inc <- interpolate_cumsurv(incratesf, 
                                         ratevar='Female.Rate.Per.100K',
                                         country='Uganda')
              mort <- interpolate_cumsurv(allmortratesf, 
                                         ratevar='Rate.Per.100K',
                                         country='Uganda')

              # Some clever test....
          }
)


test_that('format_age works for different age limits', 
          {
              data(agestructure)
              thisage <- subset(agestructure, Country=='Standard')

              # No age limits
              expect_identical(thisage, 
                               format_age(thisage, minAge=0, maxAge=100))

              # Some age limits
              thispop <- sum(subset(thisage, age>=50)$pop)
              expect_identical(thispop, 
                               sum(format_age(thisage, 
                                              minAge=50, maxAge=100)$pop))

              # Proportion sums to 1?
              expect_equal(1, sum(format_age(thisage, minAge=50,
                                             maxAge=100)$prop))
              expect_equal(1, sum(format_age(thisage, minAge=20,
                                             maxAge=39)$prop))
          }
)
