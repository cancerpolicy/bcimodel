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
          }
)
