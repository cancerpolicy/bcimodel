
################################################################################
# Tests for initialize.R
context('General worker functions')
################################################################################

test_that('sim_multinom works', 
    {
        probs <- c(0.1, 0.3, 0.6)
        simdat <- sim_multinom(100, 10, probs, names=c('a', 'b', 'c'))
        absdiff <- abs(table(simdat)/(nrow(simdat)*ncol(simdat))-probs)
        expect_true(mean(absdiff)<=0.01)
    }
)
