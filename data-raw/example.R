
################################################################################
# Create inputs for example
################################################################################

library(plyr)
#-------------------------------------------------------------------------------
# Scenarios with early detection hazard ratios
#-------------------------------------------------------------------------------

# Three scenarios: base, tamoxifen, and tamoxifen with stage shift

# Pairnum indicates, for scenarios with stage shifts, which other scenario
# is the non-stage-shifted pair. Pairs should have the same treatment 
# proportions in treatinfo below. Set pairnum to NA for the base case scenario
# and for scenarios with no early detection, only treatment changes
scenarios  <- data.frame(num=c(1:3),
                         id=c('base', 'tam', 'tamandshift'),
                         name=c('Base Case', 'Tamoxifen for ER+', 
                                'Tamoxifen for ER+ and 30% stage shift'),
                         pairnum=c(NA, NA, 2),
                         earlydetHR=c(1, 1, 0.70),
                         stringsAsFactors=FALSE)

#-------------------------------------------------------------------------------
# Natural history
#-------------------------------------------------------------------------------

naturalhist <- compile_naturalhist(prop_adv=0.85, 
                 mortrates=c(Early=0.0446, Advanced=0.21),
                 subgroup_probs=c(`ER+`=0.3, `ER-`=0.7))

stagepairs <- create_stageshift_map(naturalhist)

#-------------------------------------------------------------------------------
# Treatments
#-------------------------------------------------------------------------------
treatinfo <- data.frame(expand.grid(txSSid=c('None', 'Tamoxifen'),
                                     SSno=as.numeric(rownames(naturalhist))
                                     ))
ntreat <- nrow(treatinfo)

# Proportions sum to 1 within stage-subgroup groups
treatinfo <- transform(treatinfo,
                       SSid=c(rep('Early.ER+', 2), rep('Early.ER-', 2),
                              rep('Advanced.ER+', 2), rep('Advanced.ER-', 2)),
                        txSSno=1:ntreat,
                        txHR=c(1.0, 0.7, 1, 1, 1, 0.7, 1, 1),
                        base=c(rep(c(0.8, 0.2),4)),
                        tam=c(c(0,1), c(1,0), c(0,1), c(1,0)),
                        tamandshift=c(c(0,1), c(1,0), c(0,1), c(1,0)))

treatinfo <- treatinfo[,c('SSno', 'SSid', 'txSSno', 'txSSid', 'txHR', 'base', 
                          'tam', 'tamandshift')]

#-------------------------------------------------------------------------------
# Compile and save
#-------------------------------------------------------------------------------

ex1 <- list(pol=scenarios, nh=naturalhist, map=stagepairs, tx=treatinfo)
use_data(ex1, overwrite=TRUE)

