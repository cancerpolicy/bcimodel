
################################################################################
# Simulate a series of policies using parameters representative of 
# East Africa
################################################################################


#-------------------------------------------------------------------------------
# Setup 
#-------------------------------------------------------------------------------

library(plyr)
library(bcimodel)
data(ex1)

#-------------------------------------------------------------------------------
# Inputs
#-------------------------------------------------------------------------------
# Stats for manuscript east africa
eafr <- vector('list', length=length(ex1))

scenarios  <- data.frame(num=c(1:10),
                         id=c('adv78', 
                              'adv78.tam', 'adv78.tamchemo', 'adv78.tamchemoERpos',
                              'adv60.tam', 'adv60.tamchemo', 'adv60.tamchemoERpos',
                              'adv35.tam', 'adv35.tamchemo', 'adv35.tamchemoERpos'
                              ),
                         name=c('T0: Surgery only', 'T1: Tamoxifen for ER+', 
                                'T2: T1 plus Chemo for ER-',
                                'T3: T2 plus Chemo for Advanced ER+', 
                                'T1 plus downstaging to 60% advanced',
                                'T2 plus downstaging to 60% advanced',
                                'T3 plus downstaging to 60% advanced',
                                'T1 plus downstaging to 35% advanced',
                                'T2 plus downstaging to 35% advanced',
                                'T3 plus downstaging to 35% advanced'),
                         pairnum=c(NA, NA, NA, NA, rep(c(2,3,4), 2)),
                         earlydetHR=c(rep(1, 4), rep(0.60/0.78, 3),
                                      rep(0.35/0.78, 3)),
                         stringsAsFactors=FALSE)

eafr$nh <- compile_naturalhist(prop_adv=0.78, 
                 mortrates=c(Early=0.0446, Advanced=0.21),
                 subgroup_probs=c(`ER+`=0.41, `ER-`=0.59))

eafr$map <- create_stageshift_map(eafr$nh)

eafr$tx <- data.frame(expand.grid(txSSid=c('None', 'Tamoxifen', 'Chemo',
                                           'Tamoxifen+Chemo'),
                                     SSno=1:nrow(eafr$nh)),
                      stringsAsFactors=FALSE)
ntreat <- nrow(eafr$tx)
ntx <- length(unique(eafr$tx$txSSid))

# Proportions sum to 1 within stage-subgroup groups
eafr$tx <- transform(eafr$tx, 
                     SSid=c(rep('Early.ER+', ntx), rep('Early.ER-', ntx), 
                            rep('Advanced.ER+', ntx), rep('Advanced.ER-', ntx)), 
                     txSSno=1:ntreat)
eafr$tx <- transform(eafr$tx,
                        txHR=c(1, 0.7, 0.775, 0.775*0.7, 
                               1, 1, 0.775, 0.775,
                               1, 0.7,  0.775, 0.775*0.7,
                               1, 1, 0.775, 0.775))
adv78 <- adv78.tam <- rep(0, ntreat)
txSSid <- as.character(eafr$tx$txSSid)
SSid <- as.character(eafr$tx$SSid)
adv78[txSSid=='None'] <- 1 
adv78.tam[txSSid=='Tamoxifen' & (SSid=='Early.ER+' | SSid=='Advanced.ER+')] <- 1
adv78.tamchemo <- adv78.tam
adv78.tamchemo[txSSid=='Chemo' & (SSid=='Early.ER-' | SSid=='Advanced.ER-')] <- 1
adv78.tamchemoERpos <- adv78.tamchemo
adv78.tamchemoERpos[txSSid=='Tamoxifen' & SSid=='Advanced.ER+'] <- 0
adv78.tamchemoERpos[txSSid=='Tamoxifen+Chemo' & SSid=='Advanced.ER+'] <- 1

props <- data.frame(adv78, 
               adv78.tam, adv78.tamchemo, adv78.tamchemoERpos,
               adv78.tam, adv78.tamchemo, adv78.tamchemoERpos,
               adv78.tam, adv78.tamchemo, adv78.tamchemoERpos)
colnames(props) <- scenarios$id

eafr$tx <- data.frame(eafr$tx, props, stringsAsFactors=FALSE)
if (1==0) {
    subset(eafr$tx, SSno==1)
    subset(eafr$tx, SSno==2)
    subset(eafr$tx, SSno==3)
    subset(eafr$tx, SSno==4)
}

eafr$tx <- eafr$tx[,c('SSno', 'SSid', 'txSSno', 'txSSid', 'txHR', 'base', 
                          'tam', 'tamandshift')]
eafr$tx <- data.frame(eafr$tx, stringsAsFactors=FALSE)

#-------------------------------------------------------------------------------
# Plot incidence and mortality rates for Uganda
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Simulate policies
#-------------------------------------------------------------------------------

# Using example data
uganda_stdpop <- simpolicies(ex1$pol, ex1$nh, ex1$tx)


# manuscript_eastafrica <- simpolicies(ex1$pol, ex1$nh, ex1$tx)

