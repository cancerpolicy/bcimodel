## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(echo = TRUE)

## ----cars----------------------------------------------------------------
library(bcimodel)
data(ex1)

## ------------------------------------------------------------------------
ex1$pol

## ------------------------------------------------------------------------
ex1$nh

## ------------------------------------------------------------------------
# 85% are advanced stage
# Mortality rates are .05 for Early and 0.21 for Advanced
# 50% are ER+
(nh1 <- compile_naturalhist(prop_adv=0.85, 
                            mortrates=c(Early=0.05, Advanced=0.21), 
                            subgroup_probs=c(`ER+`=0.5, `ER-`=0.5)))

## ------------------------------------------------------------------------
# 85% are advanced stage
# Mortality rates are .05 for Early and 0.21 for Advanced
# 50% are ER+, 15% are HER2+, so uncorrelated gives 0.5*0.15 = .075 ER+HER2+
(nh2 <- compile_naturalhist(prop_adv=0.85, 
                            mortrates=c(Early=0.05, Advanced=0.21),
                            subgroup_probs=c(`ER+HER2+`=0.075,
                                             `ER-HER2+`=0.075,
                                             `ER+HER2-`=0.425,
                                             `ER-HER2-`=0.425)))

## ---- eval=FALSE---------------------------------------------------------
#  class(your_nh_data_frame) <- append(class(your_nh_data_frame), "naturalhist")

## ------------------------------------------------------------------------
ex1$map

## ------------------------------------------------------------------------
create_stageshift_map(nh2)

## ------------------------------------------------------------------------
ex1$tx

## ------------------------------------------------------------------------
set.seed(98103)
uganda_stdpop <- bcimodel::simpolicies(ex1$pol, ex1$nh, ex1$tx, popsize=10000)

## ------------------------------------------------------------------------
knitr::kable(uganda_stdpop[['5']]) # could also specify using uganda_stdpop$`5`

## ---- eval=FALSE---------------------------------------------------------
#  # Not run in vignette - three policies, 100 sims, with pop size 100,000
#  set.seed(98103)
#  uganda_stdpop <- bcimodel::simpolicies(ex1$pol, ex1$nh, ex1$tx, ncores=4)

