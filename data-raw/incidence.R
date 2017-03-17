################################################################################
# Create incidence database
################################################################################

#-------------------------------------------------------------------------------
# Specify countries to include, and their file names
# (Files stored in inst/extdata/incidence/CI5-Xd)
#-------------------------------------------------------------------------------
# 3/17/17 Update: presumes that "incfiles" was created and loaded into memory
# via data-raw/incidence_preprocess.R


# 1/10/17: We start with only Uganda, with plans to expand

#incfiles <- data.frame(filen=c('18000299_uga.csv'),
#                       country='Uganda',
#                       county='Kyadondo County',
#                       years='2003-2007')

#-------------------------------------------------------------------------------
# Setup
#-------------------------------------------------------------------------------
library(devtools)
library(bcimodel)
library(plyr)

#-------------------------------------------------------------------------------
# Create database
#-------------------------------------------------------------------------------


inclist <- lapply(1:nrow(incfiles), FUN=function(x) {
    inctable <- bcirates(system.file('extdata', 
                                     paste0('incidence/CI5-Xd/',
                                            incfiles$filen[x]),
                         package='bcimodel'), source='globocan')
    inctable <- subset(transform(inctable,
                          Country=incfiles$country[x],
                          County=incfiles$county[x],
                          Years=incfiles$years[x]),
                       select=c('Country', 'County', 'Years', 
                                'Age', 'Cases', 'Female.Rate.Per.100K'))
    return(inctable)
})


#-------------------------------------------------------------------------------
# Add SEER
#-------------------------------------------------------------------------------

inclist[[length(inclist)+1]] <- 
    subset(transform(bcirates(system.file('extdata', 
                                'incidence/SEER/bc_1975-1979_incidence.csv', 
                                package='bcimodel'), source='seer'),
                     Country='United States',
                     County='SEER 9',
                     Years='1975-1979'), 
           select=c('Country', 'County', 'Years', 
                    'Age', 'Cases', 'Female.Rate.Per.100K'))



#-------------------------------------------------------------------------------
# Compile and save
#-------------------------------------------------------------------------------

incratesf <- ldply(inclist)
use_data(incratesf, overwrite=TRUE)


