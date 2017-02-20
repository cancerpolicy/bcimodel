
################################################################################
# Create database of age structures
################################################################################

#-------------------------------------------------------------------------------
# Specify countries to include, and their file names
#-------------------------------------------------------------------------------

# Files stored in inst/extdata/incidence/CI5-Xd to access
# 1/10/17: We start with only Uganda, with plans to expand

files <- data.frame(filen=c('std_age.csv', 
                               'tza_age.csv'),
                       country=c('Standard', 
                                 'Tanzania'),
                       source=c('WHO 2000-2025 Standard via SEER 2013', 
                                'Tanzania Bureau of Statistics 2013'))


#-------------------------------------------------------------------------------
# Setup
#-------------------------------------------------------------------------------
library(devtools)
library(bcimodel)
library(plyr)

#-------------------------------------------------------------------------------
# Create database
#-------------------------------------------------------------------------------


dlist <- lapply(1:nrow(files), FUN=function(x) {
    raw <- system.file('extdata', file.path('agestructure', files$filen[x]),
                         package='bcimodel')
    formatted <- format_age(raw, minAge=0, maxAge=100, format=TRUE)
    formatted <- transform(formatted,
                          Country=files$country[x],
                          Source=files$source[x])
    return(formatted)
})

#-------------------------------------------------------------------------------
# Compile and save
#-------------------------------------------------------------------------------

agestructure <- ldply(dlist)
use_data(agestructure, overwrite=TRUE)

