################################################################################
# allmortratesf: Create cumulative survival database from lifetables
################################################################################

#-------------------------------------------------------------------------------
# Specify countries to include
#-------------------------------------------------------------------------------

# Files stored in inst/extdata/incidence/CI5-Xd to access
# 1/10/17: We start with only Uganda, with plans to expand

countries <- c('United States', 'Uganda')

#-------------------------------------------------------------------------------
# Setup
#-------------------------------------------------------------------------------
library(devtools)
library(bcimodel)

#-------------------------------------------------------------------------------
# Unzip files into working dir
#-------------------------------------------------------------------------------
ihme_file <- system.file('extdata', 'mortality/IHME_lifet/lifetable_IHME_GBD_2013_LIFE_TABLE_1990_2013_Y2014M12D17.zip', package='bcimodel')
unzip(ihme_file, overwrite=TRUE)

#-------------------------------------------------------------------------------
# Format
#-------------------------------------------------------------------------------
file <- 'IHME_GBD_2013_LIFE_TABLE_1990_2013_Y2014M12D17.CSV'

allmortrates <- allmortrates(file)

#-------------------------------------------------------------------------------
# Delete unzipped files
#-------------------------------------------------------------------------------
file.remove(grep('IHME_GBD_2013', dir(), value=TRUE))

#-------------------------------------------------------------------------------
# Compile and save
#-------------------------------------------------------------------------------

allmortratesf <- subset(allmortrates, Male==0)
use_data(allmortratesf, overwrite=TRUE)


################################################################################
# cohortltf: Load cohort mortality from BMD, used in cantrance/Annals
################################################################################

library(cantrance)
data(life_table)
cohortltf <- subset(life_table, Male==0)
use_data(cohortltf, overwrite=TRUE)


