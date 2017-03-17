
################################################################################
# Preprocessing for incidence database
################################################################################

#-------------------------------------------------------------------------------
# Create matrix of incidence file names and descriptions
#-------------------------------------------------------------------------------

# Read in hand-prepped file based on the CI5-X registry guide
incfiles_loc <- system.file('extdata/incidence/CI5-Xd/files_to_select.csv', 
                       package='bcimodel')

incfiles <- read.csv(incfiles_loc)

# Set old and new folder location of the raw CI5-X registry files
old_loc <- '~/Dropbox/BCModel/data/incidence/incidence_CI5-Xd'
new_loc <- '~/Documents/jbirnbau/bcimodel/inst/extdata/incidence/CI5-Xd'

# According to files in guide, extract and rename files and save them into
# new_loc folder

for (i in 1:nrow(incfiles)) {
    f <- incfiles$regno[i]
    newname <- incfiles$filen[i]
    file.copy(file.path(old_loc, paste0(f, '.csv')), file.path(new_loc, newname))
}


#-------------------------------------------------------------------------------
# Now run data-raw/incidence.R
#-------------------------------------------------------------------------------
