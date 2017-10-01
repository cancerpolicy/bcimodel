
################################################################################
# Functions to format data for the incidence, lifetable and age structure
# databases embedded in the package
################################################################################

# See the data creation files in data-raw

#-------------------------------------------------------------------------------
# bcirates
#-------------------------------------------------------------------------------
#' Format female breast cancer incidence rates from SEER or IARC/GLOBOCAN
#' Specifies age by the midpoint of the age category
#' @param incfile Path to file
#' @param source One of 'seer' or 'globocan'
#' @return Data frame
#' @examples
#' IARC/Globocan data
#' bcirates(system.file('extdata', 
#'                     'incidence/CI5-Xd/18000299_uga.csv', package='bcimodel'),
#'         source='globocan')
#'
#' IARC/Globocan data
#' bcirates(system.file('extdata', 
#'                     'incidence/SEER/bc_1975-1979_incidence.csv', 
#'                      package='bcimodel'),
#'         source='globocan')
#'
#' @export

bcirates <- function(incfile, source='seer') {
    switch(source, 
           seer = {
                # Read in data
                incdat <- read.csv(incfile, stringsAsFactors=FALSE)
                incdat <- within(incdat, {
                               Age <- gsub('[0-9]+=', '', Ages)
                               Age <- as.character(gsub(' years', '', Age))
                               Age[Age=='85+'] <- 87
                               Age <- suppressWarnings(sapply(strsplit(Age, '-'),
                                             function(x) mean(as.numeric(x))))
                               Cases=Count
                               Female.Rate.Per.100K=Crude.Rate
                        })
                incdat <- incdat[!is.nan(incdat$Age),]
                incdat <- incdat[!is.na(incdat$Age),]

           },
           globocan = {
                # Read in data
                incdat <- read.csv(incfile, stringsAsFactors=FALSE, header=FALSE)
                # Add the source column if it doesn't exist
                if (ncol(incdat)<=5) {
                    incdat$source=''
                    incdat$source[1] <- "http://ci5.iarc.fr/CI5I-X/Pages/download.aspx"
                }
                colnames(incdat)  <- c('sex', 'site', 'ageGroup', 
                                       'cases', 'persYears', 'source')
                # Select female breast cancer, construct rate and 
                # Age (transform 1-19 into age lower bound; ignore 19=unknown age)
                incdat <- subset(incdat, site==113 & sex==2 & ageGroup<19)
                incdat <- within(incdat, {
                                    Female.Rate.Per.100K=100000*(cases/persYears)
                                    Age=(ageGroup-1)*5 + 2
                                    Age[Age==87.5] <- 87
                                    Cases=cases
                                })

                # Simply reconcile cases but no person-years in older ages:
                # hold the rate constant at the last observed rate
                zero <- which(incdat$persYears==0)
                incdat[zero,'Female.Rate.Per.100K'] <- 
                    incdat[min(zero)-1,'Female.Rate.Per.100K']
                warning('Graph age-specific incidence to check assumptions')
           }
        )

    # Min and max age - include 0 rate at age 0
    if (min(incdat$Age)!=0) {
        if (length(unique(incdat$Age))!=nrow(incdat)) {
            stop('More than one group (sex or country) not allowed in bcirates')
        }
        incdat <- rbind(transform(incdat[1,], Age=0, Cases=0, 
                                  Female.Rate.Per.100K=0),
                        incdat)
        if ('Crude.Rate'%in%colnames(incdat)) {
            incdat$Crude.Rate[incdat$Age==0] <- 0
        }
    }

    return(incdat)
}

#-------------------------------------------------------------------------------
# allmortrates
#-------------------------------------------------------------------------------
#' Format life table data from IHME or other source 
#' 
#' @param lifefile Data frame with lifetable data
#' @param source Data source, only 'ihme' supported right now
#' @return Formatted all-cause mortality rates by age
#'
#' @export

allmortrates <- function(lifefile, source='ihme') {

    # Format into binned rates, for input to cumsurv_by_country
    switch(source,
           ihme = {
                lt <- read.csv(lifefile)
                lt <- subset(lt, year==max(lt$year) & metric=='Deaths' & 
                                   age_group_name!='All Ages' & 
                                   sex_name!='Both sexes',
                              select=c(location_name, age_group_name, 
                                       sex_id, year, mean))
                lt <- within(lt, {
                                  Country=location_name
                                  Year=year
                                  Age=gsub(" to [0-9]*", '', age_group_name)
                                  Age[Age=="<1 year"] <- '0'
                                  Age[Age=="80 plus"] <- '85'
                                  Age=as.numeric(Age) + 2
                                  Age[Age==2] <- 0
                                  Age[Age==3] <- 2.5
                                  Male=sex_id
                                  Male[sex_id==2] <- 0
                                  Rate.Per.100K=mean
                              })
            })

    # Order
    lt <- lt[order(lt$Country, lt$Male, lt$Age),]
    lt <- lt[,c('Country', 'Year', 'Male', 'Age', 
                'Rate.Per.100K', 'age_group_name')]

    return(lt)
}

#-------------------------------------------------------------------------------
# plot.rates
#-------------------------------------------------------------------------------
#' Plot the data in the incidence database
#'
#' @param x Probably age group indicator (must be numeric)
#' @param y Probably age-specific rates (per 100,000)
#' @param psize Variable to scale point size
#' @param group Variable indicating groups to plot using different colors
#' @return ggplot
#' @examples
#' # Plot incidence
#' data(incratesf)
#' plot_rates(incratesf$Age, incratesf$Female.Rate.Per.100K,
#'         psize=incratesf$Cases, group=paste(incratesf$Country,
#'                                            incratesf$Year))
#' 
#' # Plot mortality
#' data(allmortratesf)
#' df <- subset(allmortratesf, 
#'                   Country=='United States' | Country=='Uganda')
#' plot_rates(df$Age, df$Rate.Per.100K, psize=NULL, group=df$Country)
#'
#' @export

plot_rates <- function(x, y, psize=NULL, group=NULL) {
    if (length(x)!=length(y)) stop('x and y must be same length')
    if (is.null(psize)) psize <- rep(1, length(x))
    if (is.null(group)) group <- rep(1, length(x))

    df <- data.frame(x, y, psize, group)

    p <- ggplot(df, aes(x,y,group=group)) + geom_line(aes(color=group)) + 
        geom_point(aes(size=psize, color=group)) +
        scale_x_continuous(name='Age Group') +
        scale_y_continuous(name='Rate per 100,000') +
        scale_size(name='Number of Cases') + 
        scale_color_discrete(name='Population') + 
        theme_bw() 
    return(p)
}

#-------------------------------------------------------------------------------
# running_product
#-------------------------------------------------------------------------------
#' Return running product of a vector, for computing cumulative survival
#' @param x Vector of conditional survival probabilities
#' @return Vector of running products
#'
#' @export
running_product <- function(x) sapply(1:length(x), FUN=function(y) prod(x[1:y]))

#-------------------------------------------------------------------------------
# rate_to_cumsurv
#-------------------------------------------------------------------------------
#' Given event rates, compute cumulative survival
#' @param x Vector of event rates (interpreted as discrete rates, not instantaneous hazards)
#' @return Vector of cumulative survivals
#'
#' @export
rate_to_cumsurv <- function(x) {
    conditional_surv <- 1-x
    running_product(conditional_surv)
}

#-------------------------------------------------------------------------------
# interpolate_cumsurv
#-------------------------------------------------------------------------------
#' Compute cumulative event-free survival curve for single-year ages
#' 
#' Takes age-binned event rates per 100,000 and uses interpolation to estimate
#' cumulative event-free survival for single years
#' 
#' @param binnedrates Data frame with Age for age groups 
#' @param ratevar Name of rate variable
#' @param country Country to select from the data frame
#' @param maxage Maximum age desired. If data are incomplete, the 
#' rate for the oldest age will be held constant until maxage
#' @return Data frame of cumulative survival probability by single-year ages
#' @examples
#' data(incratesf)
#' data(allmortratesf)
#' 
#' inc <- interpolate_cumsurv(incratesf, 
#'                          ratevar='Female.Rate.Per.100K',
#'                          country='Uganda')
#' mort <- interpolate_cumsurv(allmortratesf, 
#'                          ratevar='Rate.Per.100K',
#'                          country='Uganda')
#'
#' @export

interpolate_cumsurv <- function(binnedrates, ratevar, maxage=NULL,
                                country=NULL) {

    if (!is.null(country)) binnedrates <- subset(binnedrates, 
                                                 Country==country)
    age_max <- max(binnedrates$Age)

    # Interpolate for single year, converting from rate to probability
    singleyears <- 
        data.frame(age = 0:age_max, 
                   singleyears_rate = approx(binnedrates$Age, 
                                             binnedrates[,ratevar], 
                                                 xout=0:age_max)$y/100000)

    if (!is.null(maxage)) {
        # Hold oldest rate constant until max age
        if (maxage>age_max) {
            addyears <- data.frame(age=(age_max+1):maxage, 
                                   singleyears_rate=
                                       subset(singleyears, 
                                              age==age_max)$singleyears_rate)
            singleyears <- rbind(singleyears, addyears)
        }
    }
    
    singleyears <- transform(singleyears,
                           cumsurv=rate_to_cumsurv(singleyears_rate))

    return(singleyears)
}


#-------------------------------------------------------------------------------
# format_age
#-------------------------------------------------------------------------------
#' Convert age-group population sizes to population proportions
#' 
#' Takes population counts by age group and returns the population proportions
#' for single-year age groups
#' 
#' @param age_data Population age structure data frame OR file
#' @param gsize Number of years covered by each age group
#' @param minAge Optional lower age limit for the population
#' @param maxAge Optional upper age limit for the population
#' @param format Does the input need to be formatted? If so, 
#' age_data can be a .csv file
#' @note See data-raw/agestructure.R
#' 
#' @export

format_age <- function(age_data, gsize=5, minAge=20, maxAge=85, 
                       format=FALSE) {

    if (format) {
        grouped <- read.csv(age_data, header=TRUE)

        if (!is.numeric(grouped$Population)) {
            stop('Error in format_age: character pop sizes')
        }

        # Split into single-year age groups
        split <- lapply(grouped$Age, function(a,data=grouped,g=gsize) {
                            NewPop <- data[which(data$Age==a), 'Population']/g
                            data.frame(age=a:(a+4), pop=rep(NewPop,g))
                             })
        age_data <- do.call('rbind', split)
        
        # Verify
        if (sum(age_data$pop)!=sum(grouped$Population)) {
            stop('Division error in format_age')
        }
    }

    # Return proportions after limiting the minimum age
    age_data <- subset(age_data, age>=minAge & age<=maxAge)
    age_data <- transform(age_data, prop=pop/sum(pop))

    return(age_data)

}


