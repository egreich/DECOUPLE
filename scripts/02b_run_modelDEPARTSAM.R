#### Run model 2

# Set run params
args<-commandArgs(TRUE)
print(args)
print("varname:")
(varname <- as.numeric(args[1]))
print("sitename:")
(sitename <- as.numeric(args[2]))
print("seed:")
(SEED <- as.numeric(args[3]))

# Set defined R seed
set.seed(SEED, kind = NULL, normal.kind = NULL)
# Generate "random" seed for jags
JAGS.seed<-ceiling(runif(1,1,10000000))

# to test the function
test=F
if(test==T){
  varname <- 4
  sitename <- 6
  scale <- 2
}

# key for which response variable corresponds to which index
if(varname == 1){ # WUE is either estimated using GPP or SIF
  varname = "WUE_GPP"
} else if(varname == 2){
  varname = "WUE_SIF"
}
# key for which response site corresponds to which index, will load a dataframe called "dat"
if(sitename == 1){
  sitename = "lae"
  load("./data_clean/laedat.RData")
} else if(sitename == 2){
  sitename = "crk"
  load("./data_clean/crkdat.RData")
} else if(sitename == 3){
  sitename = "geb"
  load("./data_clean/gebdat.RData")
} else if(sitename == 4){
  sitename = "lnf"
  load("./data_clean/lnfdat.RData")
} else if(sitename == 5){
  sitename = "yat"
  load("./data_clean/yatdat.RData")
} else if(sitename == 6){
  sitename = "lm1g"
  load("data_clean/lm1gdat.RData")
}

# Load libraries
# if(!"dplyr" %in% installed.packages()) {
#   install.packages("dplyr", repos = 'https://mirror.las.iastate.edu/CRAN/') # pick appropriate CRAN mirror
# }
# if(!"tibble" %in% installed.packages()) {
#   install.packages("tibble", repos = 'https://mirror.las.iastate.edu/CRAN/') # pick appropriate CRAN mirror
# }
# if(!"jagsUI" %in% installed.packages()) {
#   install.packages("jagsUI", repos = 'https://mirror.las.iastate.edu/CRAN/') # pick appropriate CRAN mirror
# }
# if(!"mcmcplots" %in% installed.packages()) {
#   install.packages("mcmcplots", repos = 'https://mirror.las.iastate.edu/CRAN/') # pick appropriate CRAN mirror
# }
# if(!"gsubfn" %in% installed.packages()) {
#   install.packages("gsubfn", repos = 'https://mirror.las.iastate.edu/CRAN/') # pick appropriate CRAN mirror
# }
library(tidyverse) # for arranging dataframes
library(tibble) # rowid_to_column
library(jagsUI) # for running the jags model
library(mcmcplots) # convergence plots
library(gsubfn) # for gsub for table org
library(bigleaf) # for conversions

# Load self-made functions
source("./scripts/functions.R")
source("./scripts/run_modDEPARTSAM_function.R")

# Aggregate according to scale
# if(scale == "hour"){ # aggregate by hour
#   dat <- dat %>%
#     select(-X) %>%
#     mutate(hour = format(dat$TIMESTAMP, format ="%H"), DOY = floor(DOY)) %>%
#     group_by(DOY,hour) %>%
#     summarise(across(everything(), mean, na.rm = T)) %>%
#     ungroup()
# }
# if(scale == "daily"){ # aggregate by day
#   dat <- dat %>%
#     select(-X) %>%
#     mutate(DOY = floor(DOY), TIMESTAMP = as.Date(TIMESTAMP, format= "%Y-%m-%d")) %>%
#     group_by(DOY)  %>%
#     summarise(across(everything(), mean, na.rm = T))
# }

# Load site metadata for variables needed for process-based evaporation model
df_soil <- read.csv("./data_misc/soildata.csv")

df_siteinfo <- df_soil %>% # filter for the site we're running the model for
  filter(sitecode == sitename)

Z = df_siteinfo$ref_height_m # reference height (m)
h = 	df_siteinfo$elev_m # elevation
fc = df_siteinfo$FC # field capacity
fclay = df_siteinfo$clay/100 # clay fraction
fsand = df_siteinfo$sand/100 # sand fraction
dat$SWC_shall <- dat$SWC_shall
dat$SWC_deep <- dat$SWC_deep

# Calculate intermediate parameters needed for evaporation model
# Need dataframe with columns for PA, TA, TS, WS, RH, SWC_shall
# Need separate data for Z, fsand, fclay, fc
conv.fact.time <- 3600 # we need to convert seconds to hours
#conv.fact.time <- ifelse(scale=="day", 86400, conv.fact.time) # if we need to convert seconds to days
dat <- get_evap(dat, Z = Z, h = h, fc = fc, fclay = fclay, fsand = fsand, conv.fact.time = conv.fact.time)

lowdev = F
newinits = F


# Run model 2
run_modDEPARTSAM(dat, varname = varname, sitename = sitename, # data and naming conventions
         Z = Z, h = h, fc = fc, fclay = fclay, fsand = fsand, # soil model data
         newinits = newinits, overwrite = T, post_only = F, lowdev=lowdev) # model specifics