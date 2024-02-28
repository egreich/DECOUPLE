#### Run model 1

# Set run params
args<-commandArgs(TRUE)
print(args)
print("varname:")
(varname <- as.numeric(args[1]))
print("sitename:")
(sitename <- as.numeric(args[2]))
print("seed:")
(SEED <- as.numeric(args[3]))
print("scale:")
(scale <- as.numeric(args[4]))

# Set defined R seed
set.seed(SEED, kind = NULL, normal.kind = NULL)
# Generate "random" seed for jags
JAGS.seed<-ceiling(runif(1,1,10000000))

# to test the function, switch test=T
test=F
if(test==T){
  varname <- 1
  sitename <- 1
  scale = 2
}

# key for which response variable corresponds to which index
if(varname == 1){
  varname = "ET"
} else if(varname == 2){
  varname = "GPP"
} else if(varname == 3){
  varname = "SIF_O2A"
} else if(varname == 4){
  varname = "WUE_GPP"
} else if(varname == 5){
  varname = "WUE_SIF"
} else if(varname == 6){
  varname = "T_TEA"
} else if(varname == 7){
  varname = "WUE_GPP_TEA"
} else if(varname == 8){
  varname = "WUE_SIF_TEA"
}
# key for which response site corresponds to which index
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
  sitename = "jrs"
  load("./data_clean/jrsdat.RData")
} else if(sitename == 7){
  sitename = "sq"
  load("./data_clean/sqdat.RData")
}

if(scale == 1){
  scale = "halfhour"
} else if(scale == 2){
  scale = "hour"
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

# Load self-made functions
source("./scripts/functions.R")
source("./scripts/run_mod1_function.R")

# Aggegrate according to scale
if(scale == "hour"){ # aggregate by hour
  dat <- dat %>%
    select(-X) %>%
    mutate(hour = format(dat$TIMESTAMP, format ="%H"), DOY = floor(DOY)) %>%
    group_by(DOY,hour) %>%
    summarise(across(everything(), mean, na.rm = T))
}
# if(scale == "daily"){ # aggregate by day
#   dat <- dat %>%
#     select(-X) %>%
#     mutate(DOY = floor(DOY), TIMESTAMP = as.Date(TIMESTAMP, format= "%Y-%m-%d")) %>%
#     group_by(DOY) %>%
#     summarise(across(everything(), mean, na.rm = T))
# }

lowdev = F

# Run mod1
run_mod1(dat, varname = varname, sitename = sitename, scale,  # data and naming conventions
         newinits = T, overwrite = T, lowdev=lowdev) # model specifics



