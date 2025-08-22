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


# Set defined R seed
set.seed(SEED, kind = NULL, normal.kind = NULL)
# Generate "random" seed for jags
JAGS.seed<-ceiling(runif(1,1,10000000))

# to test the function, switch test=T
test=F
if(test==T){
  varname <- 4
  sitename <- 1
  scale = 1
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
source("./scripts/run_modSAM_function.R")

# key for which response variable corresponds to which index
if(varname == 1){
  varname = "GPP"
} else if(varname == 2){
  varname = "SIF_O2A"
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
  sitename = "lm1g"
  load("data_clean/maggdat.RData")
}


lowdev = F

# Run mod1
run_modSAM(dat, varname = varname, sitename = sitename,  # data and naming conventions
         newinits = F, overwrite = T, lowdev=lowdev) # model specifics



