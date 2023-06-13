### Script to filter data
### input data: formatted data from data_input
### output data: csv dataframe for each site with SIF filtered for quality

library(tidyverse)
library(gridExtra)
source("./scripts/functions.R")

#################################################### Laegeren
# Load Laegeren data
laedatin <- read.csv("./data_formatted/CH-Lae/CH-Lae_dat.csv")

#################################################### Crk
# Read in Crk data
crkdatin <- read.csv("./data_formatted/Crk/Crk_dat.csv")

#################################################### Gebesee
# Load Gebesee data
gebdatin <- read.csv("./data_formatted/DE-Geb/DE-Geb_dat.csv")

#################################################### Leinefelde
# Load Leinefelde data
lnfdatin <- read.csv("./data_formatted/DE-Lnf/DE-Lnf_dat.csv")

 lnfdat<- lnfdatin#[5000:6000,] #%>%
#   filter(!is.na(SIF_O2B))

p1 <- ggplot(data = lnfdat) +
  geom_line(aes(x=DOY, y=SIF_O2B))+
  #xlim(5000,6000)+
  theme_bw()
p1

p2 <- ggplot(data = lnfdat) +
  geom_line(aes(x=DOY, y=Ref_750))+
  theme_bw()

p3 <- ggplot(data = lnfdat) +
  geom_line(aes(x=DOY, y=NDVI))+
  theme_bw()

grid.arrange(p1,p2,p3, nrow=3)

lnfdat$SIF_O2B <- ifelse(lnfdat$Ref_750 < mean(lnfdat$Ref_750), NA, lnfdat$SIF_O2B)

test <- lnfdatin$SIF_O2B
library(zoo)
lnfdat$test <- na.approx(lnfdat$SIF_O2B, na.rm = FALSE)


lnfdatin$SIF_O2B <- gap.fill(lnfdatin$SIF_O2B, max.gap = 12, decimals = 5, df = FALSE)


lnfdat<- lnfdatin %>%
  filter(!is.na(SIF_O2B))

ggplot(data = lnfdat) +
  geom_point(aes(x=X, y=SIF_O2B))+ 
  geom_point(aes(x=X, y=Ref_750), color = "red")+ 
  theme_bw()


#################################################### Yat
# Read in Yatir data
yatdatin <- read.csv("./data_formatted/IL-Yat/IL-Yat_dat.csv")


##################################################### Jrs

# Load JRS data

##################################################### Sq

# Load sq data




