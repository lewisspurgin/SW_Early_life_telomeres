#################################################################################*
## Early life telomere dynamics and late-life fitness in a wild bird population
## LOAD THE DATA AND PACKAGES
#################################################################################*


# Clear out R -------------------------------------------------------------
rm(list=ls())



library(MuMIn)
library(plyr)
library(ggplot2)
library(nlme)
library(magrittr)
library(rptR)
library(Rmisc)
library(lme4)
library(arm)
library(car)

# Load data ---------------------------------------------------------------
dd0 <- read.csv('Data/REL_TL_EB_EAF_24 May aft TL queries.csv') #main dataset
terr <- read.csv('Data/SW TL all territories.csv')
status <- read.csv('Data/SW TL all breedstatus.csv')
insects <- read.csv('Data/Insects.csv')
dens <- read.table('Data/PsizeFP.txt')
chickinfo <- read.csv('Data/Table chick info.csv')
mal <- read.csv('Data/Malaria.csv')

