#################################################################################*
## Early life telomere dynamics and late-life fitness in a wild bird population
## LOAD THE DATA AND PACKAGES
#################################################################################*


# Clear out R -------------------------------------------------------------
rm(list=ls())


# Load relevant libraries -------------------------------------------------
library(MuMIn)
library(arm)
library(plyr)
library(ggplot2)
library(gvlma)
library(survival)
library(nlme)
library(magrittr)
library(car)

# Load data ---------------------------------------------------------------
dd0 <- read.csv('Data/REL_TL_EB_EAF_25 Feb FULL DATA.csv') #main dataset
terr <- read.csv('Data/SW TL all territories.csv')
status <- read.csv('Data/SW TL all breedstatus.csv')
insects <- read.csv('Data/Insects.csv')
dens <- read.table('Data/PsizeFP.txt')
chickinfo <- read.csv('Data/Table chick info.csv')
