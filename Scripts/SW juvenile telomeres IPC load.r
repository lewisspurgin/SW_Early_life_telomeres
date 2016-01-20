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

# Load data ---------------------------------------------------------------
dd0 <- read.csv('Data/TL2 output 15 Dec 2015.csv') #main dataset
terr <- read.csv('Data/SW TL all territories.csv')
status <- read.csv('Data/SW TL all breedstatus.csv')
insects <- read.csv('Data/Insects.csv')
