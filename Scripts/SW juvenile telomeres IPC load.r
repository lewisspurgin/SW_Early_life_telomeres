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
library(flexsurv)

# Load data ---------------------------------------------------------------
dd <- read.csv('Data/SW TL main data for analysis 2.csv') #main dataset
terr <- read.csv('Data/SW TL all territories.csv')
status <- read.csv('Data/SW TL all breedstatus.csv')
insects <- read.csv('Data/Insects.csv')
