#################################################################################*
## Early life telomere dynamics and late-life fitness in a wild bird population
## CLEAN THE DATA
#################################################################################*

av <- ave(dd$TL,dd$BloodID)
dd$TL <- av

dd <- unique(dd)

# Telomere data -----------------------------------------------------------


#Add in some useful new TL variables

dd$LogTL <- log(dd$TL)
dd$TLKB <- dd$TL/1000


# Weird variable names ----------------------------------------------------


#Change names of weird database variables
colnames(dd)[colnames(dd) == 'OccasionDate'] <- 'CatchDate'
colnames(dd)[colnames(dd) == 'LastSighting'] <- 'DeathDate'



# Catch Year, Catch month and death year ----------------------------------


#Extract catch year, death year using substring
dd$CatchYear <- as.numeric(substr(dd$CatchDate,7,10))
dd$DeathYear <- as.numeric(substr(dd$DeathDate,7,10))





# Age data ----------------------------------------------------------------

dd$Age <- dd$CatchYear-dd$LayYear
dd$Age[dd$Ageclass]

#Sort out and order ageclass levels
dd <- droplevels(dd[dd$Ageclass != '',])
levels(dd$Ageclass) <- c('A','CH','FL','FL','FL','SA')
dd$Ageclass <- factor(dd$Ageclass,levels = c('CH','FL','FL','FL','SA','A'))

dd$Age[dd$Ageclass != 'A'] <- 0

# Survival and lifespan ---------------------------------------------------


#Add some survival variables
dd$RemainingLife <- dd$DeathYear-dd$CatchYear
dd$SurvivedNext <- ifelse(dd$RemainingLife>0,1,0)


dd$Lifespan <- dd$DeathYear-dd$LayYear

dd$Died <- ifelse(dd$DeathYear<2013,1,0)


# Sex ---------------------------------------------------------------------


#Add in an easy to understand sex variable
dd$Sex <- ifelse(dd$SexEstimate == 1,'Males','Females')








# Tarsus ------------------------------------------------------------------


#Combine Tarsus data into one column

dd$Tarsus <- NA
for(i in 1:nrow(dd))
{
  ifelse(is.na(dd$RightTarsus[i]),
         ifelse(is.na(dd$LeftTarsus[i]),
                dd$Tarsus[i] <- NA,
                dd$Tarsus[i] <- dd$LeftTarsus[i]),
         dd$Tarsus[i] <- dd$RightTarsus[i])
}

dd$RightTarsus <- NULL
dd$LeftTarsus <- NULL





# TQ and insects ----------------------------------------------------------


terr <- terr[complete.cases(terr),] #Get rid of blank rows
#terr <- subset(terr,SummerIndex>0.9) #Get rid of winter seasons

insects$Year <- as.numeric(substr(insects$SamplingDate,7,10))

# taking average for year
yearmean <- tapply(terr$TQcorrected,terr$Year,median)
terrmean <- tapply(terr$TQcorrected,terr$TerritoryID,median)

Iyearmean <- tapply(insects$InsectCount,insects$Year,mean)
Iyearsd <- tapply(insects$InsectCount,insects$Year,sd)

dd$TQspace <- NA
dd$TQtime <- NA
dd$Insect <- NA

for(i in 1:nrow(dd))
{
  if(dd$LayYear[i] %in% names(yearmean))
  {
    dd$TQtime[i] <- yearmean[names(yearmean) == dd$LayYear[i]]
  }
  
  if(dd$TerritoryID[i] %in% names(terrmean))
  {
    dd$TQspace[i] <- terrmean[names(terrmean) == dd$TerritoryID[i]]
  }
  
  if(dd$LayYear[i] %in% names(Iyearmean))
  {
    dd$Insect[i] <- Iyearmean[names(Iyearmean) == dd$LayYear[i]]  
  }
}

dd$InsectF <- ifelse(dd$Insect > 14, 'High','Low')
dd$TQspace <- log(dd$TQspace)
dd$TQtime <- log(dd$TQtime)


# Remove unwanted data ----------------------------------------------------


dd <- droplevels(subset(subset(dd,TL>1000),TL<15000))
dd <- subset(dd,BodyMass>11)


#Telomere factor variable -----------------------------------------


mymed <- median(dd$TLKB,na.rm=T)
dd$TLF <- ifelse(dd$TLKB > mymed,'Long telomeres','Short telomeres')

# Subset juveniles --------------------------------------------------------------

#dd$LayYear <- factor(dd$LayYear)
juv <- droplevels(subset(dd,Ageclass %in% c('CH','FL','SA')))



# Helpers and social group size -------------------------------------------

status <- subset(status,BreedGroupID %in% juv$BreedGroupID)

for(i in 1:nrow(juv))
{
  currentBG <- juv$BreedGroupID[i]
  currentdata <- subset(subset(status,BreedGroupID == currentBG),BirdID !=juv$BirdID[i])
  juv$GroupSize[i] <- nrow(currentdata)
  juv$Helper[i] <- nrow(subset(currentdata,Status == 'H'))
  juv$OtherJuvs[i] <- nrow(subset(currentdata,Status %in% c('CH','FL','OFL')))
  juv$NonHelper[i] <- nrow(subset(currentdata,Status %in% c('AB','ABX')))
}

#Get rid of weird group sizes
juv$HelperC <- juv$Helper-juv$NonHelper
juv$HelperF <- ifelse(juv$Helper==1,paste(juv$Helper, 'Helper'),paste(juv$Helper, 'Helpers'))


# Subset chicks and other juvs -----------------------------------------------------------


FlSA <- subset(juv,Ageclass!='CH')
FlSA <- subset(FlSA,!is.na(Tarsus))
FlSA <- subset(FlSA,!is.na(BodyMass))
FlSA$Condition <- summary(lm(BodyMass~Tarsus+Sex,data=FlSA))$resid
FlSAall <- droplevels(FlSA[complete.cases(FlSA),])


# Get rid of stuff not to be used -----------------------------------------

rm(status,helpers,hatchdate)



# Longitudinal telomere loss in adults ------------------------------------

adults <- subset(dd,Ageclass=='A')

FlSA$TLloss <- NA

for(i in 1:nrow(FlSA))
{
  currentjuv <- FlSA$BirdID[i]
  if(currentjuv %in% adults$BirdID)
  {
    currentdata <- subset(adults,BirdID==currentjuv)
    youngest <- subset(currentdata,Age==min(currentdata$Age))
    FlSA$TLloss[i] <- FlSA$TLKB[i]-youngest$TLKB[1]
  }

}


