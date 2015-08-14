#################################################################################*
## Early life telomere dynamics and late-life fitness in a wild bird population
## CLEAN THE DATA
#################################################################################*

#Average repeats of blood samples
av <- ave(dd$TL,dd$BloodID)
dd$TL <- av

dd <- unique(dd)




# Telomere data -----------------------------------------------------------

dd$LogTL <- log(dd$TL)
dd$TLKB <- dd$TL/1000
mymed <- median(dd$TLKB,na.rm=T)
dd$TLF <- ifelse(dd$TLKB > mymed,'Long telomeres','Short telomeres')


# Weird variable names ----------------------------------------------------

colnames(dd)[colnames(dd) == 'OccasionDate'] <- 'CatchDate'
colnames(dd)[colnames(dd) == 'MaxOfMaxOfSeenDate'] <- 'DeathDate'
colnames(dd)[colnames(dd) == 'MinOfFieldPeriodID'] <- 'FieldPeriodID'


# Catch Year, Catch month and death year ----------------------------------

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

dd$Agemonths <- ifelse(dd$Ageclass == 'CH',1,
                       ifelse(dd$Ageclass == 'FL',6,
                              ifelse(dd$Ageclass == 'SA',9,dd$Age*12)))

#Delta age
for(i in 1:nrow(dd))
{
  currentdata = subset(dd,BirdID == dd$BirdID[i])
  dd$MeanAge[i] <- mean(currentdata$Age)
  dd$DeltaAge[i] <- dd$Age[i]-dd$MeanAge[i]
}



# Survival and lifespan ---------------------------------------------------

dd$RemainingLife <- dd$DeathYear-dd$CatchYear
dd$SurvivedNext <- ifelse(dd$RemainingLife>0,1,0)
dd$Lifespan <- dd$DeathYear-dd$LayYear
dd$Died <- ifelse(dd$DeathYear<2013,1,0)





# Sex ---------------------------------------------------------------------

dd$Sex <- ifelse(dd$SexEstimate == 1,'Males','Females')




# Tarsus ------------------------------------------------------------------

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
terr <- subset(terr,SummerIndex>0.9) #Get rid of winter seasons

insects <- subset(insects,FieldPeriodID != 26)

# take average for year
yearmean <- tapply(terr$TQcorrected,terr$Year,mean)
terrmean <- tapply(terr$TQcorrected,terr$TerritoryID,mean)

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
  
  if(dd$FieldPeriodID[i] %in% insects$FieldPeriodID)
  {
    dd$Insect[i] <- insects$MeanInsects[insects$FieldPeriodID == dd$FieldPeriodID[i]] 
  } 
}

dd$InsectF <- ifelse(dd$Insect > 14, 'High','Low')
dd$TQspace <- log(dd$TQspace)
dd$TQtime <- log(dd$TQtime)


# Remove unwanted data ----------------------------------------------------


dd <- droplevels(subset(subset(dd,TL>1000),TL<15000))
dd <- subset(dd,BodyMass>11)
dd <- subset(dd,Insect<8)


# Subset juveniles --------------------------------------------------------------


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



# Subset Fledglings and subadults -----------------------------------------------------------


juv <- subset(juv,!is.na(Tarsus))
juv <- subset(juv,!is.na(BodyMass))


FlSA <- subset(juv,Ageclass!='CH')
chicks <- subset(juv,Ageclass == 'CH')

#FlSA$LayYear <- factor(FlSA$LayYear)
FlSAall <- droplevels(FlSA[complete.cases(FlSA),])
chickall <- droplevels(chicks[complete.cases(chicks),])


# Get rid of stuff not to be used -----------------------------------------

rm(status,helpers,hatchdate)





