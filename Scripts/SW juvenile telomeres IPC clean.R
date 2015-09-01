#################################################################################*
## Early life telomere dynamics and late-life fitness in a wild bird population
## CLEAN THE DATA
#################################################################################*


#Average repeats of blood samples
av <- ave(dd$TL,c(dd$BloodID,dd$Status))
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


# Catch Year, Catch date and death year ----------------------------------


dd$CatchYear <- as.numeric(substr(dd$CatchDate,7,10))
dd$DeathYear <- as.numeric(substr(dd$DeathDate,7,10))
dd$CatchDate <- as.Date(dd$CatchDate,"%d/%m/%Y")




# Age data ----------------------------------------------------------------

dd$Age <- dd$CatchYear-dd$LayYear


#Sort out and order ageclass levels
dd <- droplevels(dd[dd$Ageclass != '',])
levels(dd$Ageclass) <- c('A','CH','FL','FL','FL','SA')
dd$Ageclass <- factor(dd$Ageclass,levels = c('CH','FL','FL','FL','SA','A'))



dd$Agemonths <- ifelse(dd$Ageclass == 'CH',1,
                       ifelse(dd$Ageclass == 'FL',3,
                              ifelse(dd$Ageclass == 'SA',9,dd$Age*12)))
dd$AgemonthF <- ifelse(dd$Ageclass == 'CH','<1',
                       ifelse(dd$Ageclass == 'FL','1-9',
                              ifelse(dd$Ageclass == 'SA','9-12','>12')))
dd$AgemonthF <- factor(dd$AgemonthF,levels = c('<1','1-9','9-12','>12'))

dd <- subset(dd,Agemonths>0)

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


dd <- droplevels(subset(subset(dd,TL>1000),TL<12000))
#dd <- subset(dd,BodyMass>11)
dd <- subset(dd,Insect<8)


# Subset juveniles --------------------------------------------------------------


juv <- droplevels(subset(dd,Ageclass %in% c('CH','FL','OFL','SA')))
adults <- droplevels(subset(dd,Ageclass == 'A'))



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

xf <- subset(juv,Status == 'XF')

juv <- subset(juv,!is.na(Tarsus))
juv <- subset(juv,!is.na(BodyMass))


FlSA <- subset(juv,Ageclass!='CH')
chicks <- subset(juv,Ageclass == 'CH')
chicks <- subset(chicks,Status != 'XF')

#FlSA$LayYear <- factor(FlSA$LayYear)
FlSAall <- droplevels(FlSA[complete.cases(FlSA),])
chickall <- droplevels(chicks[complete.cases(chicks),])




# Look at telomere loss ---------------------------------------------------

x1 <- aggregate(CatchDate~BirdID,adults,min)
x2 <- merge(adults,x1)
x2 <- x2[!(duplicated(x2$BirdID)),]


x3 <- aggregate(TLKB~BirdID,juv,max)
x4 <- merge(x3,juv)
x4 <- x4[!(duplicated(x4$BirdID)),]

x2 <- x2[x2$BirdID %in% x4$BirdID,]
x4 <- x4[x4$BirdID %in% x2$BirdID,]

x3 <- x3[order(x3$BirdID),]
x4 <- x4[order(x4$BirdID),]

xx1 <- x4$TLKB-mean(x4$TLKB)
xx2 <- x2$TLKB-mean(x2$TLKB)          
         
rho <- cor(x4$TLKB,x2$TLKB)

(rho*xx1)-xx2

Loss <- data.frame(x4,
                   Loss = x4$TL-x2$TL,
                   TimeDiff = as.numeric(x2$CatchDate-x4$CatchDate),
                   D=(rho*xx1)-xx2)
Loss$TROC <- with(Loss,D/TimeDiff)

Loss <- subset(subset(Loss,TROC > -0.005),TROC < 0.01)



# Get rid of stuff not to be used -----------------------------------------

rm(status,helpers,hatchdate,x1,x2,x3,x4)