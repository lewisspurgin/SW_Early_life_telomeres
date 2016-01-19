#################################################################################*
## Early life telomere dynamics and late-life fitness in a wild bird population
## CLEAN THE DATA
#################################################################################
#Average repeats of blood samples
av <- ave(dd$TL,c(dd$BloodID,dd$Status,dd$PlateID))
dd$TL <- av
dd <- dd[!(duplicated(dd$BloodID)),]


# Weird variable names ----------------------------------------------------

colnames(dd)[colnames(dd) == 'OccasionDate'] <- 'CatchDate'
colnames(dd)[colnames(dd) == 'MaxOfMaxOfSeenDate'] <- 'DeathDate'
colnames(dd)[colnames(dd) == 'MinOfFieldPeriodID'] <- 'FieldPeriodID'


# Catch Year, Catch date and death year ----------------------------------


dd$CatchYear <- as.numeric(substr(dd$CatchDate,7,10))
dd$DeathYear <- as.numeric(substr(dd$DeathDate,7,10))
dd$CatchDate <- as.Date(dd$CatchDate,"%d/%m/%Y")

dd$Season <- ifelse(as.numeric(format(dd$CatchDate,'%m')) %in% c(4:10),
                    'Summer','Winter')


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



# Remove unwanted data/outliers ----------------------------------------------------


dd <- droplevels(subset(subset(dd,TL>1000),TL<30000))
dd <- subset(dd,BodyMass>5)
dd <- subset(dd,Tarsus>17)


# Telomere data -----------------------------------------------------------

dd$TLKB <- dd$TL/1000
dd$LogTL <- log10(dd$TL)




# Age data ----------------------------------------------------------------

dd$Age <- dd$CatchYear-dd$LayYear


#Sort out and order ageclass levels
dd <- droplevels(dd[dd$Ageclass != '',])
levels(dd$Ageclass) <- c('A','CH','FL','FL','FL','SA')
dd$Ageclass <- factor(dd$Ageclass,levels = c('CH','FL','FL','FL','SA','A'))
dd$Fledged <- ifelse(dd$Ageclass == 'CH','Nestlings',
                     ifelse(dd$Ageclass == 'A','Adults',
                            'Fledglings'))
dd$Fledged <- factor(dd$Fledged,levels = c('Nestlings','Fledglings','Adults'))


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







# Subset juveniles --------------------------------------------------------------


juv <- droplevels(subset(dd,Ageclass %in% c('CH','FL','OFL','SA')))
adults <- droplevels(subset(dd,Ageclass == 'A'))

juv$cenTL <- NA

#Centre telomere length by birth year
for(i in 1:nrow(juv))
{
  currentdata <- subset(juv,FieldPeriodID == FieldPeriodID[i])
  juv$cenTL[i] <- (juv$TLKB[i] - mean(currentdata$TLKB))/sd(currentdata$TLKB)
}

mymed <- mean(juv$TLKB,na.rm=T)
juv$TLF <- ifelse(juv$TLKB > mymed,'Long telomeres','Short telomeres')


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

juv <- subset(juv,GroupSize>0)





# TQ and insects ----------------------------------------------------------


terr <- terr[complete.cases(terr),] #Get rid of blank rows


insects <- subset(insects,FieldPeriodID != 26)
insects$Insectcen <- (insects$MeanInsects-mean(insects$MeanInsects))/sd(insects$MeanInsects)

# take average for year
terrmean <- tapply(terr$TQcorrected,terr$TerritoryID,mean)

juv$TQ <- NA
juv$TQI <- NA
juv$cenTQ <- NA
juv$Insect <- NA

for(i in 1:nrow(juv))
{
  if(juv$TerritoryID[i] %in% names(terrmean))
  {
    juv$TQ[i] <- terrmean[names(terrmean) == juv$TerritoryID[i]]
    juv$TQI[i] <- juv$TQ[i]/juv$GroupSize[i]
  }
  
  if(juv$FieldPeriodID[i] %in% insects$FieldPeriodID)
  {
    juv$Insect[i] <- insects$Insectcen[insects$FieldPeriodID == juv$FieldPeriodID[i]] 
  } 
}





# Subset Fledglings and subadults -----------------------------------------------------------

xf <- subset(juv,Status == 'XF')

juv <- subset(juv,!is.na(Tarsus))
juv <- subset(juv,!is.na(BodyMass))
juv <- subset(juv,Status!='XF')


for(i in 1:nrow(juv))
{
  currentdata <- subset(juv,FieldPeriodID == juv$FieldPeriodID[i])
  juv$cenTQ[i] <- (juv$TQ[i] - mean(currentdata$TQ))/sd(currentdata$TQ)
  juv$cenTarsus[i] <- (juv$Tarsus[i] -  mean(currentdata$Tarsus))/sd(currentdata$Tarsus)
  juv$cenHelper[i] <- (juv$Helper[i] -  mean(currentdata$Helper))/sd(currentdata$Helper)
  juv$cenNonHelper[i] <- (juv$NonHelper[i] -  mean(currentdata$NonHelper))/sd(currentdata$NonHelper)
}

juv <- subset(juv,cenTQ<4)


FlSA <- subset(juv,Ageclass!='CH')
chicks <- subset(juv,Ageclass == 'CH')

juvall <- droplevels(juv[complete.cases(juv),])
FlSAall <- droplevels(FlSA[complete.cases(FlSA),])
chickall <- droplevels(chicks[complete.cases(chicks),])



# Look at telomere loss ---------------------------------------------------

x1 <- aggregate(CatchDate~BirdID,adults,min)
x2 <- merge(adults,x1)
x2 <- x2[!(duplicated(x2$BirdID)),]


x3 <- aggregate(CatchDate~BirdID,juv,min)
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
for(i in 1:nrow(Loss))
{
  currentdata <- Loss[Loss$FieldPeriodID == Loss$FieldPeriodID[i],]
  Loss$cenTROC[i] <- (Loss$TROC[i]-mean(currentdata$TROC)/sd(currentdata$TROC))
}

# Get rid of stuff not to be used -----------------------------------------

rm(status,helpers,hatchdate,x1,x2,x3,x4)




# Field period average data -----------------------------------------------

juvseason <- ddply(juv,
                   .(FieldPeriodID,Season),
                   summarize,
                   TLKBmean = mean(LogTL),
                   TLKBse = se(LogTL),
                   Tarsus = mean(Tarsus),
                   Insect = mean(Insect),
                   Lifespan = mean(RemainingLife),
                   Lifespanse = se(RemainingLife),
                   CatchYear = mean(CatchYear),
                   n = length(TLKB))
juvseason$cenTL <- juvseason$TLKBmean-mean(juvseason$TLKBmean)
juvseason <- subset(juvseason,n>5)

chickseason <- ddply(chicks,
                   .(FieldPeriodID,Season),
                   summarize,
                   TLKBmean = mean(TLKB),
                   TLKBse = se(TLKB),
                   Tarsus = mean(Tarsus),
                   Insect = mean(Insect),
                   Lifespan = mean(RemainingLife),
                   Lifespanse = se(RemainingLife),
                   CatchYear = mean(CatchYear),
                   n = length(TLKB))
chickseason$cenTL <- chickseason$TLKBmean-mean(chickseason$TLKBmean)
chickseason <- subset(chickseason,n>5)


flseason <- ddply(FlSA,
                     .(FieldPeriodID,Season),
                     summarize,
                     TLKBmean = mean(TLKB),
                     TLKBse = se(TLKB),
                     Tarsus = mean(Tarsus),
                     Insect = mean(Insect),
                     Lifespan = mean(RemainingLife),
                     Lifespanse = se(RemainingLife),
                     CatchYear = mean(CatchYear),
                     n = length(TLKB))
flseason$cenTL <- flseason$TLKBmean-mean(flseason$TLKBmean)
flseason <- subset(flseason,n>5)

loss.season <- ddply(Loss,
                   .(FieldPeriodID,Season),
                   summarize,
                   TROCmean = mean(TROC),
                   TROCse = se(TROC),
                   Insect = mean(Insect),
                   Lifespan = mean(RemainingLife),
                   Lifespanse = se(RemainingLife),
                   n = length(TROC),
                   CatchYear = mean(CatchYear))
loss.season <- subset(loss.season,TROCse<0.01)
