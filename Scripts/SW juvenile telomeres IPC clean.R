#################################################################################*
## Early life telomere dynamics and late-life fitness in a wild bird population
## CLEAN THE DATA
#################################################################################

#Only use Ellie's samples
dd0 <- subset(dd0,Whodunnit == 'EAF')

#Average repeats of blood samples
av <- ave(dd0$RTL,c(dd0$BloodID,dd0$Status,dd0$PlateID))
dd0$RTL2 <- av
dd <- dd0[!(duplicated(dd0$BloodID)),]
dd$RTL <- dd$RTL2
dd$RTL2 <- NULL

# Weird variable names ----------------------------------------------------

colnames(dd)[colnames(dd) == 'OccasionDate'] <- 'CatchDate'
colnames(dd)[colnames(dd) == 'MaxOfMaxOfSeenDate'] <- 'DeathDate'
colnames(dd)[colnames(dd) == 'MinOfFieldPeriodID'] <- 'FieldPeriodID'


# Catch Year, Catch date and death year ----------------------------------


dd$CatchYear <- as.numeric(substr(dd$CatchDate,7,10))
dd$DeathYear <- as.numeric(substr(dd$DeathDate,7,10))
dd$CatchDate <- as.Date(dd$CatchDate,"%d/%m/%Y")
dd$DeathDate <- as.Date(dd$DeathDate,"%d/%m/%Y")



dd$Season <- ifelse(as.numeric(format(dd$CatchDate,'%m')) %in% c(4:10),
                    'Major','Minor')
dd <- subset(dd,Season == 'Major')


# Age data ----------------------------------------------------------------

dd$Age <- (dd$CatchYear-dd$LayYear)+1

#Sort out and order ageclass levels
dd <- droplevels(dd[dd$Ageclass != '',])
levels(dd$Ageclass) <- c('A','CH','FL','FL','FL','SA')
dd$Ageclass <- factor(dd$Ageclass,levels = c('CH','FL','FL','FL','SA','A'))
dd$Fledged <- ifelse(dd$Ageclass == 'CH','Nestlings',
                     ifelse(dd$Ageclass == 'A','Adults',
                            'Fledglings'))
dd$Fledged <- factor(dd$Fledged,levels = c('Nestlings','Fledglings','Adults'))


dd$Agemonths <- ifelse(dd$Ageclass == 'CH',1,
                       ifelse(dd$Ageclass == 'FL',6,
                              ifelse(dd$Ageclass == 'SA',10,dd$Age*12)))

dd <- subset(dd,Agemonths!=12)

dd$Age[dd$Fledged != 'Adults'] <- 1




# Tarsus and delta age ------------------------------------------------------------------

dd$Tarsus <- NA
dd$DeltaAge <- NA
dd$Malaria <- NA

malaria <- subset(malaria,!duplicated(BloodID))
malaria <- subset(malaria,Consensus <2)

for(i in 1:nrow(dd))
{
  ifelse(is.na(dd$RightTarsus[i]),
         ifelse(is.na(dd$LeftTarsus[i]),
                dd$Tarsus[i] <- NA,
                dd$Tarsus[i] <- dd$LeftTarsus[i]),
         dd$Tarsus[i] <- dd$RightTarsus[i])
  
  dd$DeltaAge[i] <- dd$Agemonths[i] - mean(subset(dd,BirdID == dd$BirdID[i])$Agemonths)
  if(dd$BloodID[i] %in% malaria$BloodID)
     {
       dd$Malaria[i] <- subset(malaria,BloodID == dd$BloodID[i])$Consensus
     }
}

dd$RightTarsus <- NULL
dd$LeftTarsus <- NULL



# Remove unwanted data/outliers ----------------------------------------------------

dd <- subset(dd,RTL<2)
dd <- subset(dd,RTL > 0.04)
dd <- subset(dd,BodyMass>5)
dd <- subset(dd,Tarsus>17)




# Survival and lifespan ---------------------------------------------------

dd$RemainingLife <- round(as.numeric(dd$DeathDate-dd$CatchDate)/365,0)
dd$SurvivedNext <- ifelse(dd$RemainingLife<1,0,1)

dd$SurvivedNext <- ifelse(dd$RemainingLife>1,1,0)
dd$Lifespan <- (dd$DeathYear-dd$LayYear)+1
dd$Died <- ifelse(dd$DeathYear<2013,1,0)



# Exclude weird seasons and early years ----------------------------------

dd <- subset(dd,LayYear>1997)
#dd <- subset(dd,Season == 'Major')
dd <- subset(dd,FieldPeriodID != 27)

# Sex ---------------------------------------------------------------------

dd$Sex <- ifelse(dd$SexEstimate == 1,'Males','Females')







# Subset juveniles --------------------------------------------------------------


juv <- droplevels(subset(dd,Ageclass %in% c('CH','FL','OFL','SA')))
adults <- droplevels(subset(dd,Ageclass == 'A'))



#Centre telomere length by birth year
juv$cenTL <- NA
juv$cohortTL <- NA
for(i in 1:nrow(juv))
{
  currentdata <- subset(juv,FieldPeriodID == FieldPeriodID[i])
  if(nrow(currentdata>4))
  {
    juv$cenTL[i] <- (juv$RTL[i] - mean(currentdata$RTL))/sd(currentdata$RTL)
    juv$cohortTL[i] <- median(currentdata$RTL)  
  }

}

mymed <- median(juv$cenTL,na.rm=T)
juv$cenTLF <- ifelse(juv$cenTL < mymed,'Short telomeres','Long telomeres')

mymed <- median(juv$cohortTL,na.rm=T)
juv$coTLF <- ifelse(juv$cohortTL < mymed,'Short telomeres','Long telomeres')


rm(mymed)

juv <- subset(juv,!(is.na(cenTL)))

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
juv$Helper[juv$Helper>1] <- 1




# TQ and insects ----------------------------------------------------------


terr <- terr[complete.cases(terr),] #Get rid of blank rows
insects <- subset(insects,FieldPeriodID != 26)
dens$SPsize <- with(dens,(PsizeFP-mean(PsizeFP))/sd(PsizeFP))


juv$TQ <- NA
juv$TQI <- NA
juv$cenTQ <- NA
juv$Insect <- NA
juv$Density <- NA
juv$SDensity <- NA
juv$HatchDate <- NA


for(i in 1:nrow(juv))
{
  if(juv$TerritoryID[i] %in% terr$TerritoryID)
  {
    cd <- subset(terr,TerritoryID == juv$TerritoryID[i])
    if(juv$FieldPeriodID[i] %in% cd$FieldPeriodID)
    {
      juv$TQ[i] <- log10(subset(cd,FieldPeriodID == juv$FieldPeriodID[i])$TQcorrected)
    } else
    {
      log10(mean(cd$TQcorrected))
    }
    juv$TQI[i] <- juv$TQ[i]/juv$GroupSize[i]
    
  }
  
  
  
  if(juv$FieldPeriodID[i] %in% insects$FieldPeriodID)
  {
    juv$Insect[i] <- insects$MeanInsects[insects$FieldPeriodID == juv$FieldPeriodID[i]] 
  }
  
  
  
  if(juv$FieldPeriodID[i] %in% dens$FieldPeriodID)
  {
    juv$Density[i] <- dens$PsizeFP[dens$FieldPeriodID == juv$FieldPeriodID[i]]
    juv$SDensity[i] <- dens$SPsize[dens$FieldPeriodID == juv$FieldPeriodID[i]]
    
  }
  
  
  if(juv$BirdID[i] %in% chickinfo$BirdID)
  {
    juv$HatchDate[i] <- paste(chickinfo$HatchDate[chickinfo$BirdID == juv$BirdID[i]])
  }
  
  
}


juv$HatchDate <- as.Date(juv$HatchDate,"%d/%m/%Y")
juv$Chickage <- as.numeric(juv$CatchDate-juv$HatchDate)
juv$O_HatchDate <- as.POSIXlt(juv$HatchDate,format='%yyyy-%mm-%dd')$yday
juv$O_CatchDate <- as.POSIXlt(juv$CatchDate,format='%yyyy-%mm-%dd')$yday
juv <- subset(juv,!(is.na(Density)))



#Get rid of NAs and cross-fostered birds
juv <- subset(juv,!is.na(Tarsus))
juv <- subset(juv,!is.na(BodyMass))
juv <- subset(juv,!is.na(TQ))
juv <- subset(juv,Status!='XF')


juv$Condition <- lm(BodyMass~Tarsus+Age+Sex,data=juv)$residuals

# Look at telomere loss ---------------------------------------------------
#Get earliest catch for each juvenile
earlies <- merge(juv,aggregate(CatchDate~BirdID,juv,min)) %>% 
subset(!(duplicated(BirdID)))

#Get earliest subsequent catch
laters <- merge(aggregate(CatchDate~BirdID,subset(dd,!(BloodID %in% earlies$BloodID)),
                          mean),
                subset(dd,!(BloodID %in% earlies$BloodID))) %>% 
subset(!(duplicated(BirdID)))


earlies <- subset(earlies,BirdID %in% laters$BirdID)
laters <- subset(laters,BirdID %in% earlies$BirdID)

earlies <- earlies[order(earlies$BirdID),]
laters <- laters[order(laters$BirdID),]


#Creat data frame with TROC
Loss <- data.frame(earlies,
                   Loss = laters$RTL-earlies$RTL,
                   RTL1 = laters$RTL,
                   Agemonths1 = laters$Agemonths,
                   TimeDiff = as.numeric(laters$CatchDate-earlies$CatchDate),
                   RemainingLife2 = laters$RemainingLife,
                   SurvivedNext2 = laters$SurvivedNext)
Loss$TROC <- with(Loss,Loss/TimeDiff)

#
Loss <- subset(Loss,TimeDiff>(365/2))
Loss$TimeDiff <- Loss$TimeDiff/365

# Get rid of stuff not to be used -----------------------------------------

rm(status,helpers,hatchdate,x1,x2,x3,x4)


# Field period average data for plots -----------------------------------------------

juvseason <- ddply(juv,
                   .(factor(LayYear),Season),
                   summarize,
                   RTL = median(RTL),
                   TLse = se(RTL),
                   Helper = mean(Helper,na.rm=T),
                   Insect = mean(Insect,na.rm=T),
                   Density = mean(Density),
                   Lifespan = median(RemainingLife),
                   LayYear = mean(CatchYear),
                   Age = mean(Agemonths),
                   n = length(TQ),
                   TQ = mean(TQ,na.rm=T),
                   Survived=mean(SurvivedNext)*100)
juvseason <- subset(juvseason,n>4)

with(juv,tapply(Died,LayYear,mean))
juv9 <- subset(juv,LayYear<2008)
juv13 <- subset(juv,LayYear<2013)


juvseason9 <- ddply(juv9,
                   .(factor(LayYear),Season),
                   summarize,
                   RTL = median(RTL),
                   TLse = se(RTL),
                   Helper = mean(Helper,na.rm=T),
                   Insect = mean(Insect,na.rm=T),
                   Density = mean(Density),
                   Lifespan = median(RemainingLife),
                   LayYear = mean(CatchYear),
                   Age = mean(Agemonths),
                   n = length(TQ),
                   TQ = mean(TQ,na.rm=T),
                   Survived=mean(SurvivedNext)*100)
juvseason9 <- subset(juvseason9,n>4)







# Longitudinal data for whole dataset -------------------------------------

#First calculate Delta RTL 
temp <- subset(dd,!duplicated(paste0(BirdID,Agemonths)))
counts <- with(temp,tapply(RTL,BirdID,length))
c3 <- names(counts[counts >= 2])
dd3 <- subset(temp,BirdID %in% c3)
DeltaRTL <- NA

for(i in 1:nrow(dd3))
{
  cd <- subset(dd3,BirdID == dd3$BirdID[i])
  cd <- cd[order(cd$Agemonths),]
  
  birdpos <- which(cd$BloodID == dd3$BloodID[i])
  
  if(birdpos == 1)
  {
    dd3$DeltaRTL[i] <- 0
  } else
  {
    prevbird <- cd[birdpos-1,'RTL']
    dd3$DeltaRTL[i] <- dd3$RTL[i]-prevbird
    dd3$MidAge[i] <- mean(c(dd3$Agemonths[i],cd[birdpos-1,'Agemonths']))
    dd3$MinAge[i] <- cd[birdpos-1,'Agemonths']
    dd3$TimeDiff[i] <- dd3$Agemonths[i] - cd[birdpos-1,'Agemonths']
  }
  
  
}


temp <- subset(dd0,!duplicated(paste0(BloodID,RTL)))
temp <- subset(temp,RTL<2)
temp <- subset(temp,RTL > 0.04)

counts <- with(temp,tapply(RTL,BloodID,length))
c3 <- names(counts[counts %in% c(2:5)])


dd3_2 <- subset(temp,BloodID %in% c3)
dd3_2 <- dd3_2[order(dd3_2$BloodID),]
counts <- with(dd3_2,tapply(RTL,BloodID,length))

dd3_2$birdpos <- unlist(sapply(counts,function(x) c(1:x)))

dd3_2$DeltaRTL <- NA

for(i in 1:nrow(dd3_2))
{
  cd <- subset(dd3_2,BloodID == dd3_2$BloodID[i])
  
  if(dd3_2$birdpos[i] == 1)
  {
    dd3_2$DeltaRTL[i] <- 0
  } else
  {
    prevbird <- subset(cd,birdpos == dd3_2$birdpos[i]-1)$RTL
    dd3_2$DeltaRTL[i] <- dd3_2$RTL[i]-prevbird
  }
  
  
}


dd3 <- subset(dd3,DeltaRTL !=0)
dd3_2 <- subset(dd3_2,DeltaRTL!=0)
dd3_2 <- dd3_2[1:nrow(dd3),]


ddL <- data.frame(DeltaRTL = c(dd3$DeltaRTL,dd3_2$DeltaRTL),Group = rep(c('Among samples','Within samples'),c(nrow(dd3),nrow(dd3_2))))


rm(temp,dd2,dd3_2)
