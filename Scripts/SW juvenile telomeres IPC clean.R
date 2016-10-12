#################################################################################*
## Early life telomere dynamics and late-life fitness in a wild bird population
## CLEAN THE DATA
#################################################################################

#Only use Ellie's samples
dd0 <- subset(dd0,Whodunnit == 'EAF')

#Average repeats of blood samples
dd0$RTL2 <- ave(dd0$RTL,c(dd0$BloodID,dd0$Status,dd0$PlateID))
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
dd$CatchMonth <- as.numeric(substr(dd$CatchDate,6,7))


dd$Season <- ifelse(as.numeric(format(dd$CatchDate,'%m')) %in% c(4:10),
                    'Major','Minor')
dd <- subset(dd,Season == 'Major')
dd <- subset(dd,CatchMonth %in% c(6:9))
dd$yday <- as.POSIXlt(dd$CatchDate)$yday - min(as.POSIXlt(dd$CatchDate)$yday)

# Age data ----------------------------------------------------------------

dd$Age <- (dd$CatchYear-dd$LayYear)

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
                              ifelse(dd$Ageclass == 'SA',10,dd$Age*12+6)))


  #dd$Age[dd$Fledged != 'Adults'] <- 1




# Tarsus and delta age ------------------------------------------------------------------

dd$Tarsus <- NA
dd$DeltaAge <- NA

for(i in 1:nrow(dd))
{
  ifelse(is.na(dd$RightTarsus[i]),
         ifelse(is.na(dd$LeftTarsus[i]),
                dd$Tarsus[i] <- NA,
                dd$Tarsus[i] <- dd$LeftTarsus[i]),
         dd$Tarsus[i] <- dd$RightTarsus[i])
  
  dd$DeltaAge[i] <- dd$Agemonths[i] - mean(subset(dd,BirdID == dd$BirdID[i])$Agemonths)
}

dd$RightTarsus <- NULL
dd$LeftTarsus <- NULL



# Remove unwanted data/outliers ----------------------------------------------------


dd <- subset(dd,RTL > 0.05)
dd <- subset(dd,CqTelomere <28)
dd <- subset(dd,CqGAPDH < 27)
dd <- subset(dd,RTL<3)

dd <- subset(dd,BodyMass>5)
dd <- subset(dd,Tarsus>17)




# Survival and lifespan ---------------------------------------------------

dd$RemainingLife <- round(as.numeric(dd$DeathDate-dd$CatchDate)/365,0)
dd$SurvivedNext <- ifelse(dd$RemainingLife<1,0,1)

dd$Lifespan <- (dd$DeathYear-dd$LayYear)+1
dd$Died <- ifelse(dd$DeathYear<2014,1,0)


# Sex ---------------------------------------------------------------------

dd$Sex <- ifelse(dd$SexEstimate == 1,'Males','Females')






#Centre telomere length by birth year
dd$cenTL <- NA
dd$cohortTL <- NA
for(i in 1:nrow(dd))
{
  currentdata <- subset(dd,FieldPeriodID == FieldPeriodID[i])
  if(nrow(currentdata>4))
  {
    dd$cenTL[i] <- (dd$RTL[i] - mean(currentdata$RTL))/sd(currentdata$RTL)
    dd$cohortTL[i] <- median(currentdata$RTL)  
  }

}

mymed <- median(dd$cenTL,na.rm=T)
dd$cenTLF <- ifelse(dd$cenTL < mymed,'Short telomeres','Long telomeres')

mymed <- median(dd$cohortTL,na.rm=T)
dd$coTLF <- ifelse(dd$cohortTL < mymed,'Short telomeres','Long telomeres')


rm(mymed)

# Helpers and social group size -------------------------------------------

status <- subset(status,BreedGroupID %in% dd$BreedGroupID)

for(i in 1:nrow(dd))
{
  currentBG <- dd$BreedGroupID[i]
  currentdata <- subset(subset(status,BreedGroupID == currentBG),BirdID !=dd$BirdID[i])
  dd$GroupSize[i] <- nrow(currentdata)
  dd$Helper[i] <- nrow(subset(currentdata,Status == 'H'))
  dd$NonHelper[i] <- nrow(subset(currentdata,Status %in% c('AB','ABX')))
}

dd <- subset(dd,GroupSize>0)
dd$Helper[dd$Helper>1] <- 1




# TQ and insects ----------------------------------------------------------


terr <- terr[complete.cases(terr),] #Get rid of blank rows
insects <- subset(insects,FieldPeriodID != 26)
dens$SPsize <- with(dens,(PsizeFP-mean(PsizeFP))/sd(PsizeFP))
terr$LogTQ <- log10(terr$TQcorrected)

dd$TQ <- NA
dd$TQI <- NA
dd$Insect <- NA
dd$BirthInsect
dd$Density <- NA
dd$cenTQ <- NA
dd$Malaria <- NA

for(i in 1:nrow(dd))
{
  if(dd$TerritoryID[i] %in% terr$TerritoryID)
  {
    cd <- subset(terr,TerritoryID == dd$TerritoryID[i])
    if(dd$FieldPeriodID[i] %in% cd$FieldPeriodID)
    {
      dd$TQ[i] <- subset(cd,FieldPeriodID == dd$FieldPeriodID[i])$LogTQ
      cy <- subset(terr,FieldPeriodID == dd$FieldPeriodID[i])$LogTQ
      
      if(!is.na(dd$TQ[i]))
      {
        dd$cenTQ[i] <- (dd$TQ[i]-mean(cy))/sd(cy)
      }

              
    }
    dd$TQI[i] <- dd$TQ[i]/dd$GroupSize[i]
    
  }
  
  
  
  if(dd$FieldPeriodID[i] %in% insects$FieldPeriodID)
  {
    dd$Insect[i] <- insects$MeanInsects[insects$FieldPeriodID == dd$FieldPeriodID[i]] 
  }
  
  
  
  if(dd$FieldPeriodID[i] %in% dens$FieldPeriodID)
  {
    dd$Density[i] <- dens$PsizeFP[dens$FieldPeriodID == dd$FieldPeriodID[i]]
    dd$SDensity[i] <- dens$SPsize[dens$FieldPeriodID == dd$FieldPeriodID[i]]
    
  }
  
  if(dd$BloodID[i] %in% mal$BloodID)
  {
    dd$Malaria[i] <- mal$Consensus[which(mal$BloodID == dd$BloodID[i])]
  }
}


#Get rid of NAs and cross-fostered birds
dd <- subset(dd,!is.na(Tarsus))
dd <- subset(dd,!is.na(BodyMass))
dd <- subset(dd,Status!='XF')
dd$Condition <- lm(BodyMass~Tarsus,data = dd)$residuals 


# Get rid of stuff not to be used -----------------------------------------

rm(status,helpers,hatchdate,x1,x2,x3,x4)



# Longitudinal data for whole dataset -------------------------------------

#First calculate Delta RTL 
temp <- subset(dd,!duplicated(paste0(BirdID,Agemonths)))
counts <- with(temp,tapply(RTL,BirdID,length))
c3 <- names(counts[counts >= 2])
dd3 <- subset(temp,BirdID %in% c3)
dd3$TimeDiff <- NA
dd3$DeltaRTL <- NA
dd3$DeltaTL <- NA
dd3$DeltaGAP <- NA
dd3$BloodID1 <- NA
dd3$SurvivedNext2 <- NA

for(i in 1:nrow(dd3))
{
  cd <- subset(dd3,BirdID == dd3$BirdID[i])
  cd <- cd[order(cd$Agemonths),]
  
  birdpos <- which(cd$BloodID == dd3$BloodID[i])
  
  if(birdpos == nrow(cd))
  {
    dd3$DeltaRTL[i] <- NA
  } else
  {
    nextbird <- cd[birdpos+1,]
    dd3$DeltaRTL[i] <- nextbird$RTL - dd3$RTL[i]
    dd3$RTL1[i] <- nextbird$RTL
    dd3$DeltaGAP[i] <- nextbird$CqGAPDH - dd3$CqGAPDH[i]
    dd3$DeltaTL[i] <- nextbird$CqTelomere - dd3$CqTelomere[i]
    dd3$TimeDiff[i] <- cd[birdpos+1,'Agemonths'] - dd3$Agemonths[i]
    dd3$BloodID1[i] <- nextbird$BloodID
    dd3$SurvivedNext1[i] <- nextbird$SurvivedNext
    dd3$RemainingLife1[i] <- nextbird$RemainingLife
    dd3$Died1[i] <- nextbird$Died
    dd3$DeltaCondition[i] <- nextbird$Condition - dd3$Condition[i]
    dd3$DeltaAge[i] <- nextbird$Agemonths - dd3$Agemonths[i]
    dd3$DeltaInsect[i] <- nextbird$Insect - dd3$Insect[i]
    dd3$Agemonths1[i] <- nextbird$Agemonths
  }
  
  
}



temp <- subset(dd0,!duplicated(paste0(BloodID,RTL)))
temp <- subset(temp,RTL<2)
temp <- subset(temp,RTL > 0.05)
#temp <- subset(temp,CqTelomere <22)
temp <- subset(temp,CqGAPDH < 25)

counts <- with(temp,tapply(RTL,BloodID,length))
c3 <- names(counts[counts %in% c(2:5)])


dd3_2 <- subset(temp,BloodID %in% c3)
dd3_2 <- dd3_2[order(dd3_2$BloodID),]
counts <- with(dd3_2,tapply(RTL,BloodID,length))

dd3_2$birdpos <- unlist(sapply(counts,function(x) c(1:x)))

dd3_2$DeltaRTL <- NA
dd3_2$DeltaTL <- NA
dd3_2$DeltaGAP <- NA
dd3_2$BloodID1 <- NA

for(i in 1:nrow(dd3_2))
{
  cd <- subset(dd3_2,BloodID == dd3_2$BloodID[i])
  
  if(dd3_2$birdpos[i] == nrow(cd))
  {
    dd3_2$DeltaRTL[i] <- NA
  } else
  {
    nextbird <- subset(cd,birdpos == dd3_2$birdpos[i]+1)
    dd3_2$DeltaRTL[i] <- nextbird$RTL - dd3_2$RTL[i]
    dd3_2$RTL1[i] <- nextbird$RTL
    dd3_2$BloodID1[i] <- nextbird$BloodID
  }
  
  
}


dd3 <- subset(dd3,!is.na(DeltaRTL))
dd3_2 <- subset(dd3_2,!is.na(DeltaRTL))

dd3$DeltaRTLF <- ifelse(dd3$DeltaRTL>0,1,0)


ddL <- data.frame(BirdID = c(dd3$BirdID,dd3_2$BirdID),
                  BloodID = c(dd3$BloodID,dd3_2$BloodID),
                  BloodID1 = c(dd3$BloodID1,dd3_2$BloodID1),
                  RTL = c(dd3$RTL,dd3_2$RTL),
                  RTL1 = c(dd3$RTL1,dd3_2$RTL1),
                  DeltaRTL = c(dd3$DeltaRTL,dd3_2$DeltaRTL),
                  DeltaTL = c(dd3$DeltaTL,dd3_2$DeltaTL),
                  DeltaGAP = c(dd3$DeltaGAP,dd3_2$DeltaGAP),
                  Group = rep(c('Among samples','Within samples'),c(nrow(dd3),nrow(dd3_2))))

rm(temp,dd3_2)





# Subset birds with juvenile samples --------------------------------------

juvdata <- subset(dd,Agemonths<12)
addata <- subset(dd,Agemonths>11)
juvdata <- juvdata[order(juvdata$Agemonths),]


addata <- subset(addata, BirdID %in% juvdata$BirdID)


toreplace <- c('TQ','Helper','NonHelper','GroupSize','Insect','Density','BodyMass','Tarsus')

for(i in 1:nrow(addata))
{
  cd <- subset(juvdata,BirdID == addata$BirdID[i])[1,toreplace]
  addata[i,toreplace] <- cd
}

juv_r <- rbind(juvdata,addata)
juv <- juvdata

juv14 <- subset(juv_r,LayYear<2014)
juv12 <- subset(juv_r,LayYear<2012)

