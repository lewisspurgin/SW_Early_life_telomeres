#################################################################################*
## Early life telomere dynamics and late-life fitness in a wild bird population
## CLEAN THE DATA
#################################################################################

#Only use Ellie's samples
dd0 <- subset(dd0,Whodunnit == 'EAF')

#Average repeats of blood samples
av <- ave(dd0$TL,c(dd0$BloodID,dd0$Status,dd0$PlateID))
dd0$TL <- av
dd <- dd0[!(duplicated(dd0$BloodID)),]


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


dd <- droplevels(subset(subset(dd,TL>50),TL<30000))
dd <- subset(dd,BodyMass>5)
dd <- subset(dd,Tarsus>17)


# Telomere data -----------------------------------------------------------

dd$TLKB <- dd$TL/1000
dd$LogTL <- log10(dd$TL)




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

dd <- subset(dd,Agemonths>0)

dd$Age[dd$Fledged != 'Adults'] <- 1

# Survival and lifespan ---------------------------------------------------

dd$RemainingLife <- dd$DeathYear-dd$CatchYear
dd$SurvivedNext <- ifelse(dd$RemainingLife>1,1,0)
dd$Lifespan <- (dd$DeathYear-dd$LayYear)+1
dd$Died <- ifelse(dd$DeathYear<2013,1,0)



# Exclude winter seasons and early years ----------------------------------

dd <- subset(dd,LayYear>1997)
dd <- subset(dd,Season == 'Summer')


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
  juv$cenTL[i] <- (juv$LogTL[i] - mean(currentdata$LogTL))/sd(currentdata$LogTL)
  juv$cohortTL[i] <- mean(currentdata$LogTL)
}

mymed <- median(juv$cenTL,na.rm=T)
juv$cenTLF <- ifelse(juv$cenTL < mymed,'Short telomeres','Long telomeres')

mymed <- median(juv$cohortTL,na.rm=T)
juv$coTLF <- ifelse(juv$cohortTL < mymed,'Short telomeres','Long telomeres')


rm(mymed)
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
subset(terr,FieldPeriodID == 105)

insects <- subset(insects,FieldPeriodID != 26)

juv$TQ <- NA
juv$TQI <- NA
juv$cenTQ <- NA
juv$Insect <- NA
juv$Density <- NA

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
      juv$TQ[i] <- log10(mean(cd$TQcorrected))
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
  }
  
  
}





# Subset Fledglings and subadults -----------------------------------------------------------



#Get rid of NAs and cross-fostered birds
juv <- subset(juv,!is.na(Tarsus))
juv <- subset(juv,!is.na(BodyMass))
juv <- subset(juv,!is.na(TQ))
juv <- subset(juv,Status!='XF')


#Subsetting
FlSA <- subset(juv,Ageclass!='CH')
chicks <- subset(juv,Ageclass == 'CH')

juvall <- droplevels(juv[complete.cases(juv),])
FlSAall <- droplevels(FlSA[complete.cases(FlSA),])
chickall <- droplevels(chicks[complete.cases(chicks),])



# Look at telomere loss ---------------------------------------------------


#Get earliest catch for each juvenile
x3 <- aggregate(CatchDate~BirdID,juv,min)
x4 <- merge(x3,juv)
x4 <- x4[!(duplicated(x4$BirdID)),]

#Get earliest subsetquent catch
laters <- subset(dd,!(BloodID %in% x4$BloodID))
x1 <- aggregate(CatchDate~BirdID,laters,min)
x2 <- merge(laters,x1)
x2 <- x2[!(duplicated(x2$BirdID)),]


x2 <- x2[x2$BirdID %in% x4$BirdID,]
x4 <- x4[x4$BirdID %in% x2$BirdID,]

x3 <- x2[order(x2$BirdID),]
x4 <- x4[order(x4$BirdID),]


#Apply Verhulst's correction for regresison to the mean
xx1 <- x4$TLKB-mean(x4$TLKB)
xx2 <- x2$TLKB-mean(x2$TLKB)          
rho <- cor(x4$TLKB,x2$TLKB)

#Creat data frame with TROC
Loss <- data.frame(x4,
                   Loss = x4$LogTL-x2$LogTL,
                   LogTL1 = x2$LogTL,
                   Agemonths1 = x2$Agemonths,
                   TimeDiff = as.numeric(x2$CatchDate-x4$CatchDate)/365,
                   RemainingLife2 = x2$RemainingLife,
                   D=(rho*xx1)-xx2)
Loss$TROC <- with(Loss,Loss/TimeDiff)

#Restrict to brids with sample within 5 years
Loss <- subset(Loss,TimeDiff<5.5)

#Calculate body condition
Loss$Condition <- lm(BodyMass~Tarsus, data=Loss)$residuals

#Remove outliers
Loss <- subset(subset(Loss,TROC> -10),TROC < 10)

# Get rid of stuff not to be used -----------------------------------------

rm(status,helpers,hatchdate,x1,x2,x3,x4)




# Field period average data for plots -----------------------------------------------

juvseason <- ddply(FlSA,
                   .(FieldPeriodID,Season),
                   summarize,
                   TL = median(LogTL),
                   TLse = se(LogTL),
                   Helper = mean(Helper),
                   Insect = median(Insect),
                   Density = mean(Density),
                   Lifespan = median(RemainingLife),
                   LayYear = mean(CatchYear),
                   Age = mean(Agemonths),
                   n = length(TLKB),
                   TQ = mean(TQ))
juvseason <- subset(juvseason,n>4)


juv9 <- subset(juv,LayYear<2009)
juv13 <- subset(juv,LayYear<2013)
