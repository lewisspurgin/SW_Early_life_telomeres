#################################################################################*
## Telomere dynamics in a wild bird population
## LOAD MODELS
#################################################################################*


# Telomere length, age and cohort -----------------------------------------

RTL_Age <- lmer(RTL~Agemonths + (1|BirdID) + (1|PlateID),data = dd)
RTL_LogAge <- lmer(RTL~LogAge + (1|BirdID) + (1|PlateID),data = dd)
RTL_LogAge_LY <- lmer(RTL~LogAge + (1|BirdID) + (1|PlateID) + (1|LayYear),data = dd)
RTL_LogAge_CY <- lmer(RTL~LogAge + (1|BirdID) + (1|PlateID) + (1|CatchYear),data = dd)
RTL_LogAge_AgeLY <- lmer(RTL~LogAge + (1|BirdID) + (1|PlateID) + (Agemonths|LayYear),data = dd)
RTL_LogAge_AgeCY <- lmer(RTL~LogAge + (1|BirdID) + (1|PlateID) + (Agemonths|CatchYear),data = dd)

RTL_DeltaAge <- lmer(RTL~DeltaAge + (1|BirdID),data = dd3)
RTL_DeltaAge_LY <- lmer(RTL~DeltaAge + (1|BirdID) + (1|LayYear),data = dd3)
RTL_DeltaAge_CY <- lmer(RTL~DeltaAge + (1|BirdID) + (1|CatchYear),data = dd3)
RTL_DeltaAge_AgeLY <- lmer(RTL~DeltaAge + (1|BirdID) + (DeltaAge|LayYear),data = dd3)
RTL_DeltaAge_AgeCY <- lmer(RTL~DeltaAge + (1|BirdID) + (DeltaAge|CatchYear),data = dd3)


# RTL and ecology ---------------------------------------------------------

ddNA <- subset(subset(dd,!is.na(Insect)),!is.na(TQ))
dd3NA <- subset(subset(dd3,!is.na(Insect)),!is.na(TQ))

RTL_Full <- lmer(RTL ~  LogAge + Tarsus + Helper + GroupSize + Sex + BodyMass + Insect + Density + TQ + (Agemonths|CatchYear) + (1|BirdID) + (1|PlateID),
                 data = ddNA,
                 REML = FALSE)
RTL_Full_NoTarsus <- lmer(RTL ~  LogAge + Helper + GroupSize + Sex + BodyMass + Insect + Density + TQ + (Agemonths|CatchYear) + (1|BirdID) + (1|PlateID),
                          data = ddNA,
                          REML = FALSE)
RTL_Full_SexTarsus <- lmer(RTL ~  LogAge + Helper + GroupSize + Sex + BodyMass + Insect + Density + Sex*Tarsus + TQ + (Agemonths|CatchYear) + (1|BirdID) + (1|PlateID),
                          data = ddNA,
                          REML = FALSE)
DeltaRTL_Full <- lmer(DeltaRTL ~  LogAge + Tarsus + Helper + GroupSize + Sex + BodyMass + Insect + Density + TQ + (1|CatchYear) + (1|BirdID),
                      data = dd3NA,
                             REML = FALSE)


