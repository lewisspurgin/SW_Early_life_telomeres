#################################################################################*
## Telomere dynamics in a wild bird population
## LOAD MODELS
#################################################################################*

# RTL and storage time
RTL_ST <- lmer(RTL ~ StorageTime + (1|BirdID) + (1|PlateID),
                data=subset(subset(dd,Agemonths > 1.5),Whodunnit == "EAF"),
                REML = FALSE)

# Telomere length, age and cohort -----------------------------------------

RTL_Age <- lmer(RTL ~ Cohort + Agemonths + LogAge + FAge + Age2 + (1|BirdID) + (1|CatchYear) + (1|PlateID),
                data=dd,
                REML = FALSE)
RTL_LogAge <- lmer(RTL ~ LogAge + Cohort + (1|BirdID) + (1|CatchYear) + (1|PlateID),
                data=dd,
                REML = FALSE)
RTL_LogAgeXcohort <- lmer(RTL ~ LogAge*Cohort + (1|BirdID) + (1|CatchYear) + (1|PlateID),
                   data=dd,
                   REML = FALSE)
RTL_DeltaAge <- lmer(RTL ~ MeanAge + Cohort + DeltaAge + DeltaLogAge + DeltaAge2 + (1|BirdID) + (1|CatchYear),
                data=dd3,
                REML = FALSE)
RTL_DeltaLogAge <- lmer(RTL ~ DeltaLogAge + (1|BirdID) + (1|CatchYear),
                   data=dd3,
                   REML = FALSE)
DeltaRTL_LogAge <- lmer(DeltaRTL ~ LogAge + (1|BirdID) + (1|CatchYear),
                        data=subset(dd3,TimeDiff < 3.5),
                        REML = FALSE)
DeltaRTLF_LogAge <- glmer(DeltaRTLF ~ LogAge + (1|BirdID) + (1|CatchYear),
                        data=subset(dd3,TimeDiff < 3.5),
                        family = "binomial")
RTL_RTL1 <- lmer(RTL1 ~ RTL + LogAge + (1|BirdID) + (1|LayYear),
                        data=dd3,
                        REML = FALSE)

# RTL and ecology ---------------------------------------------------------

ddNA <- subset(subset(dd,!is.na(Insect)),!is.na(TQ))
dd3NA <- subset(subset(dd3,!is.na(Insect)),!is.na(TQ))

RTL_Full <- lmer(RTL ~  LogAge + Tarsus + Helper + GroupSize + Sex + BodyMass + Insect + Density + TQ + (Agemonths|LayYear) + (1|BirdID),
                 data = ddNA,
                 REML = FALSE)
RTL_Full_NoTarsus <- lmer(RTL ~  LogAge + Helper + GroupSize + Sex + BodyMass + Insect + Density + TQ + (Agemonths|LayYear) + (1|BirdID),
                          data = ddNA,
                          REML = FALSE)
RTL_Full_SexTarsus <- lmer(RTL ~  LogAge + Helper + GroupSize + Sex + BodyMass + Insect + Density + Sex*Tarsus + TQ + (Agemonths|LayYear) + (1|BirdID),
                          data = ddNA,
                          REML = FALSE)
DeltaRTL_Full <- lmer(DeltaRTL ~  LogAge + Tarsus + Helper + GroupSize + Sex + BodyMass + Insect + Density + TQ + (1|BirdID),
                      data = dd3NA,
                             REML = FALSE)