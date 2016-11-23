#################################################################################*
## Telomere dynamics in a wild bird population
## LOAD MODELS
#################################################################################*


# Telomere length, age and cohort -----------------------------------------

RTL_Age <- lmer(RTL ~ Cohort*(Agemonths + LogAge + FAge) + (1|BirdID) + (1|CatchYear) + (1|PlateID),
                data=dd,
                REML = FALSE)

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


