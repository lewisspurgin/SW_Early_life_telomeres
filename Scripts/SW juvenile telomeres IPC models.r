
RTL_Age <- lmer(RTL~Agemonths + (1|BirdID),data = dd)
RTL_LogAge <- lmer(RTL~LogAge + (1|BirdID),data = dd)
RTL_LogAge_LY <- lmer(RTL~LogAge + (1|BirdID) + (1|LayYear),data = dd)
RTL_LogAge_AgeLY <- lmer(RTL~LogAge + (1|BirdID) + (Agemonths|LayYear),data = dd)
RTL_LogAge_CY <- lmer(RTL~LogAge + (1|BirdID) + (1|CatchYear),data = dd)
RTL_LogAge_AgeCY <- lmer(RTL~LogAge + (1|BirdID) + (Agemonths|CatchYear),data = dd)