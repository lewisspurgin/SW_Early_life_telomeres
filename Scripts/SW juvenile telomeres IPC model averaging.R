#################################################################################*
## Telomere dynamics in a wild bird population
## PERFORM MODEL AVERAGING
#################################################################################*


RTL_Full <- standardize(RTL_Full)
RTL_Full_Dredge <- subset(dredge(RTL_Full),delta <= 6)

DeltaRTL_Full <- standardize(DeltaRTL_Full)
DeltaRTL_Full_Dredge <- subset(dredge(DeltaRTL_Full),delta <= 6)