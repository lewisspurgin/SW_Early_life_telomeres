

pdf('Figures/Figs for Dave.pdf',
    width = 5,
    height = 5,
    paper = 'a4')


temp <- dd3
ggplot(ddL,aes(col=Group,fill = Group,y = ..scaled..,x = DeltaRTL))+
  geom_density(alpha = 0.2)+
  theme_classic()+
  ylab ('Scaled density')+
  xlab(expression(Delta*'RTL'))+
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 14,vjust=0.8),
        axis.title.x = element_text(size = 14),
        legend.position = 'none',
        legend.title = element_blank(),
        axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black')) +
  scale_colour_manual(values = c('darkgrey','gold'))+
  scale_fill_manual(values = c('darkgrey','gold'))+
  geom_vline(xintercept = 0,lty = 2)


temp$TimeDiff2 <- temp$TimeDiff/12
ggplot(temp,aes(y = DeltaRTL,x = TimeDiff2))+
  geom_point(col = 'grey')+
  theme_classic()+ 
  xlab ('Follow up time (years)')+
  ylab(expression(Delta*'RTL'))+
  ylim(c(-2,2))+
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 14,vjust=0.8),
        axis.title.x = element_text(size = 14),
        legend.position = c('none'),
        legend.title = element_blank(),
        axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black')) +
  geom_abline(slope = 0,intercept = 0,lty = 2) +
  stat_smooth(method = 'lm',col = 'darkgrey')





temp <- subset(dd,Age<13)
reps <- tapply(temp$RTL,temp$BirdID,length)
keep <- names(reps[reps>0])
temp <- subset(temp,BirdID %in% keep)
temp <- subset(temp,!duplicated(paste0(BirdID,Agemonths)))

mins <- tapply(temp$Agemonths,temp$BirdID,min)
keep2 <- names(mins[mins <400])
temp <- subset(temp,BirdID %in% keep2)

temp$Agemonths <- temp$Agemonths/12


ggplot(temp,aes(x = Agemonths,y = RTL,group = BirdID,col=factor(BirdID)))+
  geom_point(alpha = 1,col='grey')+
  geom_line(col = 'darkgrey',alpha = 0.2)+
  xlab('Age (years)') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 14,vjust=0.8),
        axis.title.x = element_text(size = 14),
        axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black'),
        legend.position = 'none')+
  geom_smooth(col = 'black',fill = 'black',lwd = 0.5,aes(group = NULL),method = 'lm',formula = y~log(x))


ggplot(temp,aes(x = Agemonths,y = RTL,group = BirdID,col=factor(BirdID)))+
  geom_point(alpha = 0.5,aes(col=Sex))+
  geom_line(aes(col = Sex),alpha = 0.2)+
  xlab('Age (years)') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 14,vjust=0.8),
        axis.title.x = element_text(size = 14),
        axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black'),
        legend.position = 'none')+
  geom_smooth(lwd = 1,aes(group = Sex,col = Sex,fill = Sex),method = 'lm',formula = y~log(x))


ggplot(temp,aes(x = Agemonths,y = RTL,col=Sex,fill = Sex))+
  geom_point(alpha = 0.5)+
  xlab('Age (years)') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 14,vjust=0.8),
        axis.title.x = element_text(size = 14),
        axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black'),
        legend.position = 'none')+
  geom_smooth(lwd = 1,method = 'lm',formula = y~log(x))


ggplot(temp,aes(x = Agemonths,y = RTL,col=Sex,fill=Sex))+
  xlab('Age (years)') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 14,vjust=0.8),
        axis.title.x = element_text(size = 14),
        axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black'),
        legend.position = 'none')+
  geom_smooth(lwd = 1,method = 'lm',formula = y~log(x))




ys <- names(which(table(subset(dd,Ageclass == 'CH')$LayYear)>10))
temp <- subset(dd,LayYear %in% ys)
temp <- subset(temp,LayYear != 1999)
temp$Cohort <- factor(temp$LayYear)
temp$Agemonths2 <- temp$Agemonths/12

ggplot(temp,aes(x = Agemonths2,y = RTL,col = Cohort))+
  xlab('Age (years)') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 14,vjust=0.8),
        axis.title.x = element_text(size = 14),
        axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black'),
        legend.position = 'top')+
  geom_smooth(method = 'lm',formula = y~log(x),se = F)




temp <- dd3
ggplot(temp,
       aes(x = RTL,
           y = RTL1)) +
  ylab('RTL time t+1') + 
  xlab('RTL time t') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 14,vjust=0.8),
        axis.title.x = element_text(size = 14),
        legend.position = c('none'),
        legend.title = element_text(size = 12),
        axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black')) +
  geom_point(col = 'grey')+
  geom_smooth(method = 'lm',col = 'darkgrey')


temp <- dd
ggplot(temp,
       aes(
         y = RTL,
         x=Insect)) +
  ylab('RTL') + 
  xlab('Insect abundance') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 14,vjust=0.8),
        axis.title.x = element_text(size = 14),
        axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black'),
        legend.position = 'none') +
  geom_boxplot(aes(group = Insect),outlier.shape = NA)+
  stat_smooth(method = 'lm',col = 'darkgrey')

temp <- subset(dd,DeltaAge != 0)
temp$DeltaAge2 <- temp$DeltaAge/12
ggplot(temp,
       aes(
         y = RTL,
         x=DeltaAge2)) +
  ylab('RTL') + 
  xlab('Delta Age (years)') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 14,vjust=0.8),
        axis.title.x = element_text(size = 14),
        axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black'),
        legend.position = 'none') +
  geom_point(col = 'grey')+
  stat_smooth(method = 'lm',col = 'darkgrey')


dev.off()
