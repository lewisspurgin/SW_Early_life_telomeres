

juv2 <- juv

temp <- rbind(juv2[,c('BirdID','LogTL','Agemonths','Fledged','SurvivedNext','Died')],adults[,c('BirdID','LogTL','Agemonths','Fledged','SurvivedNext','Died')])
temp$Agemonths <- factor(temp$Agemonths)
temp <- subset(temp, Died ==1)

output <- mat.or.vec(nlevels(temp$Agemonths),5)
ss <- rep(NA,nlevels(temp$Agemonths))

for(i in 1:nlevels(temp$Age))
{
  d1 <- subset(temp,Agemonths == levels(temp$Agemonths)[i+1])
  d2 <- subset(temp,Agemonths == levels(temp$Agemonths)[i])
  
  d1 <- d1[!(duplicated(d1$BirdID)),]
  d2 <- d2[!(duplicated(d2$BirdID)),]
  
  output[i,1] <- paste0(levels(temp$Agemonths)[i],' - ',levels(temp$Agemonths)[i+1])
  output[i,2] <- mean(d1$LogTL)-mean(d2$LogTL)
  
  within1 <- subset(d1,BirdID %in% d2$BirdID)
  
  within2 <- subset(d2,BirdID %in% d1$BirdID)
  ss[i] <- nrow(within2)
  
  within <- merge(within2,within1,by = 'BirdID')
  output[i,3] <- mean((within$LogTL.y-within$LogTL.x))
  
  disap <- subset(d2,!(BirdID %in% d1$BirdID))
  
  
  output[i,4] <- mean(d1$LogTL)-mean(disap$LogTL)
  
  ap <- subset(d1,!(BirdID %in% d2$BirdID))
  output[i,5] <- mean(d1$LogTL)-mean(ap$LogTL)
  
}


ss
output <- output[complete.cases(output),]
output <- output[1:7,]

plot(output[,2],type = 'l',ylim=c(-0.5,0.3),axes=F)
points(output[,3],type = 'l',col='red',lty=2)
points(output[,4],type = 'l',col='darkgrey',lty=2)
points(output[,5],type = 'l',col='blue',lty=2)
abline(h=0,lty=3)

axis(1,at=c(1:nrow(output)),labels=output[,1])
axis(2)


fa <- ggplot(juvseason,aes(x = Insect,y = TL))+
  geom_jitter(col = 'grey') +
  ylab('Median telomere length') + 
  xlab('Food availability') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 16,vjust=0.8),
        axis.title.x = element_text(size = 16),
        legend.position = 'none')+
  geom_point(col = 'grey')+
  stat_smooth(method = 'lm')

fb <- ggplot(juvseason,aes(x = TQ,y = TL))+
  geom_jitter(col = 'grey') +
  ylab('Median telomere length') + 
  xlab('Territory quality') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 16,vjust=0.8),
        axis.title.x = element_text(size = 16),
        legend.position = 'none')+
  geom_point(col = 'grey')+
  stat_smooth(method = 'lm')

fc <- ggplot(juvseason,aes(x = Helper,y = TL))+
  geom_jitter(col = 'grey') +
  ylab('Median telomere length') + 
  xlab('Number of helpers') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 16,vjust=0.8),
        axis.title.x = element_text(size = 16),
        legend.position = 'none')+
  geom_point(col = 'grey')+
  stat_smooth(method = 'lm')

fd <- ggplot(juvseason,aes(x = Density,y = TL))+
  geom_jitter(col = 'grey') +
  ylab('Median telomere length') + 
  xlab('Population density') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 16,vjust=0.8),
        axis.title.x = element_text(size = 16),
        legend.position = 'none')+
  geom_point(col = 'grey')+
  stat_smooth(method = 'lm')

multiplot(fa,fb,fc,fd,layout=matrix(1:4,2,2))
