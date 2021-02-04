# this script will explore the invert data

library(tidyverse)
if(!require(vegan))install.packages("vegan"); library(vegan)

# bring in data
source("scripts/03_reimport.R")

i2<-inverts%>%
  pivot_wider(names_from=taxa,values_from=abundance,values_fill=0)

i.env<-i2[,1:3]
i.com<-i2[,-1:-3]

i.env<-i.env%>%
  mutate(sampling=case_when(
    sampling==1~0,
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17),
    season=case_when(
      sampling==0~"summer",
      sampling==1~"summer",
      sampling==5~"winter",
      sampling==12~"summer",
      sampling==17~"winter"))
i.env$plot<-as.factor(i.env$plot)
# NMDS
i.com$dummy<-1
com.dist<-vegdist(i.com,"bray")

i.mds<-metaMDS(com.dist,trymax = 100)

# rda
i.com.hel<-decostand(i.com,"hellinger")

i.rda<-rda(i.com.hel~.,data=i.env)
i.rda2<-rda(i.com.hel~1,i.env)
summary(i.rda)
RsquareAdj(i.rda)
ordistep(i.rda2,scope = formula(i.rda),direction = "forward")

anova(i.rda)
plot(i.rda)
i.scores<-data.frame(scores(i.rda,1:3)$sites)%>%
  bind_cols(i.env)

ggplot(data=i.scores)+
  geom_point(aes(x=PC2,y=PC3,
                 color=treatment),size=2)+
  scale_color_viridis_d(option="B",end=.8)+
  facet_wrap(~sampling)



#fish
f2<-fish%>%
  filter(abundance!=0)%>%
  pivot_wider(names_from=taxa,values_from=abundance,values_fill=0)#%>%
#  mutate(dummy=1)

f.env<-f2[,1:3]
f.com<-f2[,-1:-3]

f.env<-f.env%>%
  mutate(sampling=case_when(
    sampling==1~0,
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17),
    season=case_when(
      sampling==0~"summer",
      sampling==1~"summer",
      sampling==5~"winter",
      sampling==12~"summer",
      sampling==17~"winter"))
f.env$plot<-as.factor(f.env$plot)

# NMDS
com.pa<-decostand(f.com,"pa")
com.dist<-vegdist(com.pa,"bray")

f.mds<-metaMDS(com.dist,trymax = 100)
plot(f.mds)
f.mds.scores<-data.frame(scores(f.mds))%>%
  bind_cols(f.env)

ggplot(data=f.mds.scores)+
  geom_point(aes(x=NMDS1,y=NMDS2,
                 color=treatment),size=2,alpha=.5)+
  scale_color_viridis_d(option="B",end=.8)+
  facet_wrap(~sampling)

# rda
f.com.hel<-decostand(f.com,"hellinger")

f.rda<-rda(f.com.hel)
i.rda2<-rda(i.com.hel~1,i.env)
summary(i.rda)
RsquareAdj(i.rda)
ordistep(i.rda2,scope = formula(i.rda),direction = "forward")

anova(f.rda)
plot(f.rda)
f.scores<-data.frame(scores(f.rda,1:3)$sites)%>%
  bind_cols(f.env)
plot(f.rda)
ggplot(data=f.scores)+
  geom_point(aes(x=PC1,y=PC2,
                 color=treatment),size=2)+
  scale_color_viridis_d(option="B",end=.8)#+
  facet_wrap(~sampling)
