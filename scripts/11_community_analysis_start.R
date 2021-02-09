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
      sampling==17~"winter"),
    yr=case_when(
      sampling==0~0,
      sampling==1~1,
      sampling==5~1,
      sampling==12~2,
      sampling==17~2))
#i.env$plot<-as.factor(i.env$plot)
# NMDS
i.com$dummy<-1
com.dist<-vegdist(i.com,"bray")

i.mds<-metaMDS(com.dist,trymax = 100)

plot(i.mds)
# rda
i.com.hel<-decostand(i.com,"hellinger")

i.pca<-rda(i.com.hel)
i.scores<-data.frame(scores(i.pca,1:3)$sites)%>%
  bind_cols(i.env)

ggplot(data=i.scores)+
  geom_point(aes(x=PC2,y=PC3,
                 color=treatment),size=2)+
  scale_color_viridis_d(option="B",end=.8)+
  facet_wrap(~sampling)

#start examining statistical relationship
trt<-as.factor(i.env$treatment)
seas<-i.env$season
samp<-i.env$sampling
yr<-as.factor(i.env$yr)
plts<-i.env$plot

tr.s.samp.mat3<-data.frame(model.matrix(~ trt*seas*samp+plts, 
             contrasts=list(trt="contr.helmert", seas="contr.helmert")))[,-1]
tr.s.samp.mat2<-data.frame(model.matrix(~ trt*seas+trt*samp+seas*samp+plts, 
                                       contrasts=list(trt="contr.helmert", seas="contr.helmert")))[,-1]

tr.s.yr.mat3<-data.frame(model.matrix(~ trt*seas*yr+plts, 
                            contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
tr.s.yr.mat2<-data.frame(model.matrix(~ trt*seas+trt*yr+yr*seas+plts, 
                                      contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
tr.s.yr.mat3.nop<-data.frame(model.matrix(~ trt*seas*yr, 
                                      contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]


i.rda.samp3<-rda(i.com.hel~.,data=tr.s.samp.mat3)
i.rda.samp2<-rda(i.com.hel~.,data=tr.s.samp.mat2)

i.rda.yr3<-rda(i.com.hel~.,data=tr.s.yr.mat3)#best model at the moment
i.rda.yr3nop<-rda(i.com.hel~.,data=tr.s.yr.mat3.nop)
i.rda.yr2<-rda(i.com.hel~.,data=tr.s.yr.mat2)

i.rda2<-rda(i.com.hel~1)
i.rda.samp
RsquareAdj(i.rda.samp3)
RsquareAdj(i.rda.samp2)
RsquareAdj(i.rda.yr3)
RsquareAdj(i.rda.yr3nop)
RsquareAdj(i.rda.yr2)

anova(i.rda.samp3,i.rda.samp2)
anova(i.rda.yr3,i.rda.yr3nop)
anova(i.rda.yr3,i.rda2)

plot(i.rda.yr3)
#bring in productivity data
sg<-sg_shoot%>%
  group_by(treatment,plot,sampling)%>%
  summarize(sg.sd=mean(SD))%>%
  mutate(sampling=case_when(
    sampling==1~0,
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17))

alg<-algae%>%
  group_by(treatment,plot,sampling)%>%
  summarize(alg=sum(abundance))%>%
  mutate(sampling=case_when(
    sampling==1~0,
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17))
i.env2<-left_join(i.env,sg)%>%
  left_join(alg)
#check to make sure data is still in the same order
summary(i.env[,1:4]==i.env2[,1:4])

sg.sd<-sg$sg.sd
alg.ab<-alg$alg
tr.s.yr.mat3.prod<-data.frame(model.matrix(~ trt*seas*yr+plts+sg.sd+alg.ab, 
                                      contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
tr.s.yr.mat3.sg<-data.frame(model.matrix(~ trt*seas*yr+plts+sg.sd, 
                                           contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
tr.s.yr.mat3.alg<-data.frame(model.matrix(~ trt*seas*yr+plts+alg.ab, 
                                           contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
i.rda.yr3.sg<-rda(i.com.hel~.,data=tr.s.yr.mat3.sg)
i.rda.yr3.alg<-rda(i.com.hel~.,data=tr.s.yr.mat3.alg)
i.rda.yr3.prod<-rda(i.com.hel~.,data=tr.s.yr.mat3.prod)
anova(i.rda.yr3.prod,i.rda.yr3)


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
