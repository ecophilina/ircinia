# this script will explore the invert data

library(tidyverse)
if(!require(vegan))install.packages("vegan"); library(vegan)

# bring in data
source("scripts/03_reimport.R")

i2<-inverts%>%
  pivot_wider(names_from=taxa,values_from=abundance,values_fill=0)

i2<-i2%>%
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


# NMDS
i.env<-i2 %>% select(treatment, plot, yr, sampling, season)
i.com<-i2 %>% select(-treatment, -plot, -yr, -sampling, -season)

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
plts<-as.factor(i.env$plot)

tr.s.samp.mat3<-data.frame(model.matrix(~ trt*seas*samp+plts, 
  contrasts=list(trt="contr.helmert", seas="contr.helmert")))[,-1]
tr.s.samp.mat2<-data.frame(model.matrix(~ trt*seas+trt*samp+seas*samp+plts, 
  contrasts=list(trt="contr.helmert", seas="contr.helmert")))[,-1]

tr.s.yr.mat3<-data.frame(model.matrix(~ trt*seas*yr+plts, 
  contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
tr.s.yr.mat2<-data.frame(model.matrix(~ trt*seas+trt*yr+yr*seas+plts, 
  contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
tr.s.yr.mat3.nop<-data.frame(model.matrix(~ trt*seas*yr, # without plot
  contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]


i.rda.samp3<-rda(i.com.hel~.,data=tr.s.samp.mat3)
i.rda.samp2<-rda(i.com.hel~.,data=tr.s.samp.mat2)

i.rda.yr3<-rda(i.com.hel~.,data=tr.s.yr.mat3)#best model at the moment
i.rda.yr3nop<-rda(i.com.hel~.,data=tr.s.yr.mat3.nop) # without plot
i.rda.yr2<-rda(i.com.hel~.,data=tr.s.yr.mat2)

i.rda2<-rda(i.com.hel~1)

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
i.env1<-left_join(i.env,sg)%>%
  left_join(alg)

#check to make sure data is still in the same order
summary(i.env[,1:4]==i.env1[,1:4])

# I'm worried that the helmert contrasts don't actually test difference between fake and real
# so attempt this part without blank data so that effect of cage applies to all plots
i.env2a <- i.env1 %>% filter(treatment!="blank")

# bind back the community data for just the retained plots
i3 <- left_join(i.env2a, i2)

# separate env and com again to sure order matches in this new dataset 
i.env2 <- i3 %>% select(treatment, plot, yr, sampling, season, sg.sd, alg)
i.com2 <- i3 %>% select(-treatment, -plot, -yr, -sampling, -season, -sg.sd, -alg)


i.com2$dummy<-1
com.dist2<-vegdist(i.com2,"bray")
i.mds2<-metaMDS(com.dist2,trymax = 100)
# plot(i.mds2)

# rda
i.com.hel2<-decostand(i.com2,"hellinger")


sg.sd<-i.env2$sg.sd
alg.ab<-i.env2$alg

trt<-as.factor(as.character(i.env2$treatment))
seas<-i.env2$season
samp<-i.env2$sampling
yr<-as.factor(i.env2$yr)
plts<-i.env2$plot

# rerun treatment model to make sure data matches exactly
tr.s.yr.mat3<-data.frame(model.matrix(~ trt*seas*yr+plts, 
  contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
i.rda.yr3<-rda(i.com.hel2~.,data=tr.s.yr.mat3)

# just add linear effect of primary producers

tr.s.yr.mat3.prod<-data.frame(model.matrix(~ trt*seas*yr+plts+sg.sd+alg.ab, 
  contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
# tr.s.yr.mat3.sg<-data.frame(model.matrix(~ trt*seas*yr+plts+sg.sd, 
#   contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
# tr.s.yr.mat3.alg<-data.frame(model.matrix(~ trt*seas*yr+plts+alg.ab, 
#   contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]

i.rda.yr3.prod<-rda(i.com.hel2~.,data=tr.s.yr.mat3.prod)
anova(i.rda.yr3.prod,i.rda.yr3)

# if above was sig we would check each type of producer separately
# i.rda.yr3.sg<-rda(i.com.hel~.,data=tr.s.yr.mat3.sg)
# i.rda.yr3.alg<-rda(i.com.hel~.,data=tr.s.yr.mat3.alg)

# add seasonal interaction with total of each type of primary producer
tr.s.yr.mat3.sg2<-data.frame(model.matrix(~ trt*seas*yr+plts+sg.sd*seas, 
  contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
tr.s.yr.mat3.alg2<-data.frame(model.matrix(~ trt*seas*yr+plts+alg.ab*seas, 
  contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]
# or interaction between types
tr.s.yr.mat3.prod2<-data.frame(model.matrix(~ trt*seas*yr+plts+sg.sd*alg.ab, 
  contrasts=list(trt="contr.helmert", seas="contr.helmert", yr="contr.helmert")))[,-1]

i.rda.yr3.sg2<-rda(i.com.hel2~.,data=tr.s.yr.mat3.sg2)
i.rda.yr3.alg2<-rda(i.com.hel2~.,data=tr.s.yr.mat3.alg2)
i.rda.prod.full<-rda(i.com.hel2~.,data=tr.s.yr.mat3.prod2)

anova(i.rda.yr3.sg2,i.rda.yr3)
anova(i.rda.yr3.alg2,i.rda.yr3)
anova(i.rda.prod.full,i.rda.yr3)
# doesn't help

# replace treatment with total of each type of primary producer
tr.s.yr.prod3<-data.frame(model.matrix(~ seas*yr*sg.sd+seas*yr*alg.ab+plts, 
  contrasts=list(seas="contr.helmert", yr="contr.helmert")))[,-1]

tr.s.yr.sg3<-data.frame(model.matrix(~ seas*yr*sg.sd+plts, 
  contrasts=list(seas="contr.helmert", yr="contr.helmert")))[,-1]
tr.s.yr.alg3<-data.frame(model.matrix(~ seas*yr*alg.ab+plts, 
  contrasts=list(seas="contr.helmert", yr="contr.helmert")))[,-1]

i.rda.prod3<-rda(i.com.hel2~.,data=tr.s.yr.prod3)
i.rda.sg3<-rda(i.com.hel2~.,data=tr.s.yr.sg3)
i.rda.alg3<-rda(i.com.hel2~.,data=tr.s.yr.alg3)

anova(i.rda.prod3,i.rda.sg3) 
anova(i.rda.prod3,i.rda.alg3) # algae is almost identical to full


tr.s.yr.prod2<-data.frame(model.matrix(~ seas*yr+sg.sd*seas+alg.ab*seas+plts, 
  contrasts=list(seas="contr.helmert", yr="contr.helmert")))[,-1]
i.rda.prod2<-rda(i.com.hel2~.,data=tr.s.yr.prod2)
anova(i.rda.prod3,i.rda.prod2)

tr.s.yr.alg2<-data.frame(model.matrix(~ seas*yr+alg.ab*seas+plts, 
  contrasts=list(seas="contr.helmert", yr="contr.helmert")))[,-1]
i.rda.alg2<-rda(i.com.hel2~.,data=tr.s.yr.alg2)
anova(i.rda.alg2,i.rda.prod2)

## any chance types interact with eachother? YES
tr.s.yr.prod2b<-data.frame(model.matrix(~ yr*seas+plts+sg.sd*alg.ab, 
  contrasts=list(seas="contr.helmert", yr="contr.helmert")))[,-1]
tr.s.yr.prod1<-data.frame(model.matrix(~ yr*seas+plts+sg.sd+alg.ab, 
  contrasts=list(seas="contr.helmert", yr="contr.helmert")))[,-1]

i.rda.prod2b<-rda(i.com.hel2~.,data=tr.s.yr.prod2b)
i.rda.prod1<-rda(i.com.hel2~.,data=tr.s.yr.prod1)
anova(i.rda.prod1,i.rda.prod2b)

## doesn't improve with season 3-way
# tr.s.yr.prod3b<-data.frame(model.matrix(~ yr*seas+plts+sg.sd*alg.ab*seas, 
#   contrasts=list(seas="contr.helmert", yr="contr.helmert")))[,-1]
# i.rda.prod3b<-rda(i.com.hel2~.,data=tr.s.yr.prod3b)
# anova(i.rda.prod3b,i.rda.prod2b)


# not nested so can't use anova, but R squared less than for treatment model for each alone
RsquareAdj(i.rda.yr3) #treatment model
RsquareAdj(i.rda.sg3) #seagrass model
RsquareAdj(i.rda.alg3) #algae model
RsquareAdj(i.rda.prod2b) #both producers interacting model comes closest but still pretty far off



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
