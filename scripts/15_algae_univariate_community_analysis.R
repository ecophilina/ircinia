# algae univariate analysis

library(tidyverse)
library(vegan)
library(lmerTest)
library(glmmTMB)

source("scripts/11_community_analysis_start.R") 

# visualize data
ggplot(alg.uni)+
  geom_jitter(aes(x=sampling,y=spr,color=treatment,group=treatment))

ggplot(alg.uni)+
  geom_jitter(aes(x=sampling,y=j,color=treatment,group=treatment))

ggplot(alg.uni)+
  geom_jitter(aes(x=sampling,y=div,color=treatment,group=treatment))

# organize data to use offset

alg.uni.0<-alg.uni%>%
  filter(sampling==0)%>%
  select(treatment,plot,start.spr=spr,start.div=div,start.j=j)

alg.uni2<-alg.uni%>%
  filter(sampling!=0)%>%
  left_join(alg.uni.0)%>%
  mutate(change.spr=spr-start.spr,
         change.div=div-start.div,
         change.j=j-start.j)

alg.uni2$treatment<-factor(alg.uni2$treatment)


spr.lmer<-lmer(change.spr~treatment*sampling + grow + 
                 (1|plot),
               data = alg.uni2%>%
                 mutate(treatment=relevel(treatment, ref = "real")))

summary(spr.lmer)
plot(spr.lmer)

# make dataset for plots

alg.plot<-alg.uni2%>%
  group_by(sampling,treatment)%>%
  summarize(mcspr=mean(change.spr),
            sdspr=sd(change.spr),
            mcdiv=mean(change.div),
            sddiv=sd(change.div),
            mcj=mean(change.j),
            sdj=sd(change.j))

ggplot(alg.plot)+
  geom_hline(yintercept=0,linetype="dashed",alpha=.5)+
  geom_point(aes(x=sampling,y=mcspr,color=treatment),size=5,position=position_dodge(0.5))+
  geom_errorbar(aes(x=sampling,ymax=mcspr+sdspr,ymin=mcspr-sdspr,color=treatment),width=.4,position=position_dodge(0.5))
