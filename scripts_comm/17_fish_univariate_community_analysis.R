# fish univariate analysis

library(tidyverse)
library(vegan)
library(lmerTest)
library(glmmTMB)

source("scripts/11_community_analysis_start.R") 

# visualize data
ggplot(fish.uni)+
  geom_jitter(aes(x=sampling,y=spr,color=treatment,group=treatment))

ggplot(fish.uni)+
  geom_jitter(aes(x=sampling,y=j,color=treatment,group=treatment))

ggplot(fish.uni)+
  geom_jitter(aes(x=sampling,y=div,color=treatment,group=treatment))

# organize data to use offset

fish.uni.0<-fish.uni%>%
  filter(sampling==0)%>%
  select(treatment,plot,start.spr=spr,start.div=div,start.j=j)

fish.uni2<-fish.uni%>%
  filter(sampling!=0)%>%
  left_join(fish.uni.0)%>%
  mutate(change.spr=spr-start.spr,
         change.div=div-start.div,
         change.j=j-start.j)

fish.uni2$treatment<-factor(fish.uni2$treatment)
hist(fish.uni2$change.spr)

summary(aov(spr~treatment,data=fish.uni%>%
              filter(sampling==0)))

ggplot(fish.uni%>%
         filter(sampling==0))+
  geom_jitter(aes(x=treatment,y=spr))

#species richness

fspr.glmm<-glmmTMB::glmmTMB(change.spr~treatment*as.factor(sampling) + 
                              (1|plot),
                            #family= poisson,
                            data = fish.uni2%>%
                              filter(season=="summer")%>%
                              mutate(treatment=relevel(treatment, ref = "real")))

summary(fspr.glmm)

performance::r2(fspr.glmm)

# look at residuals
fspr.glmm_simres <- simulateResiduals(fspr.glmm)
testDispersion(fspr.glmm_simres)
plot(fspr.glmm_simres)

# evenness
fj.glmm<-glmmTMB::glmmTMB(change.j~treatment*as.factor(sampling) + 
                            (1|plot),
                          #family= poisson,
                          data = fish.uni2%>%
                            filter(season=="summer")%>%
                            mutate(treatment=relevel(treatment, ref = "real")))

summary(fj.glmm)

performance::r2(fj.glmm)

# look at residuals
fj.glmm_simres <- simulateResiduals(fj.glmm)
testDispersion(fj.glmm_simres)
plot(fj.glmm_simres)

# diversity
fdiv.glmm<-glmmTMB::glmmTMB(change.div~treatment*as.factor(sampling) + 
                              (1|plot),
                            #family= poisson,
                            data = fish.uni2%>%
                              filter(season=="summer")%>%
                              mutate(treatment=relevel(treatment, ref = "real")))

summary(fdiv.glmm)

performance::r2(fdiv.glmm)

# look at residuals
fdiv.glmm_simres <- simulateResiduals(fdiv.glmm)
testDispersion(fdiv.glmm_simres)
plot(fdiv.glmm_simres)

# fspr.lmer<-lmer(change.spr~treatment*sampling + grow + 
#                   (1|plot),
#                 data = fish.uni2%>%
#                   mutate(treatment=relevel(treatment, ref = "real")))
# 
# summary(fspr.lmer)
# plot(fspr.lmer)
# 
# fj.lmer<-lmer(change.j~treatment*sampling +  grow +  
#                 (1|plot),
#               data = fish.uni2%>%
#                 mutate(treatment=relevel(treatment, ref = "real")))
# 
# summary(fj.lmer)
# plot(fj.lmer)
# 
# fdiv.lmer<-lmer(change.div~treatment*sampling + grow + 
#                   (1|plot),
#                 data = fish.uni2%>%
#                   mutate(treatment=relevel(treatment, ref = "real")))
# 
# summary(fdiv.lmer)
# plot(fdiv.lmer)

# make dataset for plots
fish.pr<-ggeffects::ggpredict(fspr.lmer,terms=c("treatment","sampling"))%>%
  rename(treatment=x,sampling=group)%>%
  mutate(sampling=as.numeric(sampling),
         sampling=case_when(
           sampling==1~1,
           sampling==2~5,
           sampling==3~12,
           sampling==4~17))

fish.pr$treatment<-factor(fish.pr$treatment,levels=c("blank","fake","real"),labels=c("Control","Structure Control","Sponge"))

fish.plot<-fish.uni2%>%
  group_by(sampling,treatment)%>%
  summarize(mcspr=mean(change.spr),
            sdspr=sd(change.spr),
            sespr=sd(change.spr)/sqrt(5),
            spr95=sespr*1.96,
            mcdiv=mean(change.div),
            sddiv=sd(change.div),
            sediv=sddiv/sqrt(5),
            div95=sediv*1.96,
            mcj=mean(change.j),
            sdj=sd(change.j),
            sej=sdj/sqrt(5),
            j95=sej*1.96)

fish.plot$treatment<-factor(fish.plot$treatment,levels=c("blank","fake","real"),labels=c("Control","Structure Control","Sponge"))


ggplot(fish.plot)+
  geom_hline(yintercept=0,linetype="dashed",alpha=.5)+
  geom_point(aes(x=sampling,y=mcspr,color=treatment),size=5,position=position_dodge(0.5))+
  geom_errorbar(aes(x=sampling,ymax=mcspr+spr95,ymin=mcspr-spr95,color=treatment),width=.4,position=position_dodge(0.5))+
  geom_line(data=fish.pr,aes(x=sampling,y=predicted,group=treatment,color=treatment))+
  geom_ribbon(data = fish.pr, 
              aes(sampling,ymin = conf.low, ymax = conf.high,
                  group = treatment, fill = treatment), alpha = 0.2)+
  scale_color_viridis_d(name="", option="A",end=0.6)+
  scale_fill_viridis_d(name="", option="A",end=0.6)+
  theme_bw()+
  ylab("Change in Species Richness")+
  xlab("Months into the Experiment")+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size=14))

ggsave("figures/fish_species_richness.jpg",dpi=300)   


fish.j<-ggeffects::ggpredict(fj.lmer,terms=c("treatment","sampling"))%>%
  rename(treatment=x,sampling=group)%>%
  mutate(sampling=as.numeric(sampling),
         sampling=case_when(
           sampling==1~1,
           sampling==2~5,
           sampling==3~12,
           sampling==4~17))

fish.j$treatment<-factor(fish.j$treatment,levels=c("blank","fake","real"),labels=c("Control","Structure Control","Sponge"))

ggplot(fish.plot)+
  geom_hline(yintercept=0,linetype="dashed",alpha=.5)+
  geom_point(aes(x=sampling,y=mcj,color=treatment),size=5,position=position_dodge(0.5))+
  geom_errorbar(aes(x=sampling,ymax=mcj +j95,ymin=mcj-j95,color=treatment),width=.4,position=position_dodge(0.5))+
  geom_line(data=fish.j,aes(x=sampling,y=predicted,group=treatment,color=treatment))+
  geom_ribbon(data = fish.j, 
              aes(sampling,ymin = conf.low, ymax = conf.high,
                  group = treatment, fill = treatment), alpha = 0.2)+
  scale_color_viridis_d(name="", option="A",end=0.6)+
  scale_fill_viridis_d(name="", option="A",end=0.6)+
  theme_bw()+
  ylab("Change in Evenness")+
  xlab("Months into the Experiment")+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size=14))

ggsave("figures/fish_species_evenness.jpg",dpi=300)

fish.div<-ggeffects::ggpredict(fdiv.lmer,terms=c("treatment","sampling"))%>%
  rename(treatment=x,sampling=group)%>%
  mutate(sampling=as.numeric(sampling),
         sampling=case_when(
           sampling==1~1,
           sampling==2~5,
           sampling==3~12,
           sampling==4~17))

fish.div$treatment<-factor(fish.div$treatment,levels=c("blank","fake","real"),labels=c("Control","Structure Control","Sponge"))

ggplot(fish.plot)+
  geom_hline(yintercept=0,linetype="dashed",alpha=.5)+
  geom_point(aes(x=sampling,y=mcdiv,color=treatment),size=5,position=position_dodge(0.5))+
  geom_errorbar(aes(x=sampling,ymax=mcdiv +div95,ymin=mcdiv-div95,color=treatment),width=.4,position=position_dodge(0.5))+
  geom_line(data=fish.div,aes(x=sampling,y=predicted,group=treatment,color=treatment))+
  geom_ribbon(data = fish.div, 
              aes(sampling,ymin = conf.low, ymax = conf.high,
                  group = treatment, fill = treatment), alpha = 0.2)+
  scale_color_viridis_d(name="", option="A",end=0.6)+
  scale_fill_viridis_d(name="", option="A",end=0.6)+
  theme_bw()+
  ylab("Change in Diversity")+
  xlab("Months into the Experiment")+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size=14))

ggsave("figures/fish_diversity.jpg",dpi=300)
