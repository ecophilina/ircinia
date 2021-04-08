# algae univariate analysis

library(tidyverse)
library(vegan)
library(lmerTest)
library(glmmTMB)
library(DHARMa)
library(performance)

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


# aspr.lmm<-lmer(spr~treatment*as.factor(sampling) + 
#     (1|plot),
#   data = alg.uni%>%
#     filter(season=="summer")%>%
#     mutate(treatment=relevel(treatment, ref = "real")))
# 
# summary(aspr.lmm)
# 
# # look at residuals
# spr.lmm_simres <- simulateResiduals(aspr.lmm)
# testDispersion(spr.lmm_simres)
# plot(spr.lmm_simres)


aspr.glmm<-glmmTMB::glmmTMB(change.spr~treatment*as.factor(sampling) + 
    (1|plot),
#   offset = start.spr,
  # family= tweedie,
  data = alg.uni2%>%
    filter(season=="summer")%>%
    mutate(treatment=relevel(treatment, ref = "real")))

summary(aspr.glmm)

performance::r2(aspr.glmm)
# look at residuals
spr.glmm_simres <- simulateResiduals(aspr.glmm)
testDispersion(spr.glmm_simres)
plot(spr.glmm_simres)


# aspr.glmm2<-glmmTMB::glmmTMB(spr~as.factor(sampling) * grow + 
#     as.factor(sampling) * sg.sd + 
#     sg.sd * grow + 
#     as.factor(sampling) * grow * sg.sd + 
#     (1|plot),
#   # family= tweedie,
#   data = alg.uni%>%
#     filter(season=="summer")%>%
#     mutate(treatment=relevel(treatment, ref = "real")))
# 
# summary(aspr.glmm2)
# performance::r2(aspr.glmm2)
# 
# # look at residuals
# spr.glmm_simres2 <- simulateResiduals(aspr.glmm2)
# testDispersion(spr.glmm_simres2)
# plot(spr.glmm_simres2)



alg.uni12 <- alg.uni %>% filter(sampling==12)
hist(alg.uni12$spr)

alg.uni0 <- alg.uni %>% filter(sampling==0)
hist(alg.uni0$spr)


aj.glmm<-glmmTMB::glmmTMB(change.j~treatment*as.factor(sampling) + 
    (1|plot),
 #   offset = start.j,
  data = alg.uni2%>%
    filter(season=="summer")%>%
    mutate(treatment=relevel(treatment, ref = "real")))

summary(aj.glmm)

performance::r2(aj.glmm)

# look at residuals
j.glmm_simres <- simulateResiduals(aj.glmm)
testDispersion(j.glmm_simres)
plot(j.glmm_simres)


adiv.glmm<-glmmTMB::glmmTMB(change.div~treatment*as.factor(sampling) + 
    (1|plot),
 #   offset = start.div,
  data = alg.uni2%>%
    filter(season=="summer")%>%
    mutate(treatment=relevel(treatment, ref = "real")))

summary(adiv.glmm)

performance::r2(adiv.glmm)
# look at residuals
div.glmm_simres <- simulateResiduals(adiv.glmm)
testDispersion(div.glmm_simres)
plot(div.glmm_simres)


# make dataset for plots
alg.pr<-ggeffects::ggpredict(aspr.lmer,terms=c("treatment","sampling"))%>%
  rename(treatment=x,sampling=group)%>%
  mutate(sampling=as.numeric(sampling),
         sampling=case_when(
           sampling==1~1,
           sampling==2~5,
           sampling==3~12,
           sampling==4~17))

alg.pr$treatment<-factor(alg.pr$treatment,levels=c("blank","fake","real"),labels=c("Control","Structure Control","Sponge"))

alg.plot<-alg.uni2%>%
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

alg.plot$treatment<-factor(alg.plot$treatment,levels=c("blank","fake","real"),labels=c("Control","Structure Control","Sponge"))


ggplot(alg.plot)+
  geom_hline(yintercept=0,linetype="dashed",alpha=.5)+
  geom_point(aes(x=sampling,y=mcspr,color=treatment),size=5,position=position_dodge(0.5))+
  geom_errorbar(aes(x=sampling,ymax=mcspr+spr95,ymin=mcspr-spr95,color=treatment),width=.4,position=position_dodge(0.5))+
  geom_line(data=alg.pr,aes(x=sampling,y=predicted,group=treatment,color=treatment))+
  geom_ribbon(data = alg.pr, 
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

ggsave("figures/algae_species_richness.jpg",dpi=300)   


alg.j<-ggeffects::ggpredict(aj.lmer,terms=c("treatment","sampling"))%>%
  rename(treatment=x,sampling=group)%>%
  mutate(sampling=as.numeric(sampling),
         sampling=case_when(
           sampling==1~1,
           sampling==2~5,
           sampling==3~12,
           sampling==4~17))

alg.j$treatment<-factor(alg.j$treatment,levels=c("blank","fake","real"),labels=c("Control","Structure Control","Sponge"))

ggplot(alg.plot)+
  geom_hline(yintercept=0,linetype="dashed",alpha=.5)+
  geom_point(aes(x=sampling,y=mcj,color=treatment),size=5,position=position_dodge(0.5))+
  geom_errorbar(aes(x=sampling,ymax=mcj +j95,ymin=mcj-j95,color=treatment),width=.4,position=position_dodge(0.5))+
  geom_line(data=alg.j,aes(x=sampling,y=predicted,group=treatment,color=treatment))+
  geom_ribbon(data = alg.j, 
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

ggsave("figures/algae_species_evenness.jpg",dpi=300)

alg.div<-ggeffects::ggpredict(adiv.lmer,terms=c("treatment","sampling"))%>%
  rename(treatment=x,sampling=group)%>%
  mutate(sampling=as.numeric(sampling),
         sampling=case_when(
           sampling==1~1,
           sampling==2~5,
           sampling==3~12,
           sampling==4~17))

alg.div$treatment<-factor(alg.div$treatment,levels=c("blank","fake","real"),labels=c("Control","Structure Control","Sponge"))

ggplot(alg.plot)+
  geom_hline(yintercept=0,linetype="dashed",alpha=.5)+
  geom_point(aes(x=sampling,y=mcdiv,color=treatment),size=5,position=position_dodge(0.5))+
  geom_errorbar(aes(x=sampling,ymax=mcdiv +div95,ymin=mcdiv-div95,color=treatment),width=.4,position=position_dodge(0.5))+
  geom_line(data=alg.div,aes(x=sampling,y=predicted,group=treatment,color=treatment))+
  geom_ribbon(data = alg.div, 
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

ggsave("figures/algae_diversity.jpg",dpi=300)