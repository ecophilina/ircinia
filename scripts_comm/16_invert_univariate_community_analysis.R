# inverts univariate analysis

library(tidyverse)
library(vegan)
library(lmerTest)
library(glmmTMB)

source("scripts_comm/02_community_data_org.R") 

# visualize data
ggplot(inv.uni)+
  geom_jitter(aes(x=sampling,y=spr,color=treatment,group=treatment))

ggplot(inv.uni)+
  geom_jitter(aes(x=sampling,y=j,color=treatment,group=treatment))

ggplot(inv.uni)+
  geom_jitter(aes(x=sampling,y=div,color=treatment,group=treatment))

# organize data to use offset

inv.uni.0<-inv.uni%>%
  filter(sampling==0)%>%
  select(treatment,plot,start.spr=spr,start.div=div,start.j=j)

inv.uni2<-inv.uni%>%
  filter(sampling!=0)%>%
  left_join(inv.uni.0)%>%
  mutate(change.spr=spr-start.spr,
         change.div=div-start.div,
         change.j=j-start.j)

inv.uni2$treatment<-factor(inv.uni2$treatment)
inv.uni$treatment<-factor(inv.uni$treatment)

hist(inv.uni2$change.spr)

summary(aov(spr~treatment,data=inv.uni%>%
      filter(sampling==0)))

ggplot(inv.uni%>%
  filter(sampling==0))+
  geom_point(aes(x=treatment,y=spr))

#models with treatment as the explanatory variable
#species richness

ispr.glmm<-glmmTMB::glmmTMB(change.spr~treatment*as.factor(sampling) + 
                              (1|plot),
                            #family= poisson,
                            data = inv.uni2%>%
                              filter(season=="summer")%>%
                              mutate(treatment=relevel(treatment, ref = "real")))

summary(ispr.glmm)

performance::r2(ispr.glmm)
# look at residuals
ispr.glmm_simres <- simulateResiduals(ispr.glmm)
testDispersion(ispr.glmm_simres)
plot(ispr.glmm_simres)

# evenness
ij.glmm<-glmmTMB::glmmTMB(change.j~treatment*as.factor(sampling) + 
                              (1|plot),
                            #family= poisson,
                            data = inv.uni2%>%
                              filter(season=="summer")%>%
                              mutate(treatment=relevel(treatment, ref = "real")))

summary(ij.glmm)

performance::r2(ij.glmm)
# look at residuals
ij.glmm_simres <- simulateResiduals(ij.glmm)
testDispersion(ij.glmm_simres)
plot(ij.glmm_simres)

# diversity
idiv.glmm<-glmmTMB::glmmTMB(change.div~treatment*as.factor(sampling) + 
                            (1|plot),
                          #family= poisson,
                          data = inv.uni2%>%
                            filter(season=="summer")%>%
                            mutate(treatment=relevel(treatment, ref = "real")))

summary(idiv.glmm)

performance::r2(idiv.glmm)
# look at residuals
idiv.glmm_simres <- simulateResiduals(idiv.glmm)
testDispersion(idiv.glmm_simres)
plot(idiv.glmm_simres)


# models with productivity as the explanatory variable
#species richness

ispr.glmm.grow<-glmmTMB::glmmTMB(change.spr~grow*as.factor(sampling) + sg.sd*as.factor(sampling) + 
                              (1|plot),
                            #family= poisson,
                            data = inv.uni2%>%
                              filter(season=="summer")%>%
                              mutate(treatment=relevel(treatment, ref = "real")))

summary(ispr.glmm.grow)

performance::r2(ispr.glmm)
performance::r2(ispr.glmm.grow)

# look at residuals
ispr.glmm_simres.grow <- simulateResiduals(ispr.glmm.grow)
testDispersion(ispr.glmm_simres.grow)
plot(ispr.glmm_simres.grow)

# evenness
ij.glmm<-glmmTMB::glmmTMB(change.j~treatment*as.factor(sampling) + 
                            (1|plot),
                          #family= poisson,
                          data = inv.uni2%>%
                            filter(season=="summer")%>%
                            mutate(treatment=relevel(treatment, ref = "real")))

summary(ij.glmm)

performance::r2(ij.glmm)
# look at residuals
ij.glmm_simres <- simulateResiduals(ij.glmm)
testDispersion(ij.glmm_simres)
plot(ij.glmm_simres)

# diversity
idiv.glmm<-glmmTMB::glmmTMB(change.div~treatment*as.factor(sampling) + 
                              (1|plot),
                            #family= poisson,
                            data = inv.uni2%>%
                              filter(season=="summer")%>%
                              mutate(treatment=relevel(treatment, ref = "real")))

summary(idiv.glmm)

performance::r2(idiv.glmm)
# look at residuals
idiv.glmm_simres <- simulateResiduals(idiv.glmm)
testDispersion(idiv.glmm_simres)
plot(idiv.glmm_simres)

# ispr.lmer<-lmer(change.spr~treatment*sampling + grow + 
#                  (1|plot),
#                data = inv.uni2%>%
#                  mutate(treatment=relevel(treatment, ref = "real")))
# 
# summary(ispr.lmer)
# plot(ispr.lmer)
# 
# ij.lmer<-lmer(change.j~treatment*sampling +  grow +  
#                (1|plot),
#              data = inv.uni2%>%
#                mutate(treatment=relevel(treatment, ref = "real")))
# 
# summary(ij.lmer)
# plot(ij.lmer)
# 
# idiv.lmer<-lmer(change.div~treatment*sampling + grow + 
#                  (1|plot),
#                data = inv.uni2%>%
#                  mutate(treatment=relevel(treatment, ref = "real")))
# 
# summary(idiv.lmer)
# plot(idiv.lmer)

# make dataset for plots
inv.pr<-ggeffects::ggpredict(ispr.lmer,terms=c("treatment","sampling"))%>%
  rename(treatment=x,sampling=group)%>%
  mutate(sampling=as.numeric(sampling),
         sampling=case_when(
           sampling==1~1,
           sampling==2~5,
           sampling==3~12,
           sampling==4~17))

inv.pr$treatment<-factor(inv.pr$treatment,levels=c("blank","fake","real"),labels=c("Control","Structure Control","Sponge"))

inv.plot<-inv.uni2%>%
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

inv.plot$treatment<-factor(inv.plot$treatment,levels=c("blank","fake","real"),labels=c("Control","Structure Control","Sponge"))


ggplot(inv.plot)+
  geom_hline(yintercept=0,linetype="dashed",alpha=.5)+
  geom_point(aes(x=sampling,y=mcspr,color=treatment),size=5,position=position_dodge(0.5))+
  geom_errorbar(aes(x=sampling,ymax=mcspr+spr95,ymin=mcspr-spr95,color=treatment),width=.4,position=position_dodge(0.5))+
  geom_line(data=inv.pr,aes(x=sampling,y=predicted,group=treatment,color=treatment))+
  geom_ribbon(data = inv.pr, 
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

ggsave("figures/inverts_species_richness.jpg",dpi=300)   


inv.j<-ggeffects::ggpredict(ij.lmer,terms=c("treatment","sampling"))%>%
  rename(treatment=x,sampling=group)%>%
  mutate(sampling=as.numeric(sampling),
         sampling=case_when(
           sampling==1~1,
           sampling==2~5,
           sampling==3~12,
           sampling==4~17))

inv.j$treatment<-factor(inv.j$treatment,levels=c("blank","fake","real"),labels=c("Control","Structure Control","Sponge"))

ggplot(inv.plot)+
  geom_hline(yintercept=0,linetype="dashed",alpha=.5)+
  geom_point(aes(x=sampling,y=mcj,color=treatment),size=5,position=position_dodge(0.5))+
  geom_errorbar(aes(x=sampling,ymax=mcj +j95,ymin=mcj-j95,color=treatment),width=.4,position=position_dodge(0.5))+
  geom_line(data=inv.j,aes(x=sampling,y=predicted,group=treatment,color=treatment))+
  geom_ribbon(data = inv.j, 
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

ggsave("figures/inverts_species_evenness.jpg",dpi=300)

inv.div<-ggeffects::ggpredict(idiv.lmer,terms=c("treatment","sampling"))%>%
  rename(treatment=x,sampling=group)%>%
  mutate(sampling=as.numeric(sampling),
         sampling=case_when(
           sampling==1~1,
           sampling==2~5,
           sampling==3~12,
           sampling==4~17))

inv.div$treatment<-factor(inv.div$treatment,levels=c("blank","fake","real"),labels=c("Control","Structure Control","Sponge"))

ggplot(inv.plot)+
  geom_hline(yintercept=0,linetype="dashed",alpha=.5)+
  geom_point(aes(x=sampling,y=mcdiv,color=treatment),size=5,position=position_dodge(0.5))+
  geom_errorbar(aes(x=sampling,ymax=mcdiv +div95,ymin=mcdiv-div95,color=treatment),width=.4,position=position_dodge(0.5))+
  geom_line(data=inv.div,aes(x=sampling,y=predicted,group=treatment,color=treatment))+
  geom_ribbon(data = inv.div, 
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

ggsave("figures/inverts_diversity.jpg",dpi=300)
