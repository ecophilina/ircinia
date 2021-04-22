# inverts univariate analysis

library(tidyverse)
library(vegan)
library(lmerTest)
library(glmmTMB)
library(DHARMa)
library(pacman)
pacman::p_load(s20x, lme4, AICcmodavg, MASS)

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

(i.sponge.model<-performance::r2(idiv.glmm))
# look at residuals
idiv.glmm_simres <- simulateResiduals(idiv.glmm)
testDispersion(idiv.glmm_simres)
plot(idiv.glmm_simres)


# models with productivity as the explanatory variable

# seagrass productivity model
idiv.prod.glmm<-glmmTMB::glmmTMB(change.div~sg.prod.c*as.factor(sampling) + 
                              (1|plot),
                            #family= poisson,
                            data = inv.uni2%>%
                              filter(season=="summer")%>%
                              mutate(treatment=relevel(treatment, ref = "real")))

summary(idiv.prod.glmm)

(i.prod.model<-performance::r2(idiv.prod.glmm))

# look at residuals
idiv.prod.glmm_simres <- simulateResiduals(idiv.prod.glmm)
testDispersion(idiv.prod.glmm_simres)
plot(idiv.prod.glmm_simres)

# # primary producer structure model
# idiv.struct.glmm<-glmmTMB::glmmTMB(change.div~pp.struct.c*as.factor(sampling) + 
#                                    (1|plot),
#                                  #family= poisson,
#                                  data = inv.uni2%>%
#                                    filter(season=="summer")%>%
#                                    mutate(treatment=relevel(treatment, ref = "real")))
# 
# summary(idiv.struct.glmm)
# 
# (i.struct.model<-performance::r2(idiv.struct.glmm))
# 
# # look at residuals
# idiv.struct.glmm_simres <- simulateResiduals(idiv.struct.glmm)
# testDispersion(idiv.struct.glmm_simres)
# plot(idiv.struct.glmm_simres)

# seagrass structure model
idiv.sg.struct.glmm<-glmmTMB::glmmTMB(change.div~sg.sd.c*as.factor(sampling) + 
                                     (1|plot),
                                   #family= poisson,
                                   data = inv.uni2%>%
                                     filter(season=="summer")%>%
                                     mutate(treatment=relevel(treatment, ref = "real")))

summary(idiv.sg.struct.glmm)

(i.sg.struct.model<-performance::r2(idiv.sg.struct.glmm))

# look at residuals
idiv.sg.struct.glmm_simres <- simulateResiduals(idiv.sg.struct.glmm)
testDispersion(idiv.sg.struct.glmm_simres)
plot(idiv.sg.struct.glmm_simres)

# # algae structure model
# idiv.alg.struct.glmm<-glmmTMB::glmmTMB(change.div~abund.c*as.factor(sampling) + 
#                                         (1|plot),
#                                       #family= poisson,
#                                       data = inv.uni2%>%
#                                         filter(season=="summer")%>%
#                                         mutate(treatment=relevel(treatment, ref = "real")))
# 
# summary(idiv.alg.struct.glmm)
# 
# (i.alg.struct.model<-performance::r2(idiv.alg.struct.glmm))
# 
# # look at residuals
# idiv.alg.struct.glmm_simres <- simulateResiduals(idiv.alg.struct.glmm)
# testDispersion(idiv.alg.struct.glmm_simres)
# plot(idiv.alg.struct.glmm_simres)


# trying an overall model with treatment and seagrass density


i.global.model<-glmmTMB::glmmTMB(change.div~treatment*as.factor(sampling)+
                                 sg.sd.c*as.factor(sampling) + 
                                 (1|plot),
                               #family= poisson,
                               data = inv.uni2%>%
                                 filter(season=="summer")%>%
                                 mutate(treatment=relevel(treatment, ref = "real")))


summary(i.global.model)

(i.global.model.r2<-performance::r2(i.global.model))

# look at residuals
i.global.model_simres <- simulateResiduals(i.global.model)
testDispersion(i.global.model_simres)
plot(i.global.model_simres)

# look at all R2
i.global.model.r2
i.sponge.model
i.sg.struct.model
i.prod.model

# look at aics
MuMIn::AICc(idiv.sg.struct.glmm)
MuMIn::AICc(i.global.model)
MuMIn::AICc(idiv.glmm)
MuMIn::AICc(idiv.prod.glmm)


# Model selection
cand.mod.names <- c("idiv.sg.struct.glmm","i.global.model","idiv.glmm","idiv.prod.glmm")
cand.mods <- list( ) 

# This function fills the list by model names
for(i in 1:length(cand.mod.names)) {
  cand.mods[[i]] <- get(cand.mod.names[i]) }

# Function aictab does the AICc-based model comparison
print(aictab(cand.set = cand.mods, 
             modnames = cand.mod.names))


# What if we start with breaking downt the effect of treatments on richness, eveness and diversity?

(i.sponge.model<-performance::r2(idiv.glmm))
(i.sponge.model<-performance::r2(ispr.glmm))
(i.sponge.model<-performance::r2(ij.glmm))

# since richness is much more related to treament than diversity, maybe the test of the imporance of prod and structure should be asked for richness instead (or in addition to for diversity)? 

# seagrass productivity model for richness
ispr.prod.glmm<-glmmTMB::glmmTMB(change.spr~sg.prod.c*as.factor(sampling) + 
    (1|plot),
  #family= poisson,
  data = inv.uni2%>%
    filter(season=="summer")%>%
    mutate(treatment=relevel(treatment, ref = "real")))

summary(ispr.prod.glmm)
(ispr.prod.model<-performance::r2(ispr.prod.glmm))
# ispr.prod.glmm_simres <- simulateResiduals(ispr.prod.glmm)
# testDispersion(ispr.prod.glmm_simres)
# plot(ispr.prod.glmm_simres)


# combine treatment and seagrass productivity model for richness
ispr.treat.prod.glmm<-glmmTMB::glmmTMB(change.spr~
    treatment*as.factor(sampling) + 
    sg.prod.c*as.factor(sampling) + 
    (1|plot),
  #family= poisson,
  data = inv.uni2%>%
    filter(season=="summer")%>%
    mutate(treatment=relevel(treatment, ref = "real")))

summary(ispr.treat.prod.glmm)
(ispr.prod.model<-performance::r2(ispr.treat.prod.glmm))
# ispr.treat.prod.glmm_simres <- simulateResiduals(ispr.treat.prod.glmm)
# testDispersion(ispr.treat.prod.glmm_simres)
# plot(ispr.treat.prod.glmm_simres)

# seagrass structure model for richness
ispr.sg.struct.glmm<-glmmTMB::glmmTMB(change.spr~sg.sd.c*as.factor(sampling) + 
    (1|plot),
  #family= poisson,
  data = inv.uni2%>%
    filter(season=="summer")%>%
    mutate(treatment=relevel(treatment, ref = "real")))

summary(ispr.sg.struct.glmm)
(ispr.sg.struct.model<-performance::r2(ispr.sg.struct.glmm))
# ispr.sg.struct.glmm_simres <- simulateResiduals(ispr.sg.struct.glmm)
# testDispersion(ispr.sg.struct.glmm_simres)
# plot(ispr.sg.struct.glmm_simres)

# combine treatment and struct for richness
ispr.sg.treat.struct.glmm<-glmmTMB::glmmTMB(change.spr~
    treatment*as.factor(sampling) + 
    sg.sd.c*as.factor(sampling) + 
    (1|plot),
  #family= poisson,
  data = inv.uni2%>%
    filter(season=="summer")%>%
    mutate(treatment=relevel(treatment, ref = "real")))

summary(ispr.sg.treat.struct.glmm)
(ispr.treat.struct.model<-performance::r2(ispr.sg.treat.struct.glmm))
# idiv.sg.treat.struct.glmm_simres <- simulateResiduals(ispr.sg.treat.struct.glmm)
# testDispersion(ispr.sg.treat.struct.glmm)
# plot(idiv.sg.treat.struct.glmm_simres)

(i.sponge.model<-performance::r2(ispr.glmm))
(ispr.prod.model<-performance::r2(ispr.prod.glmm))
(ispr.treat.prod.model<-performance::r2(ispr.treat.prod.glmm))
(ispr.sg.struct.model<-performance::r2(ispr.sg.struct.glmm))
(ispr.treat.struct.model<-performance::r2(ispr.sg.treat.struct.glmm))


# AIC(ispr.glmm, ispr.prod.glmm, ispr.treat.prod.glmm, ispr.sg.struct.glmm, ispr.sg.treat.struct.glmm)

MuMIn::AICc(ispr.glmm) 
MuMIn::AICc(ispr.prod.glmm)
MuMIn::AICc(ispr.treat.prod.glmm)
MuMIn::AICc(ispr.sg.struct.glmm)
MuMIn::AICc(ispr.sg.treat.struct.glmm)


# Model selection
cand.mod.names <- c("ispr.glmm", "ispr.prod.glmm", "ispr.treat.prod.glmm", "ispr.sg.struct.glmm", "ispr.sg.treat.struct.glmm")
cand.mods <- list( ) 

# This function fills the list by model names
for(i in 1:length(cand.mod.names)) {
  cand.mods[[i]] <- get(cand.mod.names[i]) }

# Function aictab does the AICc-based model comparison
print(aictab(cand.set = cand.mods, 
  modnames = cand.mod.names))




