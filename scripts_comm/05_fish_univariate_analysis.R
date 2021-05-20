# fish univariate analysis
library(vegan)
library(lmerTest)
library(glmmTMB)
library(DHARMa)
library(pacman)
pacman::p_load(s20x, lme4, AICcmodavg, MASS)
library(tidyverse)

source("scripts_comm/02_community_data_org.R")

# visualize data
ggplot(fish.uni) +
  geom_jitter(aes(
    x = sampling,
    y = spr,
    color = treatment,
    group = treatment))

ggplot(fish.uni) +
  geom_jitter(aes(
    x = sampling,
    y = f.abund,
    color = treatment,
    group = treatment))

ggplot(fish.uni) +
  geom_jitter(aes(
    x = sampling,
    y = j,
    color = treatment,
    group = treatment))

ggplot(fish.uni) +
  geom_jitter(aes(
    x = sampling,
    y = div,
    color = treatment,
    group = treatment))

# organize data to use offset

fish.uni.0 <- fish.uni %>%
  filter(sampling == 0) %>%
  dplyr::select(
    treatment,
    plot,
    start.spr = spr,
    start.div = div,
    start.j = j,
    start.a = f.abund)

fish.uni2 <- fish.uni %>%
  filter(sampling != 0) %>%
  left_join(fish.uni.0) %>%
  mutate(
    change.spr = spr - start.spr,
    change.div = div - start.div,
    change.j = j - start.j,
    change.a = f.abund - start.a)

fish.uni2$treatment <- factor(fish.uni2$treatment)

# look at distribution of potential response variables
hist(fish.uni2$change.spr)
hist(fish.uni2$change.a)
hist(fish.uni2$change.div)
hist(fish.uni2$change.j)

# make sure there aren't differences between treatments originally

#Species Richness
summary(aov(spr ~ treatment, data = fish.uni %>%
              filter(sampling == 0)))

# No difference

# Abundance

summary(aov(f.abund ~ treatment, data = fish.uni %>%
              filter(sampling == 0)))

# no difference

# diversity

summary(aov(div ~ treatment, data = fish.uni %>%
              filter(sampling == 0)))

# no difference

# evenness

summary(aov(j ~ treatment, data = fish.uni %>%
              filter(sampling == 0)))
# No difference

## now start with species richness and do model selection

#Treatment only model
fspr.treat <- glmmTMB(change.spr ~ treatment * as.factor(sampling) +
                        (1 | plot),
                      data = fish.uni2 %>%
                        filter(season == "summer") %>%
                        mutate(treatment = relevel(treatment, ref = "real")))

# seagrass productivity model for richness
fspr.prod <- glmmTMB(
  change.spr ~ sg.prod.c * as.factor(sampling) +
    (1 | plot),
  #family= poisson,
  data = fish.uni2 %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

# seagrass structure model for richness
fspr.struct <- glmmTMB(
  change.spr ~ sg.sd.c * as.factor(sampling) +
    (1 | plot),
  #family= poisson,
  data = fish.uni2 %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

# algal abundance model for richness
fspr.alg <- glmmTMB(
  change.spr ~ a.abund * as.factor(sampling) +
    (1 | plot),
  #family= poisson,
  data = fish.uni2 %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment and seagrass productivity model
fspr.treat.prod <- glmmTMB(change.spr ~ treatment * as.factor(sampling) +
                             sg.prod.c * as.factor(sampling) + (1 | plot),
                           data = fish.uni2 %>%
                             filter(season == "summer") %>%
                             mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment and struct for richness
fspr.treat.struct <- glmmTMB(change.spr ~ treatment * as.factor(sampling) +
                               sg.sd.c * as.factor(sampling) + (1 | plot),
                             data = fish.uni2 %>%
                               filter(season == "summer") %>%
                               mutate(treatment = relevel(treatment, ref = "real")))
# combine treatment and algal abundance for richness
fspr.treat.alg <- glmmTMB(change.spr ~ treatment * as.factor(sampling) +
                               a.abund * as.factor(sampling) + (1 | plot),
                             data = fish.uni2 %>%
                               filter(season == "summer") %>%
                               mutate(treatment = relevel(treatment, ref = "real")))
# combine productivity and algal abundance for richness
fspr.prod.alg <- glmmTMB(change.spr ~ a.abund +
                               sg.sd.c * as.factor(sampling) + (1 | plot),
                             data = fish.uni2 %>%
                               filter(season == "summer") %>%
                               mutate(treatment = relevel(treatment, ref = "real")))

# combine productivity and struct for richness
fspr.prod.struct <- glmmTMB(change.spr ~ sg.prod.c * as.factor(sampling)  +
                              sg.sd.c * as.factor(sampling) + (1 | plot),
                            data = fish.uni2 %>%
                              filter(season == "summer") %>%
                              mutate(treatment = relevel(treatment, ref = "real")))

# combine structure and algal abundance for richness
# full model
fspr.full <- glmmTMB(change.spr ~ treatment * as.factor(sampling) +
                       sg.prod.c * as.factor(sampling)  +
                       sg.sd.c * as.factor(sampling) + (1 | plot),
                     data = fish.uni2 %>%
                       filter(season == "summer") %>%
                       mutate(treatment = relevel(treatment, ref = "real")))

# check residuals
glmm.resids(fspr.treat)
glmm.resids(fspr.prod)
glmm.resids(fspr.struct)
glmm.resids(fspr.treat.prod)
glmm.resids(fspr.treat.struct)
glmm.resids(fspr.prod.struct)
glmm.resids(fspr.full)


# all look good. Now do model selection

# Model selection
spr.cand.mod.names <- c("fspr.treat", "fspr.prod", "fspr.treat.prod", "fspr.struct",
                        "fspr.treat.struct","fspr.prod.struct","fspr.full")
spr.cand.mods <- list( ) 

# This function fills the list by model names
for(i in 1:length(spr.cand.mod.names)) {
  spr.cand.mods[[i]] <- get(spr.cand.mod.names[i]) }

# Function aictab does the AICc-based model comparison
print(aictab(cand.set = spr.cand.mods, 
             modnames = spr.cand.mod.names))

# treatment is the best model, followed by treatment + structure

# now do model selection for abundance

#Treatment only model
fa.treat <- glmmTMB(change.a ~ treatment * as.factor(sampling) +
                      (1 | plot),
                    #family= poisson,
                    data = fish.uni2 %>%
                      filter(season == "summer") %>%
                      mutate(treatment = relevel(treatment, ref = "real")))

# seagrass productivity model 
fa.prod <- glmmTMB(change.a ~ sg.prod.c * as.factor(sampling) +
                     (1 | plot),
                   #family= poisson,
                   data = fish.uni2 %>%
                     filter(season == "summer") %>%
                     mutate(treatment = relevel(treatment, ref = "real")))

# seagrass structure model 
fa.struct <- glmmTMB(change.a ~ sg.sd.c * as.factor(sampling) +
                       (1 | plot),
                     #family= poisson,
                     data = fish.uni2 %>%
                       filter(season == "summer") %>%
                       mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment and seagrass productivity model
fa.treat.prod <- glmmTMB(change.a ~ treatment * as.factor(sampling) +
                           sg.prod.c * as.factor(sampling) + (1 | plot),
                         data = fish.uni2 %>%
                           filter(season == "summer") %>%
                           mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment and struct for richness
fa.treat.struct <- glmmTMB(change.a ~ treatment * as.factor(sampling) +
                             sg.sd.c * as.factor(sampling) + (1 | plot),
                           data = fish.uni2 %>%
                             filter(season == "summer") %>%
                             mutate(treatment = relevel(treatment, ref = "real")))

# combine productivity and struct for richness
fa.prod.struct <- glmmTMB(change.a ~ sg.prod.c * as.factor(sampling)  +
                            sg.sd.c * as.factor(sampling) + (1 | plot),
                          data = fish.uni2 %>%
                            filter(season == "summer") %>%
                            mutate(treatment = relevel(treatment, ref = "real")))
#model doesn't converge

# full model
fa.full <- glmmTMB(change.a ~ treatment * as.factor(sampling) +
                     sg.prod.c * as.factor(sampling)  +
                     sg.sd.c * as.factor(sampling) + (1 | plot),
                   data = fish.uni2 %>%
                     filter(season == "summer") %>%
                     mutate(treatment = relevel(treatment, ref = "real")))

# check residuals
glmm.resids(fa.treat)
glmm.resids(fa.prod)
glmm.resids(fa.struct)
glmm.resids(fa.treat.prod)
glmm.resids(fa.treat.struct)
glmm.resids(fa.prod.struct)
glmm.resids(fa.full)


# Residuals are ok. Now do model selection

# Model selection
abund.cand.mod.names <- c("fa.treat", "fa.prod", "fa.treat.prod", "fa.struct",
                          "fa.treat.struct","fa.full")
abund.cand.mods <- list( ) 

# This function fills the list by model names
for(i in 1:length(abund.cand.mod.names)) {
  abund.cand.mods[[i]] <- get(abund.cand.mod.names[i]) }

# Function aictab does the AICc-based model comparison
print(aictab(cand.set = abund.cand.mods, 
             modnames = abund.cand.mod.names))

# treatment is the best model - followed by treatment + productivity