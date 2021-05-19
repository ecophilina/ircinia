# inverts univariate analysis
library(vegan)
library(lmerTest)
library(glmmTMB)
library(DHARMa)
library(pacman)
pacman::p_load(s20x, lme4, AICcmodavg, MASS)
library(tidyverse)

source("scripts_comm/02_community_data_org.R")

# visualize data
ggplot(inv.uni) +
  geom_jitter(aes(
    x = sampling,
    y = spr,
    color = treatment,
    group = treatment))

ggplot(inv.uni) +
  geom_jitter(aes(
    x = sampling,
    y = i.abund,
    color = treatment,
    group = treatment))

ggplot(inv.uni) +
  geom_jitter(aes(
    x = sampling,
    y = j,
    color = treatment,
    group = treatment))

ggplot(inv.uni) +
  geom_jitter(aes(
    x = sampling,
    y = div,
    color = treatment,
    group = treatment))

# organize data to use offset

inv.uni.0 <- inv.uni %>%
  filter(sampling == 0) %>%
  dplyr::select(
    treatment,
    plot,
    start.spr = spr,
    start.div = div,
    start.j = j,
    start.a = i.abund)

inv.uni2 <- inv.uni %>%
  filter(sampling != 0) %>%
  left_join(inv.uni.0) %>%
  mutate(
    change.spr = spr - start.spr,
    change.div = div - start.div,
    change.j = j - start.j,
    change.a = i.abund - start.a)

inv.uni2$treatment <- factor(inv.uni2$treatment)

# look at distribution of potential response variables
hist(inv.uni2$change.spr)
hist(inv.uni2$change.a)
hist(inv.uni2$change.div)
hist(inv.uni2$change.j)

# make sure there aren't differences between treatments originally

#Species Richness
summary(aov(spr ~ treatment, data = inv.uni %>%
              filter(sampling == 0)))

# No difference

# Abundance

summary(aov(i.abund ~ treatment, data = inv.uni %>%
              filter(sampling == 0)))

# no difference

# diversity

summary(aov(div ~ treatment, data = inv.uni %>%
              filter(sampling == 0)))

# no difference

# evenness

summary(aov(j ~ treatment, data = inv.uni %>%
              filter(sampling == 0)))
# No difference

## now start with species richness and do model selection

#Treatment only model
ispr.treat <- glmmTMB(change.spr ~ treatment * as.factor(sampling) +
    (1 | plot),
  data = inv.uni2 %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

# seagrass productivity model for richness
ispr.prod <- glmmTMB(
  change.spr ~ sg.prod.c + as.factor(sampling) +
    (1 | plot),
  #family= poisson,
  data = inv.uni2 %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

# seagrass structure model for richness
ispr.struct <- glmmTMB(
  change.spr ~ sg.sd.c + as.factor(sampling) +
    (1 | plot),
  #family= poisson,
  data = inv.uni2 %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment and seagrass productivity model
ispr.treat.prod <- glmmTMB(change.spr ~ treatment * as.factor(sampling) +
    sg.prod.c +(1 | plot),
  data = inv.uni2 %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment and struct for richness
ispr.treat.struct <- glmmTMB(change.spr ~ treatment * as.factor(sampling) +
    sg.sd.c +  (1 | plot),
  data = inv.uni2 %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

# combine productivity and struct for richness
# Note that productivity and struct are correlated with each other so their 
# individual contributions can't be assessed in the full model, just their combined effect:

plot(sg.prod.c~sg.sd.c, data = inv.uni2 %>% filter(season == "summer") %>% mutate(treatment = relevel(treatment, ref = "real")) )


ispr.prod.struct <- glmmTMB(change.spr ~ sg.prod.c + as.factor(sampling)  +
                               sg.sd.c + (1 | plot),
                             data = inv.uni2 %>%
                               filter(season == "summer") %>%
                               mutate(treatment = relevel(treatment, ref = "real")))

# full model
ispr.full <- glmmTMB(change.spr ~ treatment * as.factor(sampling) +
                              sg.prod.c + sg.sd.c + as.factor(sampling) + (1 | plot),
                            data = inv.uni2 %>%
                              filter(season == "summer") %>%
                              mutate(treatment = relevel(treatment, ref = "real")))


# check residuals
glmm.resids(ispr.treat)
glmm.resids(ispr.prod)
glmm.resids(ispr.struct)
glmm.resids(ispr.treat.prod)
glmm.resids(ispr.treat.struct)
glmm.resids(ispr.prod.struct)
glmm.resids(ispr.full)


# all look good. Now do model selection

# Model selection
spr.cand.mod.names <- c("ispr.treat", "ispr.prod", "ispr.treat.prod", "ispr.struct",
                    "ispr.treat.struct","ispr.prod.struct","ispr.full")
spr.cand.mods <- list( ) 

# This function fills the list by model names
for(i in 1:length(spr.cand.mod.names)) {
  spr.cand.mods[[i]] <- get(spr.cand.mod.names[i]) }

# Function aictab does the AICc-based model comparison
print(aictab(cand.set = spr.cand.mods, 
             modnames = spr.cand.mod.names))

# results changed. Now the full model is the best model - 

#look at this model

summary(ispr.treat.prod)

# now do model selection for abundance

#Treatment only model
ia.treat <- glmmTMB(change.a ~ treatment * as.factor(sampling) + 
                        (1 | plot),
                      data = inv.uni2 %>%
                        filter(season == "summer") %>%
                        mutate(treatment = relevel(treatment, ref = "real")))

# seagrass productivity model 
ia.prod <- glmmTMB(change.a ~ sg.prod.c + as.factor(sampling) +
    (1 | plot),
  data = inv.uni2 %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

# seagrass structure model 
ia.struct <- glmmTMB(change.a ~ sg.sd.c + as.factor(sampling) +
    (1 | plot),
  #family= poisson,
  data = inv.uni2 %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment and seagrass productivity model
ia.treat.prod <- glmmTMB(change.a ~ treatment * as.factor(sampling) +
                             sg.prod.c + (1 | plot),
                           data = inv.uni2 %>%
                             filter(season == "summer") %>%
                             mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment and struct for abundance
ia.treat.struct <- glmmTMB(change.a ~ treatment * as.factor(sampling) +
                               sg.sd.c+  (1 | plot),
                             data = inv.uni2 %>%
                               filter(season == "summer") %>%
                               mutate(treatment = relevel(treatment, ref = "real")))

# combine productivity and struct for abundance
ia.prod.struct <- glmmTMB(change.a ~ sg.prod.c + as.factor(sampling)  +
                              sg.sd.c  + (1 | plot),
                            data = inv.uni2 %>%
                              filter(season == "summer") %>%
                              mutate(treatment = relevel(treatment, ref = "real")))

# full model
ia.full <- glmmTMB(change.a ~ treatment * as.factor(sampling) +
                       sg.prod.c  +
                       sg.sd.c  + (1 | plot),
                     data = inv.uni2 %>%
                       filter(season == "summer") %>%
                       mutate(treatment = relevel(treatment, ref = "real")))

# check residuals
glmm.resids(ia.treat)
glmm.resids(ia.prod) 
glmm.resids(ia.struct)
glmm.resids(ia.treat.prod)
glmm.resids(ia.treat.struct)
glmm.resids(ia.prod.struct) 
glmm.resids(ia.full)
# pretty terrible


# Try a count distribution and include the before condition instead of modeling differences

ia.treat.a <- glmmTMB(i.abund ~ treatment * as.factor(sampling) + 
    # offset(-log(start.a+1)) +
    (1 | plot),
  family =  poisson,
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

glmm.resids(ia.treat.a) #looks pretty good
summary(ia.treat.a)

ia.prod.a <- glmmTMB(i.abund ~  as.factor(sampling) + 
    sg.prod.c  +
    # offset(-log(start.a+1)) + 
    (1 | plot),
  family = poisson, #nbinom2,
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))
# convergence error: 
# In (function (start, objective, gradient = NULL, hessian = NULL,  :NA/NaN function evaluation)
# but gives result

ia.treat.prod.a <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
    sg.prod.c  +
    # offset(-log(start.a+1)) +
    (1 | plot),
  family =  poisson, #nbinom2,
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))
# convergence error: 
# In (function (start, objective, gradient = NULL, hessian = NULL,  :NA/NaN function evaluation)
# but gives result

# ia.treat.prod.a <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
#     sg.prod.c  +
#     # offset(-log(start.a+1)) +
#     (1 | plot),
#   family =  nbinom2,
#   data = inv.uni %>%
#     filter(season == "summer") %>%
#     mutate(treatment = relevel(treatment, ref = "real")))
# converges, but nbinoms don't work for just treatment modal

ia.struct.a <- glmmTMB(i.abund ~ as.factor(sampling) +
    sg.sd.c +     
    # offset(-log(start.a+1)) + 
    (1 | plot),
  family = poisson,
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

ia.treat.struct.a <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
    sg.sd.c +     
    # offset(-log(start.a+1)) + 
    (1 | plot),
  family = poisson,
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

ia.prod.struct.a <- glmmTMB(i.abund ~  as.factor(sampling) + 
    sg.prod.c +
    sg.sd.c  +     
    # offset(-log(start.a+1)) + 
    (1 | plot),
  family =  poisson, #nbinom2,
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))
# convergence error: 
# In (function (start, objective, gradient = NULL, hessian = NULL,  :NA/NaN function evaluation)
# but gives result

ia.full.a <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
    sg.prod.c +
    sg.sd.c +   
    # offset(-log(start.a+1)) + 
    (1 | plot),
  family = poisson, #nbinom2,
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))
# convergence error: 
# In (function (start, objective, gradient = NULL, hessian = NULL,  :NA/NaN function evaluation)
# but gives result

# ia.full.a <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
#     sg.prod.c  +
#     sg.sd.c  +
#     # offset(-log(start.a+1)) +
#     (1 | plot),
#   family = nbinom1,
#   data = inv.uni %>%
#     filter(season == "summer") %>%
#     mutate(treatment = relevel(treatment, ref = "real")))

summary(ia.full.a)

# compare residuals btw model configurations
glmm.resids(ia.treat)
glmm.resids(ia.treat.a)
glmm.resids(ia.prod)
glmm.resids(ia.prod.a)
glmm.resids(ia.struct)
glmm.resids(ia.struct.a)
glmm.resids(ia.treat.prod)
glmm.resids(ia.treat.prod.a)
glmm.resids(ia.treat.struct)
glmm.resids(ia.treat.struct.a)
glmm.resids(ia.prod.struct)
glmm.resids(ia.prod.struct.a)
glmm.resids(ia.full)
glmm.resids(ia.full.a)


# summary(ia.full.o)

# Residuals aren't a lot better. Now do model selection


# Model selection using linear model of difference from start
abund.cand.mod.names <- c("ia.treat", "ia.prod", "ia.treat.prod", "ia.struct",
  "ia.treat.struct","ia.full")
abund.cand.mods <- list( ) 

# This function fills the list by model names
for(i in 1:length(abund.cand.mod.names)) {
  abund.cand.mods[[i]] <- get(abund.cand.mod.names[i]) }

# Function aictab does the AICc-based model comparison
print(aictab(cand.set = abund.cand.mods, 
  modnames = abund.cand.mod.names))


# Model selection using poisson model of raw abunadance
abund.cand.mod.names <- c("ia.treat.a", "ia.prod.a", "ia.treat.prod.a", #
  "ia.struct.a",  "ia.full.a",
  "ia.prod.struct.a", "ia.treat.struct.a")
abund.cand.mods <- list( ) 
for(i in 1:length(abund.cand.mod.names)) {
  abund.cand.mods[[i]] <- get(abund.cand.mod.names[i]) }
print(aictab(cand.set = abund.cand.mods, 
  modnames = abund.cand.mod.names))

# treatment * productivity is best model in both cases but equal in support to the just treatment model. 
# the linear model suggests:
summary(ia.treat) 
# by first month all treatments differ and difference increases over time
summary(ia.treat.prod)
# most of change between 1 and 12 months can be explained by seagrass productivity


summary(ia.treat.a) 
# no diff at start, large increase in abundance in first month
# fake differs at 1 month, both differ by 12 months
summary(ia.treat.prod.a)
# increases in sponge plots not entirely correlated with seagrass productivity, but it does eliminate sig of differences btw treatments 