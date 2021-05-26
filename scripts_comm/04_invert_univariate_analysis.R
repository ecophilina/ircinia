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

#algal abundance model for richness
ispr.alg <- glmmTMB(
  change.spr ~ a.abund.c + as.factor(sampling) +
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

# combine treatment and algal abundance for richness
ispr.treat.alg <- glmmTMB(change.spr ~ treatment * as.factor(sampling) +
    a.abund.c + (1 | plot),
    data = inv.uni2 %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

#combine productivity and algal abundance for richness
ispr.prod.alg <- glmmTMB(change.spr ~ sg.prod.c + as.factor(sampling)  +
    a.abund.c + (1 | plot),
    data = inv.uni2 %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

# combine seagrass productivity and struct for richness
# Note that productivity and struct are correlated with each other so their 
# individual contributions can't be assessed in the full model, just their combined effect:

plot(sg.prod.c ~ sg.sd.c, data = inv.uni2 %>% filter(season == "summer") %>% 
    mutate(treatment = relevel(treatment, ref = "real")) )


ispr.prod.struct <- glmmTMB(change.spr ~ sg.prod.c + as.factor(sampling)  +
      sg.sd.c + (1 | plot),
      data = inv.uni2 %>%
      filter(season == "summer") %>%
      mutate(treatment = relevel(treatment, ref = "real")))

#combine structure and algal abundance for richness
ispr.struct.alg <- glmmTMB(change.spr ~ a.abund.c + as.factor(sampling)  +
      sg.sd.c + (1 | plot),
      data = inv.uni2 %>%
      filter(season == "summer") %>%
      mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment, productivity, and algal abundance for richenss
ispr.treat.prod.alg <- glmmTMB(change.spr ~ treatment * as.factor(sampling)  +
      sg.prod.c + a.abund.c + (1 | plot),
      data = inv.uni2 %>%
      filter(season == "summer") %>%
      mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment, productivity, and seagrass structure for richness
ispr.treat.prod.struct <- glmmTMB(change.spr ~ treatment * as.factor(sampling) +
      sg.prod.c + sg.sd.c + (1 | plot),
      data = inv.uni2 %>%
      filter(season == "summer") %>%
      mutate(treatment = relevel(treatment, ref = "real")))

#combine treatment, seagrass structure, and algal abundance for richness
ispr.treat.struct.alg <- glmmTMB(change.spr ~ treatment * as.factor(sampling) +
      sg.sd.c + a.abund.c+ (1| plot),
      data = inv.uni2 %>%
      filter(season == "summer") %>%
      mutate(treatment = relevel(treatment, ref = "real")))

# full model for richness (treatment, seagrass structure, seagrass productivity, and algal abundance)
ispr.full <- glmmTMB(change.spr~ treatment * as.factor(sampling) +
      sg.sd.c + sg.prod.c + a.abund.c+ (1|plot),
      data = inv.uni2 %>%
      filter(season == "summer") %>%
      mutate(treatment = relevel(treatment, ref = "real")))

# check residuals
glmm.resids(ispr.treat)
glmm.resids(ispr.prod)
glmm.resids(ispr.struct)
glmm.resids(ispr.alg) # messy
glmm.resids(ispr.treat.prod)
glmm.resids(ispr.treat.struct)
glmm.resids(ispr.treat.alg)
glmm.resids(ispr.prod.alg)
glmm.resids(ispr.prod.struct)
glmm.resids(ispr.struct.alg)
glmm.resids(ispr.treat.prod.alg)
glmm.resids(ispr.treat.prod.struct)
glmm.resids(ispr.treat.struct.alg)
glmm.resids(ispr.full)


# ispr.alg and ispr.prod.struct are kind of iffy. Now do model selection

# Model selection
spr.cand.mod.names <- c(
"ispr.treat",
"ispr.prod",
"ispr.struct",
"ispr.alg",
"ispr.treat.prod",
"ispr.treat.struct",
"ispr.treat.alg",
"ispr.prod.alg",
"ispr.prod.struct",
"ispr.struct.alg",
"ispr.treat.prod.alg",
"ispr.treat.prod.struct",
"ispr.treat.struct.alg",
"ispr.full")
spr.cand.mods <- list( ) 

# This function fills the list by model names
for(i in 1:length(spr.cand.mod.names)) {
  spr.cand.mods[[i]] <- get(spr.cand.mod.names[i]) }

# Function aictab does the AICc-based model comparison
print(aictab(cand.set = spr.cand.mods, 
             modnames = spr.cand.mod.names))

# it looks like treatment and productivity is now the best model for species richness.....- 

#look at top model
summary(ispr.treat.prod)


# model selection for abundance ####

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

#algal abundance model for abundance
ia.alg <- glmmTMB(change.a ~ a.abund.c+ as.factor(sampling) +
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

#combine treatment and algal abundance for abundance
ia.treat.alg <- glmmTMB(change.a ~ treatment * as.factor(sampling) +
                             a.abund.c+  (1 | plot),
                           data = inv.uni2 %>%
                             filter(season == "summer") %>%
                             mutate(treatment = relevel(treatment, ref = "real")))

#combine seagrass productivity and algal abundance for abundance
ia.prod.alg <- glmmTMB(change.a ~ sg.prod.c + as.factor(sampling)  +
                          a.abund.c + (1 | plot),
                          data = inv.uni2 %>%
                            filter(season == "summer") %>%
                            mutate(treatment = relevel(treatment, ref = "real")))

# combine productivity and struct for abundance
ia.prod.struct <- glmmTMB(change.a ~ sg.prod.c + as.factor(sampling)  +
                              sg.sd.c  + (1 | plot),
                            data = inv.uni2 %>%
                              filter(season == "summer") %>%
                              mutate(treatment = relevel(treatment, ref = "real")))

#combine structure and algal abundance for abundance
ia.struct.alg <- glmmTMB(change.a ~ sg.sd.c + as.factor(sampling)  +
                    a.abund.c+ (1 | plot),
                          data = inv.uni2 %>%
                            filter(season == "summer") %>%
                            mutate(treatment = relevel(treatment, ref = "real")))

#combine treatment, productivity, and algal abundance for abundance
ia.treat.prod.alg <- glmmTMB(change.a ~ treatment * as.factor(sampling) +
                     sg.prod.c  +
                     a.abund.c  + (1 | plot),
                   data = inv.uni2 %>%
                     filter(season == "summer") %>%
                     mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment, productivity and structure for abundance
ia.treat.prod.struct <- glmmTMB(change.a ~ treatment * as.factor(sampling) +
                       sg.prod.c  +
                       sg.sd.c  + (1 | plot),
                     data = inv.uni2 %>%
                       filter(season == "summer") %>%
                       mutate(treatment = relevel(treatment, ref = "real")))

#combine treatment, structure, and algal abundance for abundance
ia.treat.struct.alg <- glmmTMB(change.a ~ treatment * as.factor(sampling) +
                     a.abund.c  +
                     sg.sd.c  + (1 | plot),
                   data = inv.uni2 %>%
                     filter(season == "summer") %>%
                     mutate(treatment = relevel(treatment, ref = "real")))
#full model for abundance
ia.full <- glmmTMB(change.a ~ treatment * as.factor(sampling) +
                     sg.prod.c  +
                     sg.sd.c  + a.abund.c+ (1 | plot),
                   data = inv.uni2 %>%
                     filter(season == "summer") %>%
                     mutate(treatment = relevel(treatment, ref = "real")))

# check residuals
glmm.resids(ia.treat)
glmm.resids(ia.prod)
glmm.resids(ia.struct)
glmm.resids(ia.alg) # ugly residuals
glmm.resids(ia.treat.prod)
glmm.resids(ia.treat.struct)
glmm.resids(ia.treat.alg)
glmm.resids(ia.prod.alg)
glmm.resids(ia.prod.struct)
glmm.resids(ia.struct.alg)
glmm.resids(ia.treat.prod.alg)
glmm.resids(ia.treat.prod.struct)
glmm.resids(ia.treat.struct.alg)
glmm.resids(ia.full)


# ia.alg has ugly residuals. Now do model selection

# Model selection
abund.cand.mod.names <- c(
  "ia.treat",
  "ia.prod",
  "ia.struct",
  "ia.alg",
  "ia.treat.prod",
  "ia.treat.struct",
  "ia.treat.alg",
  "ia.prod.alg",
  "ia.prod.struct",
  "ia.struct.alg",
  "ia.treat.prod.alg",
  "ia.treat.prod.struct",
  "ia.treat.struct.alg",
  "ia.full")
abund.cand.mods <- list( ) 

# This function fills the list by model names
for(i in 1:length(abund.cand.mod.names)) {
  abund.cand.mods[[i]] <- get(abund.cand.mod.names[i]) }

# Function aictab does the AICc-based model comparison
print(aictab(cand.set = abund.cand.mods, 
             modnames = abund.cand.mod.names))

# we have three models at the top - treat * sg prod, treat, treat * sg struct

summary(ia.treat.prod)

summary(ia.treat)

summary(ia.treat.struct)
