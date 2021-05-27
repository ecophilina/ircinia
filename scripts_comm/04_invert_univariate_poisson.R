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


## species richness model selection ####

#Treatment only model
ispr.treat <- glmmTMB(spr ~ treatment * as.factor(sampling) +
    (1 | plot),
  family = compois(link = "log"),
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

# seagrass productivity model for richness
ispr.prod <- glmmTMB(spr ~ as.factor(sampling) + sg.prod.c + 
    (1 | plot),
  family = compois(link = "log"),
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

# seagrass structure model for richness
ispr.struct <- glmmTMB(spr ~ as.factor(sampling) + sg.sd.c + 
    (1 | plot),
  family = compois(link = "log"),
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

#algal abundance model for richness
ispr.alg <- glmmTMB(spr ~ as.factor(sampling) + a.abund.c + 
    (1 | plot),
  family = compois(link = "log"),
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment and seagrass productivity model
ispr.treat.prod <- glmmTMB(spr ~ treatment * as.factor(sampling) +
    sg.prod.c +(1 | plot),
  family = compois(link = "log"),
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment and struct for richness
ispr.treat.struct <- glmmTMB(spr ~ treatment * as.factor(sampling) +
    sg.sd.c + (1 | plot),
  family = compois(link = "log"),
    data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment and algal abundance for richness
ispr.treat.alg <- glmmTMB(spr ~ treatment * as.factor(sampling) +
    a.abund.c + (1 | plot),
  family = compois(link = "log"),
    data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

#combine productivity and algal abundance for richness
ispr.prod.alg <- glmmTMB(spr ~ as.factor(sampling)  +
    a.abund.c + sg.prod.c + (1 | plot),
  family = compois(link = "log"),
    data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))

# combine seagrass productivity and struct for richness
# Note that productivity and struct are correlated with each other so their 
# individual contributions can't be assessed in the full model, just their combined effect:

# plot(sg.prod.c~sg.sd.c, data = inv.uni %>% filter(season == "summer") %>% mutate(treatment = relevel(treatment, ref = "real")) )


ispr.prod.struct <- glmmTMB(spr ~ as.factor(sampling)  +
      sg.sd.c + sg.prod.c + (1 | plot),
  family = compois(link = "log"),
      data = inv.uni %>%
      filter(season == "summer") %>%
      mutate(treatment = relevel(treatment, ref = "real")))

#combine structure and algal abundance for richness
ispr.struct.alg <- glmmTMB(spr ~ as.factor(sampling)  +
      sg.sd.c + a.abund.c + (1 | plot),
  family = compois(link = "log"),
      data = inv.uni %>%
      filter(season == "summer") %>%
      mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment, productivity, and algal abundance for richenss
ispr.treat.prod.alg <- glmmTMB(spr ~ treatment * as.factor(sampling)  +
      sg.prod.c + a.abund.c + (1 | plot),
  family = compois(link = "log"),
      data = inv.uni %>%
      filter(season == "summer") %>%
      mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment, productivity, and seagrass structure for richness
ispr.treat.prod.struct <- glmmTMB(spr ~ treatment * as.factor(sampling) +
      sg.prod.c + sg.sd.c + (1 | plot),
  family = compois(link = "log"),
      data = inv.uni %>%
      filter(season == "summer") %>%
      mutate(treatment = relevel(treatment, ref = "real")))

#combine treatment, seagrass structure, and algal abundance for richness
ispr.treat.struct.alg <- glmmTMB(spr ~ treatment * as.factor(sampling) +
      sg.sd.c + a.abund.c + (1| plot),
  family = compois(link = "log"),
      data = inv.uni %>%
      filter(season == "summer") %>%
      mutate(treatment = relevel(treatment, ref = "real")))

# full model for richness (treatment, seagrass structure, seagrass productivity, and algal abundance)
ispr.full <- glmmTMB(spr ~ treatment * as.factor(sampling) +
      sg.sd.c + sg.prod.c + a.abund.c + (1|plot),
  family = compois(link = "log"),
      data = inv.uni %>%
      filter(season == "summer") %>%
      mutate(treatment = relevel(treatment, ref = "real")))

# check residuals
glmm.resids(ispr.treat)
glmm.resids(ispr.prod)
glmm.resids(ispr.struct)
glmm.resids(ispr.alg) 
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

# #algal abundance model without time variable
# ispr.alg.notime <- glmmTMB(spr ~ a.abund.c + 
#     (1 | plot),
#   family = poisson,
#   data = inv.uni %>%
#     filter(season == "summer") %>%
#     mutate(treatment = relevel(treatment, ref = "real")))
# 
# glmm.resids(ispr.alg.notime)
# 

# Model selection
spr.cand.mod.names <- c(
"ispr.treat",
"ispr.prod",
"ispr.struct",
"ispr.alg",
# "ispr.alg.notime", 
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


# top models of invert spr ####
summary(ispr.treat)
# top model is algae
summary(ispr.alg)


# treatment is delta 4.42 AIC worse, this appears to be because sampling*treatment is highly correlated with algae abundance but takes up way more degrees of freedom! 

ggplot(data = inv.uni %>% filter(season == "summer") %>% 
    mutate(treatment = factor(treatment, 
      levels = c("blank", "fake", "real"), 
      labels = c("Control", "Structure Control", "Sponge")))) + 
  geom_jitter(aes(a.abund, spr, colour = treatment, alpha = as.factor(sampling)), 
    width = 0.2, height = 0.2, size = 2.5) + 
  # geom_line(aes(a.abund, spr, group = as.factor(plot))) +
  geom_smooth(method = "lm", aes(a.abund, spr), colour = "black") + 
  scale_alpha_discrete(range = c(0.3,1), name = "Months") + 
  scale_color_viridis_d(option="A", begin=0, end=0.6,"")+
  xlab("Algae abundance") +
  ylab("Invertebrate species richness") +
  # coord_cartesian(expand = F) +
  ggsidekick::theme_sleek()
ggsave("invert-spr-by-algae.png", width = 5, height = 3)


# treatment model for comparison
summary(ispr.treat)
# confirm that the control plots don't change sig with time
ispr.treat.f <- glmmTMB(spr ~ treatment *  as.factor(sampling) +
    (1 | plot),
  family = poisson,
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "fake")))

ispr.treat.b <- glmmTMB(spr ~ treatment * as.factor(sampling) +
    (1 | plot),
  family = poisson,
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "blank")))

summary(ispr.treat.f)
summary(ispr.treat.b) #blank at 1 month probably increased by chance, diff gone by 12 months



# model selection for abundance ####


# Treatment only model
ia.treat <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
  (1 | plot),
family = poisson,
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# seagrass productivity model
ia.prod <- glmmTMB(i.abund ~ as.factor(sampling) +
  sg.prod.c + (1 | plot),
family = poisson,
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# seagrass structure model
ia.struct <- glmmTMB(i.abund ~ as.factor(sampling) +
  sg.sd.c + (1 | plot),
family = poisson,
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# algal abundance model for abundance
ia.alg <- glmmTMB(i.abund ~ a.abund.c + as.factor(sampling) +
  (1 | plot),
family = poisson,
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)
# combine treatment and seagrass productivity model
ia.treat.prod <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
  sg.prod.c + (1 | plot),
family = poisson,
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment and struct for abundance
ia.treat.struct <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
  sg.sd.c + (1 | plot),
family = poisson,
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment and algal abundance for abundance
ia.treat.alg <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
  a.abund.c + (1 | plot),
family = poisson,
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine seagrass productivity and algal abundance for abundance
ia.prod.alg <- glmmTMB(i.abund ~ as.factor(sampling) +
  sg.prod.c + a.abund.c + (1 | plot),
family = poisson,
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine productivity and struct for abundance
ia.prod.struct <- glmmTMB(i.abund ~ as.factor(sampling) +
  sg.prod.c + sg.sd.c + (1 | plot),
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine structure and algal abundance for abundance
ia.struct.alg <- glmmTMB(i.abund ~ as.factor(sampling) +
  sg.sd.c + a.abund.c + (1 | plot),
family = poisson,
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment, productivity, and algal abundance for abundance
ia.treat.prod.alg <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
  sg.prod.c + a.abund.c + (1 | plot),
family = poisson,
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment, productivity and structure for abundance
ia.treat.prod.struct <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
  sg.prod.c + sg.sd.c + (1 | plot),
family = poisson,
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment, structure, and algal abundance for abundance
ia.treat.struct.alg <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
  a.abund.c + sg.sd.c + (1 | plot),
family = poisson,
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)
# full model for abundance
ia.full <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
  sg.prod.c + sg.sd.c + a.abund.c + (1 | plot),
family = poisson,
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)
# check residuals
glmm.resids(ia.treat)
glmm.resids(ia.prod)
glmm.resids(ia.struct)
glmm.resids(ia.alg) # ugly residuals
glmm.resids(ia.treat.prod)
glmm.resids(ia.treat.struct)
glmm.resids(ia.treat.alg)
glmm.resids(ia.prod.alg)
glmm.resids(ia.prod.struct)# odd resids
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

# confirm that the control plots don't change sig with time
ia.treat.prod.f <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
    sg.prod.c + (1 | plot),
  family = poisson,
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "fake"))
)
ia.treat.prod.b <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
    sg.prod.c + (1 | plot),
  family = poisson,
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "blank"))
)

summary(ia.treat.prod.f)
summary(ia.treat.prod.b)
