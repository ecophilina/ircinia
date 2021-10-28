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
    group = treatment
  ))

ggplot(fish.uni) +
  geom_jitter(aes(
    x = sampling,
    y = f.abund,
    color = treatment,
    group = treatment
  ))

ggplot(fish.uni) +
  geom_jitter(aes(
    x = sampling,
    y = j,
    color = treatment,
    group = treatment
  ))

ggplot(fish.uni) +
  geom_jitter(aes(
    x = sampling,
    y = div,
    color = treatment,
    group = treatment
  ))

# organize data to use offset

fish.uni.0 <- fish.uni %>%
  filter(sampling == 0) %>%
  dplyr::select(
    treatment,
    plot,
    start.spr = spr,
    start.div = div,
    start.j = j,
    start.a = f.abund
  )

fish.uni2 <- fish.uni %>%
  filter(sampling != 0) %>%
  left_join(fish.uni.0) %>%
  mutate(
    change.spr = spr - start.spr,
    change.div = div - start.div,
    change.j = j - start.j,
    change.a = f.abund - start.a
  )

fish.uni2$treatment <- factor(fish.uni2$treatment)
fish.uni$treatment <- factor(fish.uni$treatment)
# look at distribution of potential response variables
hist(fish.uni2$change.spr)
hist(fish.uni2$change.a)
hist(fish.uni2$change.div)
hist(fish.uni2$change.j)

# make sure there aren't differences between treatments originally

# Species Richness
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

# Treatment only model
fspr.treat <- glmmTMB(spr ~ treatment * as.factor(sampling) +
  (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# seagrass productivity model for richness
fspr.prod <- glmmTMB(
  spr ~ as.factor(sampling) + sg.prod.c +  (1 | plot),
  family = compois(link = "log"),
  data = fish.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real"))
)

# seagrass structure model for richness
fspr.struct <- glmmTMB(
  spr ~ as.factor(sampling) + sg.sd.c +  (1 | plot),
  family = compois(link = "log"),
  data = fish.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real"))
)

# algal abundance model for richness
fspr.alg <- glmmTMB(
  spr ~ as.factor(sampling) + a.abund.c + (1 | plot),
  family = compois(link = "log"),
  data = fish.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment and seagrass productivity model
fspr.treat.prod <- glmmTMB(spr ~ treatment * as.factor(sampling) +
  sg.prod.c + (1 | plot),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment and struct for richness
fspr.treat.struct <- glmmTMB(spr ~ treatment * as.factor(sampling) +
  sg.sd.c + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)
# combine treatment and algal abundance for richness
fspr.treat.alg <- glmmTMB(spr ~ treatment * as.factor(sampling) +
  a.abund.c + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)
# combine productivity and algal abundance for richness
fspr.prod.alg <- glmmTMB(spr ~ as.factor(sampling) +
    a.abund.c + sg.prod.c + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine productivity and struct for richness
fspr.prod.struct <- glmmTMB(spr ~ as.factor(sampling) + 
    sg.prod.c + sg.sd.c + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine structure and algal abundance for richness
fspr.struct.alg <- glmmTMB(spr ~ as.factor(sampling) + 
  a.abund.c + sg.sd.c + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment, productivity, and algal abundance for richness
fspr.treat.prod.alg <- glmmTMB(spr ~ treatment * as.factor(sampling) +
  sg.prod.c + a.abund.c + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment, productivity, and structure for richness
fspr.treat.prod.struct <- glmmTMB(spr ~ treatment * as.factor(sampling) +
  sg.prod.c + sg.sd.c + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment, structure, and algal abundance for richness
fspr.treat.struct.alg <- glmmTMB(spr ~ treatment * as.factor(sampling) +
  sg.sd.c + a.abund.c + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# full model for richness
fspr.full <- glmmTMB(spr ~ treatment * as.factor(sampling) +
  sg.sd.c + a.abund.c + sg.prod.c + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# check residuals
glmm.resids(fspr.treat) 
glmm.resids(fspr.prod)
glmm.resids(fspr.struct) 
glmm.resids(fspr.alg)
glmm.resids(fspr.treat.prod)
glmm.resids(fspr.treat.struct)
glmm.resids(fspr.treat.alg)
glmm.resids(fspr.prod.alg)
glmm.resids(fspr.prod.struct)
glmm.resids(fspr.struct.alg)
glmm.resids(fspr.treat.prod.alg)
glmm.resids(fspr.treat.prod.struct) 
glmm.resids(fspr.treat.struct.alg)
glmm.resids(fspr.full) 

# not quite as good as change version. Now do model selection

# Model selection
spr.cand.mod.names <- c(
  "fspr.treat",
  "fspr.prod",
  "fspr.struct",
  "fspr.alg",
  "fspr.treat.prod",
  "fspr.treat.struct",
  "fspr.treat.alg",
  "fspr.prod.alg",
  "fspr.prod.struct",
  "fspr.struct.alg",
  "fspr.treat.prod.alg",
  "fspr.treat.prod.struct",
  "fspr.treat.struct.alg",
  "fspr.full"
)
spr.cand.mods <- list()

# This function fills the list by model names
for (i in 1:length(spr.cand.mod.names)) {
  spr.cand.mods[[i]] <- get(spr.cand.mod.names[i])
}

# Function aictab does the AICc-based model comparison
(fspr.aic<-data.frame(aictab(
  cand.set = spr.cand.mods,
  modnames = spr.cand.mod.names
)))
write_rds(fspr.aic,"working_data/FishSprAIC.rds")

# top model: I'm a bit confused as to which groups differ when, but only real interacts with sampling 
#to reference in manuscript
(fspr.treat.sum<-summary(fspr.treat))
write_rds(fspr.treat.sum,"working_data/FishSprTreatSum.rds")

# confirm that the control plots don't change sig with time
fspr.treat.f <- glmmTMB(spr ~ treatment * as.factor(sampling) +
    (1 | plot),
  family = compois(link = "log"),
  data = fish.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "fake"))
)
fspr.treat.b <- glmmTMB(spr ~ treatment * as.factor(sampling) +
    (1 | plot),
  family = compois(link = "log"),
  data = fish.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "blank"))
)
summary(fspr.treat.f)
summary(fspr.treat.b)

# algae isn't sig. so we can conclude that it's inclusion is spurious
summary(fspr.treat.alg)



## equally supported model with non-sig effect of algae; it's pretty messy
f.alg.pr<-ggeffects::ggpredict(fspr.treat.alg,terms=c("a.abund.c", "sampling", "treatment"))%>%
  filter(group == 12) %>% # plotting prediction for 12 months in
  mutate(treatment = factor(facet,
    levels = c("blank", "fake", "real"),
    labels = c("Control", "Structure", "Sponge")))
f.alg.pr$x <- (f.alg.pr$x * fish.uni$abund.global.sd[1]*2 )+ fish.uni$a.abund.global[1]

ggplot(data = fish.uni %>% filter(season == "summer") %>% 
    mutate(treatment = factor(treatment, 
      levels = c("blank", "fake", "real"), 
      labels = c("Control", "Structure", "Sponge")))) + 
  geom_ribbon(data = f.alg.pr,
    aes(x, ymin = conf.low, ymax = conf.high
      , fill = treatment
      ), alpha = 0.2)+
  geom_jitter(aes(a.abund, spr, colour = treatment, alpha = as.factor(sampling)), 
    width = 0.3, height = 0.1, size = 2.5) + 
  geom_line(data=f.alg.pr,aes(x=x,y=predicted, colour = treatment))+
  scale_alpha_discrete(range = c(0.3,1), name = "Months") + 
  scale_color_viridis_d(option="A", begin=0, end=0.65,"")+
  scale_fill_viridis_d(option="A", begin=0, end=0.65,"")+
  xlab("Algae Abundance") +
  ylab("Fish Species Richness") +
  # coord_cartesian(expand = F, xlim = c(-0.1,35.35)) +
  coord_cartesian(expand = F, ylim = c(-0.1,7.5), xlim = c(-0.1,35.35)) +
  ggsidekick::theme_sleek()

ggsave("figures/fish_spr_by_algae.jpg", width = 5, height = 3)
# 
### same thing for other well supported model with sg density
f.struct.pr<-ggeffects::ggpredict(fspr.treat.struct,terms=c("sg.sd.c", "sampling", "treatment"))%>%
  filter(group == 12) %>% # plotting prediction for 12 months in
  mutate(treatment = factor(facet,
    levels = c("blank", "fake", "real"),
    labels = c("Control", "Structure", "Sponge")))
f.struct.pr$x <- (f.struct.pr$x * fish.uni$sg.sd.global.sd[1]*2 )+ fish.uni$sg.sd.global[1]

ggplot(data = fish.uni %>% filter(season == "summer") %>% 
    mutate(treatment = factor(treatment, 
      levels = c("blank", "fake", "real"), 
      labels = c("Control", "Structure", "Sponge")))) + 
  geom_ribbon(data = f.struct.pr,
    aes(x, ymin = conf.low, ymax = conf.high
      , fill = treatment
    ), alpha = 0.2)+
  geom_jitter(aes(sg.sd, spr, colour = treatment, alpha = as.factor(sampling)), 
    width = 0.3, height = 0.1, size = 2.5) + 
  geom_line(data=f.struct.pr,aes(x=x,y=predicted, colour = treatment))+
  scale_alpha_discrete(range = c(0.3,1), name = "Months") + 
  scale_color_viridis_d(option="A", begin=0, end=0.65,"")+
  scale_fill_viridis_d(option="A", begin=0, end=0.65,"")+
  xlab("Seagrass Density") +
  ylab("Fish Species Richness") +
  coord_cartesian(expand = F) +
  # coord_cartesian(expand = F, ylim = c(-0.1,7.5), xlim = c(-0.1,35.35)) +
  ggsidekick::theme_sleek()

ggsave("figures/fish_spr_by_sgsd.jpg", width = 5, height = 3)


## probably not interesting but for comparison with invert one and less busy
alg.pr<-ggeffects::ggpredict(fspr.alg,terms=c("a.abund.c", "sampling"))%>%
  filter(group == 12)  # plotting prediction for 12 months in
alg.pr$x <- (alg.pr$x * fish.uni$abund.global.sd[1]*2 )+ fish.uni$a.abund.global[1]

ggplot(data = fish.uni %>% filter(season == "summer") %>% 
    mutate(treatment = factor(treatment, 
      levels = c("blank", "fake", "real"), 
      labels = c("Control", "Structure", "Sponge")))) + 
  geom_ribbon(data = alg.pr,
    aes(x, ymin = conf.low, ymax = conf.high
      # , fill = treatment
    ), alpha = 0.2)+
  geom_jitter(aes(a.abund, spr, colour = treatment, alpha = as.factor(sampling)), 
    width = 0.3, height = 0.1, size = 2.5) + 
  # geom_line(data=alg.pr,aes(x=x,y=predicted, colour = treatment))+
  geom_line(data=alg.pr,aes(x=x,y=predicted))+
  scale_alpha_discrete(range = c(0.3,1), name = "Months") + 
  scale_color_viridis_d(option="A", begin=0, end=0.65,"")+
  scale_fill_viridis_d(option="A", begin=0, end=0.65,"")+
  xlab("Algae Abundance") +
  ylab("Fish Species Richness") +
  # coord_cartesian(expand = F, xlim = c(-0.1,35.35)) +
  coord_cartesian(expand = F, ylim = c(-0.1,7.5), xlim = c(-0.1,35.35)) +
  ggsidekick::theme_sleek()

ggsave("figures/fish_spr_by_algae2.jpg", width = 5, height = 3)


# now do model selection for abundance

# Treatment only model
fa.treat <- glmmTMB(f.abund ~ treatment * as.factor(sampling) +
  (1 | plot),
family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# seagrass productivity model
fa.prod <- glmmTMB(f.abund ~ sg.prod.c + as.factor(sampling) +
  (1 | plot),
family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# seagrass structure model
fa.struct <- glmmTMB(f.abund ~ sg.sd.c + as.factor(sampling) +
  (1 | plot),
family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# algal abundance for abundance model
fa.alg <- glmmTMB(f.abund ~ a.abund.c + as.factor(sampling) +
  (1 | plot),
family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)
# combine treatment and seagrass productivity model
fa.treat.prod <- glmmTMB(f.abund ~ treatment * as.factor(sampling) +
  sg.prod.c + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment and struct for abundance
fa.treat.struct <- glmmTMB(f.abund ~ treatment * as.factor(sampling) +
  sg.sd.c + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment and algal abundance for abundance
fa.treat.alg <- glmmTMB(f.abund ~ treatment * as.factor(sampling) +
  a.abund.c + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine productivity and algal abundance for abundance
fa.prod.alg <- glmmTMB(f.abund ~ sg.prod.c + as.factor(sampling) +
  a.abund.c + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine productivity and struct for abundance
fa.prod.struct <- glmmTMB(f.abund ~ sg.prod.c + sg.sd.c +
  as.factor(sampling) + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine algal abundance and structure for abundance
fa.struct.alg <- glmmTMB(f.abund ~ sg.sd.c + a.abund.c +
  as.factor(sampling) + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment, productivity, and algal abundance for abundance
fa.treat.prod.alg <- glmmTMB(f.abund ~ treatment * as.factor(sampling) +
  sg.prod.c + a.abund.c + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)
# combine treatment, productivity, and structure for abundance
fa.treat.prod.struct <- glmmTMB(f.abund ~ treatment * as.factor(sampling) +
  sg.prod.c + sg.sd.c + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment, structure, and algal abundance for abundance
fa.treat.struct.alg <- glmmTMB(f.abund ~ treatment * as.factor(sampling) +
  a.abund.c + sg.sd.c + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# full model
fa.full <- glmmTMB(f.abund ~ treatment * as.factor(sampling) +
  a.abund.c + sg.sd.c + sg.prod.c + (1 | plot),
  family = compois(link = "log"),
data = fish.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# check residuals
glmm.resids(fa.treat)
glmm.resids(fa.prod) 
glmm.resids(fa.struct)
glmm.resids(fa.alg) 
glmm.resids(fa.treat.prod)
glmm.resids(fa.treat.struct)
glmm.resids(fa.treat.alg)
glmm.resids(fa.prod.alg)
glmm.resids(fa.prod.struct)
glmm.resids(fa.struct.alg)
glmm.resids(fa.treat.prod.alg)
glmm.resids(fa.treat.prod.struct) 
glmm.resids(fa.treat.struct.alg)
glmm.resids(fa.full)



# Residuals are ok. Now do model selection

# Model selection
abund.cand.mod.names <- c(
  "fa.treat",
  "fa.prod",
  "fa.struct",
  "fa.alg",
  "fa.treat.prod",
  "fa.treat.struct",
  "fa.treat.alg",
  "fa.prod.alg",
  "fa.prod.struct",
  "fa.struct.alg",
  "fa.treat.prod.alg",
  "fa.treat.prod.struct",
  "fa.treat.struct.alg",
  "fa.full"
)
abund.cand.mods <- list()

# This function fills the list by model names
for (i in 1:length(abund.cand.mod.names)) {
  abund.cand.mods[[i]] <- get(abund.cand.mod.names[i])
}

# Function aictab does the AICc-based model comparison
(fabund.aic<-data.frame(aictab(
  cand.set = abund.cand.mods,
  modnames = abund.cand.mod.names
)))
write_rds(fabund.aic,"working_data/FishAbundAIC.rds")
# treatment only is the best model?

#to reference in the manuscript
(fa.treat.sum<-summary(fa.treat))
write_rds(fa.treat.sum,"working_data/FishAbundTreatSum.rds")

# confirm that the control plots don't change sig with time
fa.treat.f <- glmmTMB(f.abund ~ treatment * as.factor(sampling) +
    (1 | plot),
  family = compois(link = "log"),
  data = fish.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "fake"))
)
fa.treat.b <- glmmTMB(f.abund ~ treatment * as.factor(sampling) +
    (1 | plot),
  family = compois(link = "log"),
  data = fish.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "blank"))
)
summary(fa.treat.f)
summary(fa.treat.b)
