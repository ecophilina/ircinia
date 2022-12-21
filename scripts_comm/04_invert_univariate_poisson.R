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
inv.uni$treatment <- factor(inv.uni$treatment)

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
#   family = compois(link = "log"),
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
(inv.spraic<-data.frame(aictab(cand.set = spr.cand.mods, 
             modnames = spr.cand.mod.names)))

write_rds(inv.spraic,"working_data/InvSprAIC_Results.rds")

# top models of invert spr ####
(ispr.treat.sum<-summary(ispr.treat))
(ispr.alg.sum<-summary(ispr.alg))

write_rds(ispr.treat.sum,"working_data/InvSprTreatSum.rds")
write_rds(ispr.alg.sum,"working_data/InvSprAlgSum.rds")

# treatment is less than 2 delta AIC better than the algae model, this suggest that algae 


#algal abundance model for richness
ispr.alg2 <- glmmTMB(spr ~ a.abund.c + 
    (1 | plot),
  family = compois(link = "log"),
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "real")))
summary(ispr.alg2)

# make dataset for plots
alg.pr<-ggeffects::ggpredict(ispr.alg,terms=c("a.abund.c", "sampling"))%>%
  filter(group == 12)
alg.pr$x <- (alg.pr$x * inv.uni$abund.global.sd[1]*2 )+ inv.uni$a.abund.global[1]

ggplot(data = inv.uni %>% filter(season == "summer") %>% 
    mutate(treatment = factor(treatment, 
      levels = c("blank", "fake", "real"), 
      labels = c("Control", "Structure", "Sponge")))) + 
  geom_ribbon(data = alg.pr,
              aes(x , ymin = conf.low, ymax = conf.high),alpha=.2)+
  geom_jitter(aes(a.abund, spr, colour = treatment, alpha = as.factor(sampling)), 
    width = 0.3, height = 0.1, size = 2.5) + 
  geom_line(data=alg.pr,aes(x=x,y=predicted))+
  scale_alpha_discrete(range = c(0.3,1), name = "Months") + 

  scale_color_viridis_d(option="A", begin=0, end=0.65,"")+
  xlab("Algae Abundance") +
  ylab("Invertebrate Taxa Richness") +
  coord_cartesian(expand = F, ylim = c(-0.1,7.5), xlim = c(-0.1,35.35)) +
  ggsidekick::theme_sleek()
ggsave("figures/fig3.tiff", units=c("mm"),dpi=600,width = 129, height = 100)
ggsave("figures/invert_spr_by_algae.png", width = 5, height = 3)
ggsave("figures/invert_spr_by_algae.pdf", width = 5, height = 3)


# treatment model for comparison
summary(ispr.treat)
# confirm that the control plots don't change sig with time
ispr.treat.f <- glmmTMB(spr ~ treatment *  as.factor(sampling) +
    (1 | plot),
  family = compois(link = "log"),
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "fake")))

ispr.treat.b <- glmmTMB(spr ~ treatment * as.factor(sampling) +
    (1 | plot),
  family = compois(link = "log"),
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "blank")))

summary(ispr.treat.f)
summary(ispr.treat.b) #blank at 1 month probably increased by chance, diff gone by 12 months



# model selection for abundance ####


# Treatment only model
ia.treat <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
  (1 | plot),
family = compois(link = "log"),
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# seagrass productivity model
ia.prod <- glmmTMB(i.abund ~ as.factor(sampling) +
  sg.prod.c + (1 | plot),
family = compois(link = "log"),
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# seagrass structure model
ia.struct <- glmmTMB(i.abund ~ as.factor(sampling) +
  sg.sd.c + (1 | plot),
family = compois(link = "log"),
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# algal abundance model for abundance
ia.alg <- glmmTMB(i.abund ~ a.abund.c + as.factor(sampling) +
  (1 | plot),
family = compois(link = "log"),
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)
# combine treatment and seagrass productivity model
ia.treat.prod <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
  sg.prod.c + (1 | plot),
family = compois(link = "log"),
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment and struct for abundance
ia.treat.struct <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
  sg.sd.c + (1 | plot),
family = compois(link = "log"),
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment and algal abundance for abundance
ia.treat.alg <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
  a.abund.c + (1 | plot),
family = compois(link = "log"),
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine seagrass productivity and algal abundance for abundance
ia.prod.alg <- glmmTMB(i.abund ~ as.factor(sampling) +
  sg.prod.c + a.abund.c + (1 | plot),
family = compois(link = "log"),
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
family = compois(link = "log"),
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment, productivity, and algal abundance for abundance
ia.treat.prod.alg <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
  sg.prod.c + a.abund.c + (1 | plot),
family = compois(link = "log"),
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment, productivity and structure for abundance
ia.treat.prod.struct <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
  sg.prod.c + sg.sd.c + (1 | plot),
family = compois(link = "log"),
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)

# combine treatment, structure, and algal abundance for abundance
ia.treat.struct.alg <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
  a.abund.c + sg.sd.c + (1 | plot),
family = compois(link = "log"),
data = inv.uni %>%
  filter(season == "summer") %>%
  mutate(treatment = relevel(treatment, ref = "real"))
)
# full model for abundance
ia.full <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
  sg.prod.c + sg.sd.c + a.abund.c + (1 | plot),
family = compois(link = "log"),
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
(inv.abund.aic<-data.frame(aictab(cand.set = abund.cand.mods, 
             modnames = abund.cand.mod.names)))

write_rds(inv.abund.aic,"working_data/InvAbundAIC_Results.rds")

# we have three models at the top - treat * sg prod, treat, treat * sg struct

(ia.treat.sum<-summary(ia.treat))
(ia.treat.prod.sum<-summary(ia.treat.prod))
(ia.treat.struct.sum<-summary(ia.treat.struct))

write_rds(ia.treat.sum,"working_data/InvAbundTreatSum.rds")
write_rds(ia.treat.prod.sum,"working_data/InvAbundTreatProdSum.rds")
write_rds(ia.treat.struct.sum,"working_data/InvAbundTreatStructSum.rds")

#to reference in manuscript

# confirm that the control plots don't change sig with time
ia.treat.prod.f <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
    sg.prod.c + (1 | plot),
  family = compois(link = "log"),
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "fake"))
)
ia.treat.prod.b <- glmmTMB(i.abund ~ treatment * as.factor(sampling) +
    sg.prod.c + (1 | plot),
  family = compois(link = "log"),
  data = inv.uni %>%
    filter(season == "summer") %>%
    mutate(treatment = relevel(treatment, ref = "blank"))
)

summary(ia.treat.prod.f)
summary(ia.treat.prod.b)


# calculate % changes

# abundance of inverts considering only treatment
(ia.treat.prop <- as.data.frame(ia.treat.sum$coefficients$cond))
(ia.0 <- exp(ia.treat.prop[1,1]))
(ia.1 <- exp(ia.treat.prop[1,1]+ ia.treat.prop[4,1]))
(ia.12 <- exp(ia.treat.prop[1,1]+ ia.treat.prop[5,1]))

# change at 1 month
ia.1/ia.0
# change at 12 months
ia.12/ia.0

# control plot at 12 months
(ia.b.12 <- exp(ia.treat.prop[1,1]+ ia.treat.prop[2,1]+ia.treat.prop[5,1]+ia.treat.prop[8,1]))
# how much less than sponge
ia.12/ia.b.12

# structure plot at 12 months
(ia.f.12 <- exp(ia.treat.prop[1,1]+ ia.treat.prop[3,1]+ia.treat.prop[5,1]+ia.treat.prop[9,1]))
# how much less than sponge
ia.12/ia.f.12

(ia.treat.prod.prop <- as.data.frame(ia.treat.prod.sum$coefficients$cond))

(ia.prod.0 <- exp(ia.treat.prod.prop[1,1]))
(ia.prod.1 <- exp(ia.treat.prod.prop[1,1]+ ia.treat.prod.prop[4,1]))
(ia.prod.12 <- exp(ia.treat.prod.prop[1,1]+ ia.treat.prod.prop[5,1]))

# change at 1 month
ia.prod.1/ia.prod.0
# change at 12 months
ia.prod.12/ia.prod.0

# effect of productivity
# (ia.prod <- exp(ia.treat.prod.prop[1,1] + ia.treat.prod.prop[5,1] + ia.treat.prod.prop[6,1]))
# (ia.prod.2 <- exp(ia.treat.prod.prop[1,1] + ia.treat.prod.prop[5,1]))
# ia.prod/ia.prod.2
(ia.prod.only <- exp(ia.treat.prod.prop[6,1]))

# get missing B and CI f
decimalplaces <- function(x) {
  if ((x - round(x)) != 0) {
    strs <- strsplit(as.character(format(x, scientific = F)), "\\.")
    n <- nchar(strs[[1]][2])
  } else {
    n <- 0
  }
  return(n) 
}

efsize.tmb<-function(rtable,row){
  betas<-signif(rtable$coefficients$cond[row,1], 2)
  betas<-ifelse(abs(betas)<10,
                formatC(signif(betas,2), digits=2, format="fg", flag="#"), 
                round(betas))
  betas2 <- as.numeric(betas)
  decicount <- ifelse(betas2>1 & betas2 <10, 1, decimalplaces(betas2))
  se<-rtable$coefficients$cond[row,2]
  cil<-rtable$coefficients$cond[row,1]-1.96*se
  cil<-format(round(cil, digits=decicount), scientific=F)
  cih<-rtable$coefficients$cond[row,1]+1.96*se
  cih<-format(round(cih, digits=decicount), scientific=F)
  return(paste0(betas,", CI = ",cil," to ",cih))
}

efsize.tmb(ia.treat.prod.sum,5)
# efsize.tmb(fspr.treat.sum,9)

# richness of inverts considering only treatment
ispr.treat.prop <- as.data.frame(ispr.treat.sum$coefficients$cond)
(ispr.0 <- exp(ispr.treat.prop[1,1]))
(ispr.1 <- exp(ispr.treat.prop[1,1]+ ispr.treat.prop[4,1]))
(ispr.12 <- exp(ispr.treat.prop[1,1]+ ispr.treat.prop[5,1]))

ispr.1/ispr.0
ispr.12/ispr.0

(ispr.b.12 <- exp(ispr.treat.prop[1,1]+ ispr.treat.prop[2,1]+ispr.treat.prop[5,1]+ispr.treat.prop[8,1]))

ispr.12/ispr.b.12

(ispr.f.12 <- exp(ispr.treat.prop[1,1]+ ispr.treat.prop[3,1]+ispr.treat.prop[5,1]+ispr.treat.prop[9,1]))

ispr.12/ispr.f.12

# effect of algae
(ispr.alg.prop <- as.data.frame(ispr.alg.sum$coefficients$cond))
(ispr.alg <- exp(ispr.alg.prop[4,1]))



