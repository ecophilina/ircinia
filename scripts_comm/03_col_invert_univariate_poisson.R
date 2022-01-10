# Clonal inverts univariate analysis
library(vegan)
library(lmerTest)
library(glmmTMB)
library(DHARMa)
library(pacman)
pacman::p_load(s20x, lme4, AICcmodavg, MASS)
library(tidyverse)

source("scripts_comm/02_community_data_org.R")

# #view total number of invert taxa observed in plots (15)
# inverts %>% 
#   filter(abundance > 0) %>% 
#   select(taxa) %>% 
#   distinct() %>% 
#   View()

# visualize data
# ggplot(col.inv.uni) +
#   geom_jitter(aes(
#     x = sampling,
#     y = spr,
#     color = treatment,
#     group = treatment))
# 
# ggplot(col.inv.uni) +
#   geom_jitter(aes(
#     x = sampling,
#     y = i.abund,
#     color = treatment,
#     group = treatment))
# 
# ggplot(col.inv.uni) +
#   geom_jitter(aes(
#     x = sampling,
#     y = j,
#     color = treatment,
#     group = treatment))
# 
# ggplot(col.inv.uni) +
#   geom_jitter(aes(
#     x = sampling,
#     y = div,
#     color = treatment,
#     group = treatment))

# organize data to use offset

col.inv.uni.0 <- col.inv.uni %>%
  filter(sampling == 0) %>%
  dplyr::select(
    treatment,
    plot,
    start.spr = spr,
    start.div = div,
    start.j = j,
    start.a = i.abund)

col.inv.uni2 <- col.inv.uni %>%
  filter(sampling != 0) %>%
  left_join(col.inv.uni.0) %>%
  mutate(
    change.spr = spr - start.spr,
    change.div = div - start.div,
    change.j = j - start.j,
    change.a = i.abund - start.a)

col.inv.uni2$treatment <- factor(col.inv.uni2$treatment)
col.inv.uni$treatment <- factor(col.inv.uni$treatment)

# look at distribution of potential response variables
hist(col.inv.uni2$change.spr)
hist(col.inv.uni2$change.a)
hist(col.inv.uni2$change.div)
hist(col.inv.uni2$change.j)

# make sure there aren't differences between treatments originally

#Species Richness
summary(aov(spr ~ treatment, data = col.inv.uni %>%
              filter(sampling == 0)))

# No difference

# Abundance

summary(aov(i.abund ~ treatment, data = col.inv.uni %>%
              filter(sampling == 0)))

# no difference

# diversity

summary(aov(div ~ treatment, data = col.inv.uni %>%
              filter(sampling == 0)))

# no difference

# evenness

summary(aov(j ~ treatment, data = col.inv.uni %>%
              filter(sampling == 0)))
# No difference


## species richness model selection ####

#Treatment only model
cispr.treat <- glmmTMB(spr ~ treatment * as.factor(sampling) +
                        (1 | plot),
                      family = compois(link = "log"),
                      data = col.inv.uni %>%
                        filter(season == "summer") %>%
                        mutate(treatment = relevel(treatment, ref = "real")))

# seagrass productivity model for richness
cispr.prod <- glmmTMB(spr ~ as.factor(sampling) + sg.prod.c + 
                       (1 | plot),
                     family = compois(link = "log"),
                     data = col.inv.uni %>%
                       filter(season == "summer") %>%
                       mutate(treatment = relevel(treatment, ref = "real")))

# seagrass structure model for richness
cispr.struct <- glmmTMB(spr ~ as.factor(sampling) + sg.sd.c + 
                         (1 | plot),
                       family = compois(link = "log"),
                       data = col.inv.uni %>%
                         filter(season == "summer") %>%
                         mutate(treatment = relevel(treatment, ref = "real")))

#algal abundance model for richness
cispr.alg <- glmmTMB(spr ~ as.factor(sampling) + a.abund.c + 
                      (1 | plot),
                    family = compois(link = "log"),
                    data = col.inv.uni %>%
                      filter(season == "summer") %>%
                      mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment and seagrass productivity model
cispr.treat.prod <- glmmTMB(spr ~ treatment * as.factor(sampling) +
                             sg.prod.c +(1 | plot),
                           family = compois(link = "log"),
                           data = col.inv.uni %>%
                             filter(season == "summer") %>%
                             mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment and struct for richness
cispr.treat.struct <- glmmTMB(spr ~ treatment * as.factor(sampling) +
                               sg.sd.c + (1 | plot),
                             family = compois(link = "log"),
                             data = col.inv.uni %>%
                               filter(season == "summer") %>%
                               mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment and algal abundance for richness
cispr.treat.alg <- glmmTMB(spr ~ treatment * as.factor(sampling) +
                            a.abund.c + (1 | plot),
                          family = compois(link = "log"),
                          data = col.inv.uni %>%
                            filter(season == "summer") %>%
                            mutate(treatment = relevel(treatment, ref = "real")))

#combine productivity and algal abundance for richness
cispr.prod.alg <- glmmTMB(spr ~ as.factor(sampling)  +
                           a.abund.c + sg.prod.c + (1 | plot),
                         family = compois(link = "log"),
                         data = col.inv.uni %>%
                           filter(season == "summer") %>%
                           mutate(treatment = relevel(treatment, ref = "real")))

# combine seagrass productivity and struct for richness
# Note that productivity and struct are correlated with each other so their 
# individual contributions can't be assessed in the full model, just their combined effect:

# plot(sg.prod.c~sg.sd.c, data = col.inv.uni %>% filter(season == "summer") %>% mutate(treatment = relevel(treatment, ref = "real")) )


cispr.prod.struct <- glmmTMB(spr ~ as.factor(sampling)  +
                              sg.sd.c + sg.prod.c + (1 | plot),
                            family = compois(link = "log"),
                            data = col.inv.uni %>%
                              filter(season == "summer") %>%
                              mutate(treatment = relevel(treatment, ref = "real")))

#combine structure and algal abundance for richness
cispr.struct.alg <- glmmTMB(spr ~ as.factor(sampling)  +
                             sg.sd.c + a.abund.c + (1 | plot),
                           family = compois(link = "log"),
                           data = col.inv.uni %>%
                             filter(season == "summer") %>%
                             mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment, productivity, and algal abundance for richenss
cispr.treat.prod.alg <- glmmTMB(spr ~ treatment * as.factor(sampling)  +
                                 sg.prod.c + a.abund.c + (1 | plot),
                               family = compois(link = "log"),
                               data = col.inv.uni %>%
                                 filter(season == "summer") %>%
                                 mutate(treatment = relevel(treatment, ref = "real")))

# combine treatment, productivity, and seagrass structure for richness
cispr.treat.prod.struct <- glmmTMB(spr ~ treatment * as.factor(sampling) +
                                    sg.prod.c + sg.sd.c + (1 | plot),
                                  family = compois(link = "log"),
                                  data = col.inv.uni %>%
                                    filter(season == "summer") %>%
                                    mutate(treatment = relevel(treatment, ref = "real")))

#combine treatment, seagrass structure, and algal abundance for richness
cispr.treat.struct.alg <- glmmTMB(spr ~ treatment * as.factor(sampling) +
                                   sg.sd.c + a.abund.c + (1| plot),
                                 family = compois(link = "log"),
                                 data = col.inv.uni %>%
                                   filter(season == "summer") %>%
                                   mutate(treatment = relevel(treatment, ref = "real")))

# full model for richness (treatment, seagrass structure, seagrass productivity, and algal abundance)
cispr.full <- glmmTMB(spr ~ treatment * as.factor(sampling) +
                       sg.sd.c + sg.prod.c + a.abund.c + (1|plot),
                     family = compois(link = "log"),
                     data = col.inv.uni %>%
                       filter(season == "summer") %>%
                       mutate(treatment = relevel(treatment, ref = "real")))

# check residuals
glmm.resids(cispr.treat)
glmm.resids(cispr.prod)
glmm.resids(cispr.struct)
glmm.resids(cispr.alg) 
glmm.resids(cispr.treat.prod) 
glmm.resids(cispr.treat.struct)
glmm.resids(cispr.treat.alg) 
glmm.resids(cispr.prod.alg)
glmm.resids(cispr.prod.struct)
glmm.resids(cispr.struct.alg)
glmm.resids(cispr.treat.prod.alg) 
glmm.resids(cispr.treat.prod.struct) 
glmm.resids(cispr.treat.struct.alg) 
glmm.resids(cispr.full) 

# #algal abundance model without time variable
# cispr.alg.notime <- glmmTMB(spr ~ a.abund.c + 
#     (1 | plot),
#   family = compois,
#   data = col.inv.uni %>%
#     filter(season == "summer") %>%
#     mutate(treatment = relevel(treatment, ref = "real")))
# 
# glmm.resids(cispr.alg.notime)
# 

# Model selection
spr.cand.mod.names <- c(
  "cispr.treat",
  "cispr.prod",
  "cispr.struct",
  "cispr.alg",
  # "cispr.alg.notime", 
  "cispr.treat.prod",
  "cispr.treat.struct",
  "cispr.treat.alg",
  "cispr.prod.alg",
  "cispr.prod.struct",
  "cispr.struct.alg",
  "cispr.treat.prod.alg",
  "cispr.treat.prod.struct",
  "cispr.treat.struct.alg",
  "cispr.full")
spr.cand.mods <- list( ) 

# This function fills the list by model names
for(i in 1:length(spr.cand.mod.names)) {
  spr.cand.mods[[i]] <- get(spr.cand.mod.names[i]) }

# Function aictab does the AICc-based model comparison
col.aic<-data.frame(aictab(cand.set = spr.cand.mods, 
             modnames = spr.cand.mod.names))

write_rds(col.aic,"working_data/ColInvAIC_Results.rds")

# top models of clonal invert spr is structure, alg, and productivity all within 2 delta AICcfrom each other
(cispr.struct.sum<-summary(cispr.struct))
(cispr.prod.sum<-summary(cispr.prod))
(cispr.alg.sum<-summary(cispr.alg))

write_rds(cispr.struct.sum,"working_data/cispr_struct_sum.rds")
write_rds(cispr.prod.sum,"working_data/cispr_prod_sum.rds")
write_rds(cispr.alg.sum,"working_data/cispr_alg_sum.rds")

#nothing is significant

ggplot(data = col.inv.uni %>% filter(season == "summer") %>% 
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
  ylab("Clonal Invertebrate species richness") 
  # coord_cartesian(expand = F) +
  #ggsidekick::theme_sleek()
ggsave("clonal-invert-spr-by-algae.png", width = 5, height = 3)


# confirm that the control plots don't change sig with time
cispr.treat.f <- glmmTMB(spr ~ treatment *  as.factor(sampling) +
                          (1 | plot),
                        family = compois,
                        data = col.inv.uni %>%
                          filter(season == "summer") %>%
                          mutate(treatment = relevel(treatment, ref = "fake")))

cispr.treat.b <- glmmTMB(spr ~ treatment * as.factor(sampling) +
                          (1 | plot),
                        family = compois,
                        data = col.inv.uni %>%
                          filter(season == "summer") %>%
                          mutate(treatment = relevel(treatment, ref = "blank")))

summary(cispr.treat.f)
summary(cispr.treat.b) #blank at 1 month probably increased by chance, diff gone by 12 months
 #no difference