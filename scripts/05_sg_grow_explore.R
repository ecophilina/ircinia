# How does the Ircinia sponge influence seagrass growth? #
#Here I will explore the sg_growth data set#

# packages----
## @knitr loadpackages
# load packages #

if(!require(here))install.packages('here');library(here)
if(!require(tidyverse))install.packages('tidyverse');library(tidyverse)
if(!require(readxl))install.packages('readxl');library(readxl)
if(!require(lmerTest))install.packages('lmerTest');library(lmerTest)
#imports all the data sets with awesome new code!
source("scripts/03_reimport.R")

# is there a sig dif in seagrass growth per day among blank, fake, 
# and real treatment over time

sg_grow$gpd<-sg_grow$total.growth.mm2/sg_grow$days

sv<-sg_grow%>%
  filter(sampling==1)%>% group_by(plot,dist,treatment) %>%
  mutate(start_gr = mean(gpd)) %>% 
  select(plot, dist, treatment, start_gr) %>% distinct()

# add lagged values as possible predictors
sg_lag <- sg_grow %>% group_by(plot, dist, sampling, treatment) %>% 
  summarise(mean_gpd=mean(gpd, na.rm = T)) %>%
  ungroup() %>% group_by(plot,dist) %>%
  arrange(sampling, .by_group = TRUE) %>% 
  mutate(
    previous_gpd = lag(mean_gpd, order_by = sampling),
    previous_gpd2 = lag(previous_gpd, order_by = sampling)
  )
sgg1<- left_join(sg_grow, sg_lag) 

sgg2<- left_join(sgg1, sv) %>% 
  filter(sampling!=1)%>%
  mutate(season=case_when(
    sampling==2~"s",
    sampling==3~"w",
    sampling==4~"s",
    sampling==5~"w"),
    yr=case_when(
      sampling==2~1,
      sampling==3~1,
      sampling==4~2,
      sampling==5~2),
    delta_gpd = gpd - previous_gpd,
    delta_gpd_st = gpd - start_gr,
    dist_factor = case_when(
      dist >= 1~"farther",
      # dist < 1~"nearer"),
      dist < 1~"closer"),
    stake_id = paste0(plot, "-", dist)
  )



## distance pooling
# at 0 and 0.5 both real and fake have increased growth
sgmd0s<-lmer(gpd~ treatment + yr + season +  (1|plot),
  offset=start_gr,
  data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "real"))
  # data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "fake"))
  %>% filter(dist == 0)
)

(sg0<-anova(sgmd0s))

sgmd0s<-lmer(gpd~ treatment + yr + season +  (1|plot),
  offset=start_gr,
  data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "real"))
  # data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "fake"))
  %>% filter(dist == 0.5)
)

(sg05<-anova(sgmd0s))

# at >= 1: no differences
sgmd0s<-lmer(gpd~ treatment + yr + season * previous_gpd + (1|plot), 
  data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "real"))
  # data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "fake"))
  %>% filter(dist == 1)
)
(sg1<-anova(sgmd0s))



#### removing previous growth rate makes sense because seasonal difference explains most of the relationship 
sgg2$treatment<-factor(sgg2$treatment)
sgmd0s<-lmer(gpd~ 
    treatment * dist_factor + 
    treatment * yr + 
    # dist_factor * yr +
    treatment * season + 
              (1|plot) + # adding this explains so little variation that it doesn't change anything...
              (1|stake_id), 
             offset=start_gr,
             # data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "blank"))
             data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "real"))
             # data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "fake"))
)


# (sgr.aov<-glmmTMB:::Anova.glmmTMB(sgmd0s, type = "III"))
(sgr<-summary(sgmd0s))

sgmd0sb<-lmer(gpd~ treatment * dist_factor+treatment * yr + treatment*season   + (1|stake_id), 
              offset=start_gr,
              data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "blank"))
              #data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "real"))
              # data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "fake"))
)
sgb<-summary(sgmd0sb)
sgmd0sf<-lmer(gpd~ treatment * dist_factor+treatment * yr + treatment*season  + 
    (1|plot) + # adding this explains so little variation that it doesn't change anything...
    (1|stake_id), 
              offset=start_gr,
              # data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "blank"))
              #data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "real"))
              data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "fake"))
)
sgf<-summary(sgmd0sf)
(sgaov<-anova(sgmd0sf))

