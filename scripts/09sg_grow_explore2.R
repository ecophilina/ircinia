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



ggplot(sg_grow%>%
         filter(sampling %in% c(2, 3, 4, 5))%>%
         group_by(sampling,plot,dist,treatment)%>%
         summarise(gpd=mean(gpd)))+
  geom_line(aes(x=sampling,y=gpd,group=plot,color=treatment))+
  facet_grid(treatment~dist,scales = "free")

ggplot(sg_grow%>%
         filter(sampling %in% c(4,5))%>%
         group_by(sampling,plot,dist,treatment)%>%
         summarise(gpd=mean(gpd)))+
  geom_line(aes(x=dist,y=gpd,group=plot,color=treatment))+
  facet_grid(treatment~sampling)

# now write the model.
sg_grow$treatment<-as.factor(sg_grow$treatment)
sgm<-lmer(gpd~as.factor(sampling)*dist*treatment+(1|plot),data=sg_grow%>%
            mutate(treatment=relevel(treatment,ref = "fake")))
summary(sgm)

# looking at the effect sampling at dist 0
sgmd<-lmer(gpd~as.factor(sampling)*treatment+(1|plot),data=sg_grow%>%
             mutate(treatment=relevel(treatment,ref = "fake"))%>%
             filter(dist==0))
summary(sgmd)

# looking at the effect sampling at dist 0.5
sgmd05<-lmer(gpd~as.factor(sampling)*treatment+(1|plot),data=sg_grow%>%
             mutate(treatment=relevel(treatment,ref = "fake"))%>%
             filter(dist==0.5))
summary(sgmd05)

# put sampling 1 in its own column and remove sampling 1 rows.
# then include sampling 1 as an offset.

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
    quad_id = paste0(plot, dist)
      )

# looking at the effect sampling at dist 0 with delta GPD
sgmd0ds<-lmer(delta_gpd_st~yr*treatment*season+(1|plot), data=sgg2%>%
    # mutate(treatment=relevel(treatment,ref = "fake"))%>%
    mutate(treatment=relevel(treatment,ref = "real"))%>%
    filter(dist==0))
summary(sgmd0d)

# looking at the effect sampling at dist 0 with offset
sgmd0s<-lmer(gpd~yr*treatment*season+(1|plot), offset=start_gr, data=sgg2%>%
    # mutate(treatment=relevel(treatment,ref = "fake"))%>%
    mutate(treatment=relevel(treatment,ref = "real"))%>%
    filter(dist==0))
summary(sgmd0s)

# an offset does exactly the same thing as the delta calculation

## BUT if delta calculation is previous time step instead of starting condition
# looking at the effect sampling at dist 0 with delta GPD
sgmd0d<-lmer(delta_gpd~yr*treatment*season+(1|plot), data=sgg2%>%
    # mutate(treatment=relevel(treatment,ref = "fake"))%>%
    mutate(treatment=relevel(treatment,ref = "real"))%>%
    filter(dist==0))
summary(sgmd0d)

sgmd0<-lmer(gpd ~ yr*treatment*season+ previous_gpd*season + (1|plot), offset=start_gr, data=sgg2%>%
    # mutate(treatment=relevel(treatment,ref = "fake"))%>%
    mutate(treatment=relevel(treatment,ref = "real"))%>%
    filter(dist==0))
summary(sgmd0)




# at 0 both real and fake have increased growth
sgmd0s<-lmer(gpd~ treatment + yr + season * previous_gpd + (1|plot), #+ (1|quad_id)
  offset=start_gr,
  data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "real"))
  # data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "fake"))
  %>% filter(dist == 0)
)

summary(sgmd0s)

# at 0.5 not different from blank but fake worse
sgmd0s<-lmer(gpd~ treatment + yr + season * previous_gpd + (1|plot), #+ (1|quad_id)
  offset=start_gr,
  data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "real"))
  # data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "fake"))
  %>% filter(dist == 0.5)
)

summary(sgmd0s)

# at >= 1: no differences
sgmd0s<-lmer(gpd~ treatment + yr + season * previous_gpd + (1|plot), #+ (1|quad_id)
  offset=start_gr,
  data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "real"))
  # data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "fake"))
  %>% filter(dist >= 1)
)

summary(sgmd0s)

# if both close distances pooled
sgmd0s<-lmer(gpd~ treatment + yr + season * previous_gpd + (1|plot), #+ (1|quad_id)
  offset=start_gr,
  data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "real"))
  # data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "fake"))
  %>% filter(dist < 1)
)

summary(sgmd0s)

#### most powerful option but feels a bit contrived...
# if both close (< 1m ) distances pooled and compared with further distances
sgmd0s<-lmer(gpd~ treatment * dist_factor + yr + season * previous_gpd + (1|plot), #+ (1|quad_id)
  offset=start_gr,
  # data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "blank"))
  data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "real"))
  # data=sgg2 %>% mutate(treatment=relevel(treatment, ref = "fake"))
  # %>% filter(dist < 2)
)

summary(sgmd0s)
