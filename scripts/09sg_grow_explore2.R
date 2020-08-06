# How does the Ircinia sponge influence seagrass growth? #
#Here I will explore the sg_growth data set#

# packages----
## @knitr loadpackages
# load packages #

if(!require(here))install.packages('here');library(here)
if(!require(tidyverse))install.packages('tidyverse');library(tidyverse)
if(!require(readxl))install.packages('readxl');library(readxl)
#imports all the data sets with awesome new code!
source("scripts/03_reimport.R")

# is there a sig dif in seagrass growth per day among blank, fake, 
# and real treatment over time

sg_grow$gpd<-sg_grow$total.growth.mm2/sg_grow$days

ggplot(sg_grow%>%
         filter(sampling %in% c(1,4))%>%
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
  filter(sampling==1)%>%
  rename(start_gr = gpd) %>% 
  select(plot, dist, treatment,start_gr)
sgg2<-left_join(sg_grow,sv)%>%
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
      sampling==5~2))

# looking at the effect sampling at dist 0.5 with offset
sgmd05<-lmer(gpd~yr*treatment*season+(1|plot), offset=start_gr,data=sgg2%>%
               mutate(treatment=relevel(treatment,ref = "fake"))%>%
               filter(dist==0.5))
summary(sgmd05)

# looking at the effect sampling at dist 0 with offset
sgmd0<-lmer(gpd~yr*treatment*season+(1|plot), offset=start_gr,data=sgg2%>%
               mutate(treatment=relevel(treatment,ref = "fake"))%>%
               filter(dist==0))
summary(sgmd0)

# looking at the effect sampling at dist 1 with offset
sgmd0<-lmer(gpd~yr*treatment*season+(1|plot), offset=start_gr,data=sgg2%>%
              mutate(treatment=relevel(treatment,ref = "fake"))%>%
              filter(dist==1))
summary(sgmd0)
