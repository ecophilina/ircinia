# How does the Ircinia sponge influence seagrass growth? #
#Here I will explore the sg_growth data set#
# run 03_reimport.R fist #

# packages----
## @knitr loadpackages
# load packages #

if(!require(here))install.packages('here');library(here)
if(!require(tidyverse))install.packages('tidyverse');library(tidyverse)
if(!require(readxl))install.packages('readxl');library(readxl)
#imports all the data sets with awesome new code!
source("scripts/03_reimport.R")

#Before exploring, we need to make the growth per day added as a new column because collection days were not done in even intervals
sg_grow$gd <- sg_grow$total.growth.mm2/sg_grow$days


#creating an object for ggplot
bp<-ggplot(data=sg_grow)


#adding seasons to the plots
sg_grow <- sg_grow %>% 
  mutate(season=case_when(
    sampling==1~"summer",
    sampling==2~"summer",
    sampling==3~"winter",
    sampling==4~"summer",
    sampling==5~"winter" ))

#Here I created a boxplot of growth per day by treatment with four different distances. It appeareds that for both blank and fake treatments, the gd decreased as the distance increased. However, with the real sponge treatment there is little to no change in gd with increasing distance. 
bp+
  geom_boxplot(aes(x=treatment,y=gd))+ 
  facet_wrap(~dist)

#Here I created a boxplot of gd by treatment sorted into seasons and distance. As we noted previously, there is higher growth rates in summer seasons than winter. In this boxplot, it is hard to see if there is any change in gd with increasing distance. 
bp+
  geom_boxplot(aes(x=treatment,y=gd))+ 
  facet_grid(dist~season)
#Here I created a boxplot of gd by treatment incorporating distance and sampling. It looks really gross and is hard to tell if there is a trend.
bp+
  geom_boxplot(aes(x=treatment,y=gd))+ 
  facet_grid(dist~sampling)
#Here I plotted mean growth per day by sampling time separated by treatment. The blank treatments seem to start lower than the fake and real and end lower in terms of mgd as well.
#changing the data set to allow for mean gd and sd

bp<-sg_grow %>% filter(dist==0.5) %>% 
  group_by(plot,sampling,treatment) %>%
  summarise(mgd=mean(gd), sdgd=sd(gd)) %>%
  ggplot()

bp+
  geom_point(aes(x=sampling, y=mgd, color=treatment))+
  geom_line(aes(x=sampling, y=mgd, color=treatment, group=plot))
# Here I plotted mgd by sampling time creating plots for each distance. 

bp<-sg_grow %>% 
  group_by(plot,sampling,treatment,dist) %>%
  summarise(mgd=mean(gd), sdgd=sd(gd)) %>%
  ggplot()

bp+
  geom_point(aes(x=sampling, y=mgd, color=treatment))+
  geom_line(aes(x=sampling, y=mgd, color=treatment, group=plot))+
  facet_wrap(~dist)

#here I am attempting to do the same as just above (viewing mgd over time and distance) but only looking at the real sponge. I realized that changing treatment to 'real' did not do what I wanted it to do but I got stuck----

bp<-sg_grow %>% 
  filter(treatment=='real') %>%
  group_by(plot,sampling,treatment,dist) %>%
  summarise(mgd=mean(gd), sdgd=sd(gd)) %>%
  ggplot()

bp+
  geom_point(aes(x=sampling, y=mgd, color=treatment))+
  geom_line(aes(x=sampling, y=mgd, color=treatment, group=plot))+
  facet_wrap(~dist)


##here I am attempting to do the same as just above (viewing mgd over time and distance) but only looking at the fake sponge.
bp<-sg_grow %>% 
  filter(treatment=='fake') %>%
  group_by(plot,sampling,treatment,dist) %>%
  summarise(mgd=mean(gd), sdgd=sd(gd)) %>%
  ggplot()

bp+
  geom_point(aes(x=sampling, y=mgd, color=treatment))+
  geom_line(aes(x=sampling, y=mgd, color=treatment, group=plot))+
  facet_wrap(~dist)

##here I am attempting to do the same as just above (viewing mgd over time and distance) but only looking at the blank treatment.
bp<-sg_grow %>% 
  filter(treatment=='blank') %>%
  group_by(plot,sampling,treatment,dist) %>%
  summarise(mgd=mean(gd), sdgd=sd(gd)) %>%
  ggplot()

bp+
  geom_point(aes(x=sampling, y=mgd, color=treatment))+
  geom_line(aes(x=sampling, y=mgd, color=treatment, group=plot))+
  facet_wrap(~dist)

#Now I would like to review mgd over time ONLY reviewing summer time of sampling (1,2,and 4) with the real treatment

bp<-sg_grow %>% 
  filter(sampling %in% c(1,2,4) & treatment=='real') %>%
  group_by(plot,sampling,treatment,dist) %>%
  summarise(mgd=mean(gd), sdgd=sd(gd)) %>%
  ggplot()

bp+
  geom_point(aes(x=sampling, y=mgd, color=treatment))+
  geom_line(aes(x=sampling, y=mgd, color=treatment, group=plot))+
  facet_wrap(~dist)

#Now I would like to review mgd over time ONLY reviewing winter time of sampling (3 and 5) with the real treatment
bp<-sg_grow %>% 
  filter(sampling %in% c(3,5)) %>%
  group_by(plot,sampling,treatment,dist) %>%
  summarise(mgd=mean(gd), sdgd=sd(gd)) %>%
  ggplot()

bp+
  geom_point(aes(x=sampling, y=mgd, color=treatment))+
  geom_line(aes(x=sampling, y=mgd, color=treatment, group=plot))+
  facet_wrap(~dist)

#This is additional explorational code I wanted to add after our meeting on 7/9/2020. 

#comparing growth per day by distance separated into treatment groups and sampling group
sg_grow %>% ggplot(aes(as.factor(dist), gd)) + geom_point(alpha=0.5) + 
  facet_grid(treatment~sampling)

sg_grow %>% filter(treatment=="blank") %>%
  ggplot(aes(dist, gd, color=as.factor(plot)))+
  geom_point() +
 # geom_line() +
  geom_smooth(method="lm", se=F) +
  facet_grid(dist ~ as.factor(sampling))

#Using sampling one as a start value and creating a new dataframe.
start_val <- filter(sg_grow, sampling == 1) %>%  
  #group_by(plot, dist) %>% 
  group_by(plot) %>% 
  mutate(start_gd = mean(gd)) %>% 
 # select(plot, dist, start_gd) %>% unique()
  select(plot, start_gd) %>% unique()
#taking the original dataframe and adding on th eaverage start value. For every value you it will repeat the start value of mean gd. Delta gd. sampling 1 is now our baseline so it is removed here
sg_grow2 <- left_join(sg_grow, start_val) %>% 
  mutate(delta_gd = gd - start_gd) %>% filter(sampling != 1)

#Now I am plotting the change in gd from the start value across distances.
sg_grow2 %>% 
  ggplot(aes(dist, delta_gd, color=as.factor(plot)))+
  geom_point() +
  # geom_line() +
  geom_smooth(method="lm", se=F) +
  facet_grid(sampling~ treatment)

#violin plot of change in gd from the start gd across distances.Not as informative as the plot but still fun:)
sg_grow2 %>% 
  ggplot(aes(as.factor(dist), delta_gd))+
  geom_violin()+
  geom_point() +
  facet_grid(sampling~ treatment)

#messing around with mixed models
install.packages('lmerTest')

pretreat <- sg_grow %>% filter(sampling==1)

#This is showing a significant difference of real treatment and distance from the blank treatment. We need to keep this in mind for future analysis.
mod <- lmerTest::lmer(gd~treatment * dist + (1|plot) , data = pretreat) 
summary(mod)

#this plot shows two plots increased in gd with increasing distance. This probably just happened by chance.
pretreat %>% filter(treatment=='real') %>%
  ggplot(aes(dist, gd, color=as.factor(plot))) +
  geom_point() +
  geom_smooth(method='lm', se=F)

#trying to do mixed effets modeling
sgg05<-lmerTest::lmer(gd~treatment,data=pretreat%>%filter(sampling==1&dist==0.5))
summary(pn05)
TukeyHSD(pn05)

ipn<-lmer(nvalue~treatment+as.factor(dist)+(1|plot),data=sgn%>%filter(sampling==1&nut=="PN"))
  

  