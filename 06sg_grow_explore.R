# How does the Ircinia sponge influence seagrass growth? #
#Here I will explore the sg_growth data set#
# run 03_reimport.R fist #

# packages----
## @knitr loadpackages
# load packages #

if(!require(here))install.packages('here');library(here)
if(!require(tidyverse))install.packages('tidyverse');library(tidyverse)
if(!require(readxl))install.packages('readxl');library(readxl)

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
  group_by(plot,sampling,treatment='real',dist) %>%
  summarise(mgd=mean(gd), sdgd=sd(gd)) %>%
  ggplot()

bp+
  geom_point(aes(x=sampling, y=mgd, color=treatment))+
  geom_line(aes(x=sampling, y=mgd, color=treatment, group=plot))+
  facet_wrap(~dist)


##here I am attempting to do the same as just above (viewing mgd over time and distance) but only looking at the fake sponge.
bp<-sg_grow %>% 
  group_by(plot,sampling,treatment='fake',dist) %>%
  summarise(mgd=mean(gd), sdgd=sd(gd)) %>%
  ggplot()

bp+
  geom_point(aes(x=sampling, y=mgd, color=treatment))+
  geom_line(aes(x=sampling, y=mgd, color=treatment, group=plot))+
  facet_wrap(~dist)

##here I am attempting to do the same as just above (viewing mgd over time and distance) but only looking at the blank treatment.
bp<-sg_grow %>% 
  group_by(plot,sampling,treatment='blank',dist) %>%
  summarise(mgd=mean(gd), sdgd=sd(gd)) %>%
  ggplot()

bp+
  geom_point(aes(x=sampling, y=mgd, color=treatment))+
  geom_line(aes(x=sampling, y=mgd, color=treatment, group=plot))+
  facet_wrap(~dist)

#Now I would like to review mgd over time ONLY reviewing summer time of sampling (1,2,and 4)

bp<-sg_grow %>% 
  group_by(plot,sampling== (1, 2, 4 ),treatment,dist) %>%
  summarise(mgd=mean(gd), sdgd=sd(gd)) %>%
  ggplot()

bp+
  geom_point(aes(x=sampling== 1,2,4, y=mgd, color=treatment))+
  geom_line(aes(x=sampling, y=mgd, color=treatment, group=plot))+
  facet_wrap(~dist)