# first attempt at compositional vector analysis for fish

library(vegan)
library(lmerTest)
library(glmmTMB)
library(DHARMa)
library(tidyverse)

source("scripts_comm/02_community_data_org.R")
# subset down to only 0 and 12 month sampling data to start
fish.com2<-fish.com[fish.env$sampling %in% c(0,12),]
fish.env2<-fish.env[fish.env$sampling %in% c(0,12),]

# do hellinger transformation then RDA

f.com.hel<-decostand(fish.com2,"hellinger")

f.pca<-rda(f.com.hel)

# now extract coordinates for each row in the dataset and join this to environmental data

fish.env3<-bind_cols(fish.env2,data.frame(scores(f.pca,choices=c(1,2),display ="sites")))

# now organize data 

fish.env30<-fish.env3%>%
  filter(sampling==0)%>%
  dplyr::select(plot,start.x = PC1,start.y = PC2)

fish.env4<-fish.env3%>%
  filter(sampling==12)%>%
  dplyr::select(plot,treatment,end.x = PC1,end.y = PC2)%>%
  left_join(fish.env30)%>%
  mutate(dy = end.y - start.y,
         dx = end.x - start.x,
         vlength = sqrt(dx^2 +dy^2),
         vlength.rsc = scales::rescale(vlength,to = c(0,1),from = range(vlength)),
         angle = atan2(dy,dx),
         plot.start = 0,
         plot.end.x = cos(angle)*vlength.rsc,
         plot.end.y = sin(angle)*vlength.rsc)

# make plot
ggplot(data=fish.env4)+
  geom_segment(aes(x = plot.start, y = plot.start,
                   xend = plot.end.x, yend = plot.end.y,color=treatment),
               arrow = arrow(length = unit(0.5, "cm")))
