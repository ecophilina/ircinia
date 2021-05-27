# first attempt at compositional vector analysis for fish

library(vegan)
library(tidyverse)
library(patchwork)

source("scripts_comm/02_community_data_org.R")

# Fish 

# subset down to only 0 and 12 month sampling data to start
fish.com2<-fish.com[fish.env$sampling %in% c(0,1,12),]
fish.env2<-fish.env[fish.env$sampling %in% c(0,1,12),]

# do hellinger transformation then RDA

f.com.hel<-decostand(fish.com2,"hellinger")

f.pca<-rda(f.com.hel)

# now extract coordinates for each row in the dataset and join this to environmental data

fish.env3<-bind_cols(fish.env2,data.frame(scores(f.pca,choices=c(1,2),display ="sites")))%>%
  mutate(fillvar = paste0(treatment,".",sampling))

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
         plot.end.y = sin(angle)*vlength.rsc,
         plot.end.x2 = end.x - start.x,
         plot.end.y2 = end.y - start.y,
         lsize = ifelse(duplicated(vlength)&duplicated(angle),2,1))

# Look at differences in vector length between treatments
fish.vl.aov<-aov(vlength~treatment,data=fish.env4)
summary(fish.vl.aov)
TukeyHSD(fish.vl.aov)

fish.angle.aov<-aov(angle~treatment,data=fish.env4)
summary(fish.angle.aov)
TukeyHSD(fish.angle.aov)

# function to make a circle to put on plots
circleFun <- function(center = c(0, 0), diameter = 1, npoints = 100) {
  r <- diameter / 2
  tt <- seq(0, 2 * pi, length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

circ <- circleFun(center = c(0, 0), diameter = 2, npoints = 500)

# make plot
(fish.vplot<-ggplot(data=fish.env4)+
  geom_segment(aes(x = plot.start, y = plot.start,
                xend = plot.end.x, yend = plot.end.y,color=treatment,size = as.factor(lsize),alpha=lsize),
              # xend = plot.end.x2, yend = plot.end.y2,color=treatment,size = as.factor(lsize),alpha=lsize),
  arrow = arrow(length = unit(0.3, "cm")))+
  geom_path(data = circ, aes(x, y), lty = 2, alpha = 0.7) +
  theme_bw()+
  coord_fixed()+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        panel.border = element_blank())+
  scale_color_viridis_d(option="A", begin=0, end=0.6,name="Treatment",labels=c("Control","Structure","Sponge"))+
    scale_size_discrete(range=c(1,1.75),name = "Number of Plots")+
    scale_alpha_continuous(range=c(0.75,1),guide=FALSE)+
  geom_text(aes(x=-1,y=1),label="a",size=10))

fish.mid<-fish.env3%>%
  filter(sampling==1)%>%
  select(treatment, plot,mid.x=PC1,mid.y=PC2)

fish.segments<-left_join(fish.env4,fish.mid)

(fish.sp.space<-ggplot()+
  # geom_point(aes(x=PC1,y=PC2,color=treatment,shape=as.factor(sampling)),data=fish.env3,size=3,alpha=.5,position=position_dodge(0.01))+
  # geom_segment(aes(x=start.x,y=start.y,xend=mid.x,yend=mid.y,color=treatment),data=fish.segments)+
  # geom_segment(aes(x=mid.x,y=mid.y,xend=end.x,yend=end.y,color=treatment),data=fish.segments,  arrow = arrow(length = unit(0.3, "cm"))) +
    geom_point(aes(x=PC1,y=PC2,color=treatment,shape=as.factor(sampling)),data=fish.env3%>%filter(sampling!=1),size=3,alpha=.5,position=position_dodge(0.01))+
    geom_segment(aes(x=start.x,y=start.y,xend=end.x,yend=end.y,color=treatment),data=fish.segments,arrow = arrow(length = unit(0.3, "cm")))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_viridis_d(option="A", begin=0, end=0.6,name="Treatment",labels=c("Control","Structure","Sponge"))+
    scale_shape_discrete(name = "Sampling"))


# Non-colonial inverts

# subset down to only 0 and 12 month sampling data to start
inv.com2<-inv.com[inv.env$sampling %in% c(0,1,12),]
inv.env2<-inv.env[inv.env$sampling %in% c(0,1,12),]

# do hellinger transformation then RDA

i.com.hel<-decostand(inv.com2,"hellinger")

i.pca<-rda(i.com.hel)

# now extract coordinates for each row in the dataset and join this to environmental data

inv.env3<-bind_cols(inv.env2,data.frame(scores(i.pca,choices=c(1,2),display ="sites")))

# now organize data 

inv.env30<-inv.env3%>%
  filter(sampling==0)%>%
  dplyr::select(plot,start.x = PC1,start.y = PC2)

inv.env4<-inv.env3%>%
  filter(sampling==12)%>%
  dplyr::select(plot,treatment,end.x = PC1,end.y = PC2)%>%
  left_join(inv.env30)%>%
  mutate(dy = end.y - start.y,
         dx = end.x - start.x,
         vlength = sqrt(dx^2 +dy^2),
         vlength.rsc = scales::rescale(vlength,to = c(0,1),from = range(vlength)),
         angle = atan2(dy,dx),
         plot.start = 0,
         plot.end.x = cos(angle)*vlength.rsc,
         plot.end.y = sin(angle)*vlength.rsc,
         plot.end.x2 = end.x - start.x,
         plot.end.y2 = end.y - start.y,
         lsize = ifelse(duplicated(vlength)&duplicated(angle),2,1))

# Look at differences in vector length between treatments
inv.vl.aov<-aov(vlength~treatment,data=inv.env4)
summary(inv.vl.aov)


inv.angle.aov<-aov(angle~treatment,data=inv.env4)
summary(inv.angle.aov)

# make plot
(inv.plot<-ggplot(data=inv.env4%>%
                    filter(abs(plot.end.y)>0.000001 | abs(plot.end.x) >0.000001))+
  geom_segment(aes(x = plot.start, y = plot.start,
                   xend = plot.end.x, yend = plot.end.y,color=treatment),
                   # xend = plot.end.x2, yend = plot.end.y2,color=treatment),
               alpha=.75,size=1,
                  arrow = arrow(length = unit(0.3, "cm")))+
  geom_path(data = circ, aes(x, y), lty = 2, alpha = 0.7) +
  theme_bw()+
  coord_fixed()+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        panel.border = element_blank(),
        legend.position = "none")+
  scale_color_viridis_d(option="A", begin=0, end=0.6)+
  geom_text(aes(x=-1,y=1),label="b",size=10))

inv.mid<-inv.env3%>%
  filter(sampling==1)%>%
  select(treatment, plot,mid.x=PC1,mid.y=PC2)

inv.segments<-left_join(inv.env4,inv.mid)

(inv.sp.space<-ggplot()+
    # geom_point(aes(x=PC1,y=PC2,color=treatment,shape=as.factor(sampling)),data=inv.env3,size=3,alpha=.5,position=position_dodge(0.01))+
    # geom_segment(aes(x=start.x,y=start.y,xend=mid.x,yend=mid.y,color=treatment),data=inv.segments)+
    # geom_segment(aes(x=mid.x,y=mid.y,xend=end.x,yend=end.y,color=treatment),data=inv.segments,  arrow = arrow(length = unit(0.3, "cm"))) +
    geom_point(aes(x=PC1,y=PC2,color=treatment,shape=as.factor(sampling)),data=inv.env3%>%filter(sampling!=1),size=3,alpha=.5)+
    geom_segment(aes(x=start.x,y=start.y,xend=end.x,yend=end.y,color=treatment),data=inv.segments,arrow = arrow(length = unit(0.3, "cm")))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_viridis_d(option="A", begin=0, end=0.6,name="Treatment",labels=c("Control","Structure","Sponge"))+
    scale_shape_discrete(name = "Sampling"))


# colonial inverts

# subset down to only 0 and 12 month sampling data to start
col.inv.com2<-col.inv.com[col.inv.env$sampling %in% c(0,1,12),]
col.inv.env2<-col.inv.env[col.inv.env$sampling %in% c(0,1,12),]

# do hellinger transformation then RDA

ci.com.hel<-decostand(decostand(col.inv.com2,"pa"),"hellinger")

ci.pca<-rda(ci.com.hel)

# now extract coordinates for each row in the dataset and join this to environmental data

col.inv.env3<-bind_cols(col.inv.env2,data.frame(scores(ci.pca,choices=c(1,2),display ="sites")))

# now organize data

col.inv.env30<-col.inv.env3%>%
  filter(sampling==0)%>%
  dplyr::select(plot,start.x = PC1,start.y = PC2)

col.inv.env4<-col.inv.env3%>%
  filter(sampling==12)%>%
  dplyr::select(plot,treatment,end.x = PC1,end.y = PC2)%>%
  left_join(col.inv.env30)%>%
  mutate(dy = end.y - start.y,
         dx = end.x - start.x,
         vlength = sqrt(dx^2 +dy^2),
         vlength.rsc = scales::rescale(vlength,to = c(0,1),from = range(vlength)),
         angle = atan2(dy,dx),
         plot.start = 0,
         plot.end.x = cos(angle)*vlength.rsc,
         plot.end.y = sin(angle)*vlength.rsc,
         plot.end.x2 = end.x - start.x,
         plot.end.y2 = end.y - start.y,
         lsize = ifelse(duplicated(vlength)&duplicated(angle),2,1))

# Look at differences in vector length between treatments
col.inv.vl.aov<-aov(vlength~treatment,data=col.inv.env4)
summary(col.inv.vl.aov)

col.inv.angle.aov<-aov(angle~treatment,data=col.inv.env4)
summary(col.inv.angle.aov)

# make plot
(col.inv.plot<-ggplot(data=col.inv.env4%>%
                        filter(abs(plot.end.y)>0.000001 | abs(plot.end.x) >0.000001))+
    geom_segment(aes(x = plot.start, y = plot.start,
                     xend = plot.end.x, yend = plot.end.y,color=treatment,size = lsize,alpha=lsize),
                     # xend = plot.end.x2, yend = plot.end.y2,color=treatment,size = lsize,alpha=lsize),
                 arrow = arrow(length = unit(0.3, "cm")))+
    geom_path(data = circ, aes(x, y), lty = 2, alpha = 0.7) +
    theme_bw()+
    coord_fixed()+
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          line = element_blank(),
          panel.border = element_blank(),
          legend.position = "none")+
    scale_color_viridis_d(option="A", begin=0, end=0.6)+
    scale_size_continuous(range=c(1,1.75))+
    scale_alpha_continuous(range=c(0.75,1))+
    geom_text(aes(x=-1,y=1),label="c",size=10))

col.inv.mid<-col.inv.env3%>%
  filter(sampling==1)%>%
  select(treatment, plot,mid.x=PC1,mid.y=PC2)

col.inv.segments<-left_join(col.inv.env4,col.inv.mid)

(col.inv.sp.space<-ggplot()+
    # geom_point(aes(x=PC1,y=PC2,color=treatment,shape=as.factor(sampling)),data=col.inv.env3,size=3,alpha=.5,position=position_dodge(0.01))+
    # geom_segment(aes(x=start.x,y=start.y,xend=mid.x,yend=mid.y,color=treatment),data=col.inv.segments)+
    # geom_segment(aes(x=mid.x,y=mid.y,xend=end.x,yend=end.y,color=treatment),data=col.inv.segments,  arrow = arrow(length = unit(0.3, "cm"))) +
    geom_point(aes(x=PC1,y=PC2,color=treatment),data=col.inv.env3%>%filter(sampling!=1),size=3,alpha=.5)+
    geom_segment(aes(x=start.x,y=start.y,xend=end.x,yend=end.y,color=treatment),data=col.inv.segments,arrow = arrow(length = unit(0.3, "cm")))+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_viridis_d(option="A", begin=0, end=0.6,name="Treatment",labels=c("Control","Structure","Sponge"))+
    scale_shape_discrete(name = "Sampling"))


fish.plot / inv.plot / col.inv.plot +plot_layout(guides="collect")


plot(i.pca)


