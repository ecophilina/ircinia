# first attempt at compositional vector analysis for fish

library(vegan)
library(rphylopic)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggpubr)

source("scripts_comm/02_community_data_org.R")

# load animal shapes
fishpng <- image_data("0b9cdf1f-ccbc-4922-8cf6-60f90d07107e", size = 256)[[1]]
crabpng <- image_data("9958579e-5e63-4b7c-8e76-9b1a92d7f7ca", size = 256)[[1]]
clonalpng <- image_data("dbbf1325-10e5-4880-a27b-2d9afb5dc55c", size = 256)[[1]]


# Fish 

# subset down to only 0 and 12 month sampling data to start
fish.com2<-fish.com[fish.env$sampling %in% c(0,1,12),]
fish.env2<-fish.env[fish.env$sampling %in% c(0,1,12),]

# do hellinger transformation then RDA
# normally RDA is used for “constrained ordination” (ordination w/covariates or predictor)
# without predictors, RDA is the same as PCA
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

fish.env4$treatment <- as.factor(fish.env4$treatment)


# Look at just differences in vector length between treatments
fish.vl.aov<-aov(vlength~treatment,data=fish.env4%>%
                  mutate(treatment = relevel(treatment, ref = "real")))
summary(fish.vl.aov)
TukeyHSD(fish.vl.aov)

#Checking residuals 
par(mfrow=c(2,2))
plot(fish.vl.aov)

#checking the distributions
par(mfrow=c(1,1))
hist(fish.env4$vlength)

#Look at just differences in vector angle between treatments
fish.angle.aov<-aov(angle~treatment,data=fish.env4%>%
                     mutate(treatment = relevel(treatment, ref = "real")))
summary(fish.angle.aov)
TukeyHSD(fish.angle.aov)

#Checking residuals 
par(mfrow=c(2,2))
plot(fish.angle.aov)

#checking the distributions
par(mfrow=c(1,1))
hist(fish.env4$angle)


#to reference in manuscript
(fish.vl.aov.sum<-anova(fish.vl.aov))
(fish.vl.aov.tuk<-TukeyHSD(fish.vl.aov)$treatment)
  
(fish.angle.aov.sum<-anova(fish.angle.aov))
(fish.angle.aov.tuk<-TukeyHSD(fish.angle.aov)$treatment)

write_rds(fish.vl.aov.sum,"working_data/FishVLSum.rds")
write_rds(fish.vl.aov.tuk,"working_data/FishVLTuk.rds")
write_rds(fish.angle.aov.sum,"working_data/FishAngleSum.rds")
write_rds(fish.angle.aov.tuk,"working_data/FishAngleTuk.rds")


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
    scale_color_viridis_d(option="A", begin=0, end=0.65,name="Treatment",labels=c("Control","Structure","Sponge"))+
    scale_size_discrete(range=c(1,1.75),name = "Number of Plots")+
    scale_alpha_continuous(range=c(0.75,1),guide=FALSE)+
    #  geom_text(aes(x=-1,y=1),label="a",size=10)+
    add_phylopic(fishpng,x=-.8,y=.9,ysize=.25,alpha=1))

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
    scale_color_viridis_d(option="A", begin=0, end=0.65,name="Treatment",labels=c("Control","Structure","Sponge"))+
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

inv.env4$treatment <- as.factor(inv.env4$treatment)

# Look at differences in vector length between treatments
inv.vl.aov<-aov(vlength~treatment,data=inv.env4%>%
                 mutate(treatment = relevel(treatment, ref = "real")))
summary(inv.vl.aov)
# TukeyHSD(inv.vl.aov)# not sig above so not needed

#Checking residuals 
par(mfrow=c(2,2))
plot(inv.vl.aov)

#checking the distributions
par(mfrow=c(1,1))
hist(inv.env4$vlength)

#Look at differences in vector length between treatments
inv.angle.aov<-aov(angle~treatment,data=inv.env4%>%
                    mutate(treatment = relevel(treatment, ref = "real")))
summary(inv.angle.aov)
# TukeyHSD(inv.angle.aov) # not sig above so not needed

#Checking residuals 
par(mfrow=c(2,2))
plot(inv.angle.aov)

#checking the distributions
par(mfrow=c(1,1))
hist(inv.env4$angle)


#to reference in manuscript
(inv.vl.aov.sum<-anova(inv.vl.aov))

(inv.angle.aov.sum<-anova(inv.angle.aov))

write_rds(inv.vl.aov.sum,"working_data/InvVLSum.rds")
write_rds(inv.angle.aov.sum,"working_data/InvAngleSum.rds")



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
    scale_color_viridis_d(option="A", begin=0, end=0.65)+
    #  geom_text(aes(x=-1,y=1),label="b",size=10)+
    add_phylopic(crabpng,x=-.8,y=.9,ysize=.45,alpha=1))

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
    scale_color_viridis_d(option="A", begin=0, end=0.65,name="Treatment",labels=c("Control","Structure","Sponge"))+
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

col.inv.env4$treatment <- as.factor(col.inv.env4$treatment)


# Look at differences in vector length between treatments
col.inv.vl.aov<-aov(vlength~treatment,data=col.inv.env4%>%
                     mutate(treatment = relevel(treatment, ref = "real")))
summary(col.inv.vl.aov)
# TukeyHSD(col.inv.vl.aov) # not sig above so not needed

#Checking residuals 
par(mfrow=c(2,2))
plot(col.inv.vl.aov)

#checking the distributions
par(mfrow=c(1,1))
hist(col.inv.env4$vlength)

#Look at differences in vector angle between treatments
col.inv.angle.aov<-aov(angle~treatment,data=col.inv.env4%>%
                        mutate(treatment = relevel(treatment, ref = "real")))
summary(col.inv.angle.aov)
# TukeyHSD(col.inv.angle.aov) # not sig above so not needed

#Checking residuals 
par(mfrow=c(2,2))
plot(col.inv.angle.aov)

#checking the distributions
par(mfrow=c(1,1))
hist(col.inv.env4$angle)


#to reference in manuscript
(col.inv.vl.aov.sum<-anova(col.inv.vl.aov))

(col.inv.angle.aov.sum<-anova(col.inv.angle.aov))

write_rds(col.inv.vl.aov.sum,"working_data/ColInvVLSum.rds")
write_rds(col.inv.angle.aov.sum,"working_data/ColInvAngleSum.rds")
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
    scale_color_viridis_d(option="A", begin=0, end=0.65)+
    scale_size_continuous(range=c(1,1.75))+
    scale_alpha_continuous(range=c(0.75,1))+
    #    geom_text(aes(x=-1,y=1),label="c",size=10)+
    add_phylopic(clonalpng,x=-.8,y=.9,ysize=.25,alpha=1))

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
    scale_color_viridis_d(option="A", begin=0, end=0.65,name="Treatment",labels=c("Control","Structure","Sponge"))+
    scale_shape_discrete(name = "Sampling"))



## THIS IS THE PLOT FOR THE PAPER

fish.vplot / inv.plot  + plot_layout(guides="collect")

ggsave("figures/community_vector_plots.jpg",dpi=300,width=4,height=7)
ggsave("figures/community_vector_plots.png",dpi=300,width=4,height=7)
ggsave("figures/community_vector_plots.pdf",dpi=300,width=4,height=7)

# now making species association plots

#fish
fsp<-scores(f.pca,choices=c(1,2),display = "species")
fsp<-data.frame(fsp)%>%
  filter((abs(PC1)+abs(PC2))>.01)%>%
  mutate(sp=row.names(.))

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
# the denominator of the diameter is the number of principle components where the
# eigenvalue is >0
circ <- circleFun(center = c(0, 0), diameter = sqrt(2 / 7), npoints = 500)

# make one plot to see who to keep
ggplot(data=fsp)+
  geom_segment(aes(x=0,xend=PC1,y=0,yend=PC2,group=sp,color=sp),
               arrow = arrow(length = unit(0.3, "cm")))+
  geom_path(data = circ, aes(x, y), lty = 2, color = "grey", alpha = 0.7)

# keep grunt, damselfish, slippery dick, and green blotch and rescale to vector length 1



fsp2<-filter(fsp,sp %in% c("grunt", "damselfish", "slippery dick", "green blotch"))%>%
  mutate(vlength = sqrt(PC1^2 +PC2^2),
         vlength.rsc = vlength/max(vlength),
         angle = atan2(PC2,PC1),
         plot.start = 0,
         plot.end.x = cos(angle)*vlength.rsc,
         plot.end.y = sin(angle)*vlength.rsc,)

#now make pretty graph
(fish.spp<-ggplot(data=fsp2)+
    # geom_vline(aes(xintercept=0),linetype="dashed",alpha=.5)+
    # geom_hline(aes(yintercept=0),linetype="dashed",alpha=.5)+
    geom_segment(aes(x=0,xend=PC1,y=0,yend=PC2),
                 arrow = arrow(length = unit(0.3, "cm")))+
    coord_fixed(xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))+
    #geom_path(data = circ, aes(x, y), lty = 2, alpha = 0.5)+
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())+
    # ylab("PC2")+
    # xlab("PC1")+
    geom_text(aes(x=PC1,y=PC2,
                  label=c("Slippery dick","Green blotch parrotfish","Damselfish","Grunts")),
              nudge_x = c(.4,.2,.6,-.3),
              nudge_y = c(-.15,-.1,.1,.1)))#+
#add_phylopic(fishpng,x=-.8,y=.9,ysize=.25,alpha=1))


#inverts
isp<-scores(i.pca,choices=c(1,2),display = "species")
isp<-data.frame(isp)%>%
  filter((abs(PC1)+abs(PC2))>.1)%>%
  mutate(sp=row.names(.))

# the denominator of the diameter is the number of principle components where the
# eigenvalue is >0
circ <- circleFun(center = c(0, 0), diameter = sqrt(2 / 8), npoints = 500)

# make one plot to see who to keep
ggplot(data=isp)+
  geom_segment(aes(x=0,xend=PC1,y=0,yend=PC2,group=sp,color=sp),
               arrow = arrow(length = unit(0.3, "cm")))+
  geom_path(data = circ, aes(x, y), lty = 2, color = "grey", alpha = 0.7)

# keep cerith, little white snail, and ground tunicate


isp2<-filter(isp,sp %in% c("cerith", "little white snail", "ground tunicate"))%>%
  mutate(vlength = sqrt(PC1^2 +PC2^2),
         vlength.rsc = vlength/max(vlength),
         angle = atan2(PC2,PC1),
         plot.start = 0,
         plot.end.x = cos(angle)*vlength.rsc,
         plot.end.y = sin(angle)*vlength.rsc,)

#now make pretty graph
(invrt.spp<-ggplot(data=isp2)+
    # geom_vline(aes(xintercept=0),linetype="dashed",alpha=.5)+
    # geom_hline(aes(yintercept=0),linetype="dashed",alpha=.5)+
    geom_segment(aes(x=0,xend=PC1,y=0,yend=PC2),
                 arrow = arrow(length = unit(0.3, "cm")))+
    coord_fixed(xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))+
    #geom_path(data = circ, aes(x, y), lty = 2, alpha = 0.5)+
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())+
    # ylab("PC2")+
    # xlab("PC1")+
    geom_text(aes(x=PC1,y=PC2,
                  label=c("Ceriths","Whelk","Sea squirt")),
              nudge_y = c(.1,-.1,-.1)))#+
#add_phylopic(crabpng,x=-.8,y=.9,ysize=.25,alpha=1))

# #colonial inverts
# cisp<-scores(ci.pca,choices=c(1,2),display = "species")
# cisp<-data.frame(cisp)%>%
#   mutate(sp=row.names(.))
# 
# # the denominator of the diameter is the number of principle components where the
# # eigenvalue is >0
# circ <- circleFun(center = c(0, 0), diameter = sqrt(2 / 7), npoints = 500)
# 
# # make one plot to see who to keep
# ggplot(data=cisp)+
#   geom_segment(aes(x=0,xend=PC1,y=0,yend=PC2,group=sp,color=sp),
#                arrow = arrow(length = unit(0.3, "cm")))+
#   geom_path(data = circ, aes(x, y), lty = 2, color = "grey", alpha = 0.7)
# 
# # keep cerith, little white snail, and ground tunicate


# cisp2<-filter(cisp,sp %in% c("black and orange tunicate", "pinkish brown tunicate", "white and green tunicate"))%>%
#   mutate(vlength = sqrt(PC1^2 +PC2^2),
#          vlength.rsc = vlength/max(vlength),
#          angle = atan2(PC2,PC1),
#          plot.start = 0,
#          plot.end.x = cos(angle)*vlength.rsc,
#          plot.end.y = sin(angle)*vlength.rsc,)
# 
# #now make pretty graph
# (cinvrt.spp<-ggplot(data=cisp2)+
#     # geom_vline(aes(xintercept=0),linetype="dashed",alpha=.5)+
#     # geom_hline(aes(yintercept=0),linetype="dashed",alpha=.5)+
#     geom_segment(aes(x=0,xend=PC1,y=0,yend=PC2),
#                  arrow = arrow(length = unit(0.3, "cm")))+
#     coord_fixed(xlim=c(-1.5,1.5),ylim=c(-1.5,1.5))+
#     #geom_path(data = circ, aes(x, y), lty = 2, alpha = 0.5)+
#     theme(panel.grid = element_blank(),
#           panel.background = element_blank(),
#           axis.text = element_blank(),
#           axis.ticks = element_blank(),
#           axis.title = element_blank())+
#     # ylab("PC2")+
#     # xlab("PC1")+
#     geom_text(aes(x=PC1,y=PC2,
#                   label=c("Tunicate-white","Tunicate-pink","Tunciate-black + orange")),
#               nudge_x = c(-.57,0,.1),
#               nudge_y = c(.18,.15,-.15)))#+
# #add_phylopic(clonalpng,x=-.8,y=.9,ysize=.25,alpha=1))
# 

(fish.vplot +fish.spp) / (inv.plot +invrt.spp) + plot_layout(guides="collect")+plot_annotation(tag_levels="a")

ggsave("figures/fig4.tiff",dpi=600,width=174,height=130,units = c("mm"))


# Plot vector traits where amount of change (vector length) vs. direction of change (vector angle)
# all sponge plots change in much more similar ways for fish and non-clonal inverts
# comment on how for fish sponge plots change more than many of the control and structure plots


(by1 <- ggplot(fish.env4, aes(angle, vlength, colour=treatment)) + 
  geom_jitter(size = 4, alpha = 0.6, width = 0.4, height = 0.0) + 
  ylab(" ") + xlab(" ") +
  scale_color_viridis_d(option="A", begin=0, end=0.65,name="Treatment",labels=c("Control","Structure","Sponge"))+
  ggsidekick::theme_sleek()+theme(axis.title.x = element_blank()) +
  add_phylopic(fishpng,x=-2.4,y=.8,ysize=.45,alpha=1))

(by2 <- ggplot(inv.env4, aes(angle, vlength, colour=treatment)) + 
  geom_jitter(size = 4, alpha = 0.6, width = 0.1, height = 0.0) + 
  ylab("Vector length") + xlab(" ") +
  scale_color_viridis_d(option="A", begin=0, end=0.65,name="Treatment",labels=c("Control","Structure","Sponge"))+
  ggsidekick::theme_sleek()+ theme(axis.title.x = element_blank()) +
  add_phylopic(crabpng,x=-1.9,y=.95,ysize=.75,alpha=1))

# (by3 <- ggplot(col.inv.env4, aes(angle, vlength, colour=treatment)) + 
#   geom_jitter(size = 4, alpha = 0.6, width = 0.2, height = 0.0) + 
#   ylab(" ") + 
#   xlab("Vector angle") +
#   scale_color_viridis_d(option="A", begin=0, end=0.65,name="Treatment",labels=c("Control","Structure","Sponge"))+
#   ggsidekick::theme_sleek()+
#   add_phylopic(clonalpng,x=-2.6,y=1.2,ysize=.65,alpha=1))


(by1 / by2 )+ plot_layout(guides="collect") + plot_annotation(tag_levels = 'a')

# ggsave("figures/upload_vect-length-by-angle-plot.tiff", dpi=300, width = 4, height = 8)

