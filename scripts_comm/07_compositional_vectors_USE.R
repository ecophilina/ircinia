# updated attempt at compositional vector analysis - using appropriate circular statistics
# for angles

#load packages and data----
library(vegan)
library(rphylopic)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(FSA)
if(!require(circular))install.packages("circular");library(circular)
if(!require(bpnreg))install.packages("bpnreg");library(bpnreg)
if(!require(gridGraphics))install.packages("gridGraphics");library(gridGraphics)
source("scripts_comm/02_community_data_org.R")

# load animal shapes----
fishpng <- image_data("0b9cdf1f-ccbc-4922-8cf6-60f90d07107e", size = 256)[[1]]
crabpng <- image_data("9958579e-5e63-4b7c-8e76-9b1a92d7f7ca", size = 256)[[1]]
clonalpng <- image_data("dbbf1325-10e5-4880-a27b-2d9afb5dc55c", size = 256)[[1]]

# organize fish data ----
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
         dx=ifelse(abs(dx)<0.0000000001,0,dx),
         dy=ifelse(abs(dy)<0.0000000001,0,dy),
         vlength = sqrt(dx^2 +dy^2),
         vlength.rsc = scales::rescale(vlength,to = c(0,1),from = range(vlength)),
         angle = atan2(dy,dx),
         angle.d=angle*180/pi,
         angle.d=ifelse(angle.d<0,360+angle.d,angle.d),
         plot.start = 0,
         plot.end.x = cos(angle)*vlength.rsc,
         plot.end.y = sin(angle)*vlength.rsc,
         plot.end.x2 = end.x - start.x,
         plot.end.y2 = end.y - start.y,
         ang.circ=circular(angle,units="radians"),
         ang.circ.d=circular(angle.d,units="degrees"),
         lsize = ifelse(duplicated(vlength)&duplicated(angle),2,1))

fish.env4$treatment <- as.factor(fish.env4$treatment)

# organize invert data ----
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
         dx=ifelse(abs(dx)<0.0000000001,0,dx),
         dy=ifelse(abs(dy)<0.0000000001,0,dy),
         vlength = sqrt(dx^2 +dy^2),
         vlength.rsc = scales::rescale(vlength,to = c(0,1),from = range(vlength)),
         angle = atan2(dy,dx),
         angle.d=angle*(180/pi),
         angle.d=ifelse(angle.d<0,360+angle.d,angle.d),
         plot.start = 0,
         plot.end.x = cos(angle)*vlength.rsc,
         plot.end.y = sin(angle)*vlength.rsc,
         plot.end.x2 = end.x - start.x,
         plot.end.y2 = end.y - start.y,
         ang.circ=circular(angle,units="radians"),
         ang.circ.d=circular(angle.d,units="degrees"),
         lsize = ifelse(duplicated(vlength)&duplicated(angle),2,1))

inv.env4$treatment <- as.factor(inv.env4$treatment)

#summary tables and visualization----
# set global ggplot options
theme_set(theme_bw()+theme(panel.grid=element_blank()))
# fish ----
#summary table
fish.vl<-fish.env4%>%
  group_by(treatment)%>%
  summarize(vlength.mean=signif(mean(vlength),2),
            vlength.sd=signif(sd(vlength),2))
fish.ang<-fish.env4%>%
  filter(vlength>0)%>%
  group_by(treatment)%>%
  summarize(n.longerthan.0=n(),
            angle.mean=signif(mean(ang.circ.d),2),
            angle.rho=signif(rho.circular(ang.circ.d),2),
            angle.sd=signif(sd(ang.circ.d),2))

(fish.summary<-left_join(fish.vl,fish.ang))

#visualize vector lengths
fish.vlengths<-ggplot(data=fish.env4)+
  geom_jitter(aes(x=treatment,y=vlength),size=4,alpha=.3,width=.1)

# joint plot
fish.1<-fish.vlengths+~plot.circular(fish.env4$ang.circ.d[fish.env4$treatment=="real"&fish.env4$vlength>0],main="live sponge")
fish.1<-fish.1+~plot.circular(fish.env4$ang.circ.d[fish.env4$treatment=="fake"&fish.env4$vlength>0],main="structure")
fish.1+~plot.circular(fish.env4$ang.circ.d[fish.env4$treatment=="blank"&fish.env4$vlength>0],main="control")

# inverts ----

inv.vl<-inv.env4%>%
  group_by(treatment)%>%
  summarize(vlength.mean=signif(mean(vlength),2),
            vlength.sd=signif(sd(vlength),2))
inv.ang<-inv.env4%>%
  filter(vlength>0)%>%
  group_by(treatment)%>%
  summarize(n.longerthan.0=n(),
            angle.mean=signif(mean(ang.circ.d),2),
            angle.rho=signif(rho.circular(ang.circ.d),2),
            angle.sd=signif(sd(ang.circ.d),2))

(inv.summary<-left_join(inv.vl,inv.ang))

#visualize vector lengths
inv.vlengths<-ggplot(data=inv.env4)+
  geom_jitter(aes(x=treatment,y=vlength),size=4,alpha=.3,width=.1)

# joint plot
inv.1<-inv.vlengths|~plot.circular(inv.env4$ang.circ.d[inv.env4$treatment=="real"&inv.env4$vlength>0],main="live sponge")
inv.1<-inv.1+~plot.circular(inv.env4$ang.circ.d[inv.env4$treatment=="fake"&inv.env4$vlength>0],main="structure")
inv.1+~plot.circular(inv.env4$ang.circ.d[inv.env4$treatment=="blank"&inv.env4$vlength>0],main="control")

# statistics----
#fish
#vector length
(fish.vl.kw<-kruskal.test(vlength~treatment,data=fish.env4))
# significant so look at differences between groups
dunnTest(vlength~treatment,data=fish.env4)
# live sponge vectors are longer than control and structure plots
# no difference between control and structure plots

# vector angle
# insufficient data points to complete this analysis

# inverts
#vector length
(inv.vl.kw<-kruskal.test(vlength~treatment,data=inv.env4))
# no significant difference between vector lengths but we can look at angle here
# compare structure and control plots
watson.two.test(inv.env4$ang.circ.d[inv.env4$treatment=="fake"&inv.env4$vlength>0],
                inv.env4$ang.circ.d[inv.env4$treatment=="blank"&inv.env4$vlength>0])

# compare structure and live sponge plots
watson.two.test(inv.env4$ang.circ.d[inv.env4$treatment=="real"&inv.env4$vlength>0],
                inv.env4$ang.circ.d[inv.env4$treatment=="fake"&inv.env4$vlength>0])

# compare control and live sponge plots
watson.two.test(inv.env4$ang.circ.d[inv.env4$treatment=="real"&inv.env4$vlength>0],
                inv.env4$ang.circ.d[inv.env4$treatment=="blank"&inv.env4$vlength>0])

#try the watson williams test
watson.williams.test(ang.circ.d~treatment,inv.env4)

trialx<-circular(c(0,0,0,0,0,0,0,0,0,0),units="degrees")
trialy<-circular(c(180,180,180,180,180,180,180,180,180,180),units="degrees")
plot(trialx)
plot(trialy)
watson.two.test(trialx,trialy)
