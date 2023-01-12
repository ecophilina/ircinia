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
if(!require(effectsize))install.packages("effectsize");library(effectsize)
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
write.csv(fish.summary,"working_data/fish_compvector_summary.csv",row.names = FALSE)

#visualize vector lengths
(fish.vlengths<-ggplot(data=fish.env4)+
    #geom_boxplot(aes(x=treatment,y=vlength))+
    geom_jitter(aes(x=treatment,y=vlength),size=4,alpha=.3,width=.15)+
    ylab("Vector length")+
    xlab("")+
    scale_x_discrete(labels=c("Control","Structure","Live sponge"))+
    theme(axis.text = element_text(size=12),
          axis.title.y=element_text(size=14)))

# joint plot
fish.1<-fish.vlengths+~plot.circular(fish.env4$ang.circ.d[fish.env4$treatment=="real"&fish.env4$vlength>0],main="live sponge")
fish.1<-fish.1+~plot.circular(fish.env4$ang.circ.d[fish.env4$treatment=="fake"&fish.env4$vlength>0],main="structure")
fish.1+~plot.circular(fish.env4$ang.circ.d[fish.env4$treatment=="blank"&fish.env4$vlength>0],main="control")
# save individual plots to make multipanel figure in gimp
png(filename="fish_livesponge_angleplot.png",width=500,height=500, units="px")
plot.circular(fish.env4$ang.circ.d[fish.env4$treatment=="real"&fish.env4$vlength>0],main="Live sponge",cex=2,stack=TRUE,sep=.05,bins=360)
dev.off()

png(filename="fish_structure_angleplot.png",width=500,height=500, units="px")
plot.circular(fish.env4$ang.circ.d[fish.env4$treatment=="fake"&fish.env4$vlength>0],main="Structure",cex=2,stack=TRUE,sep=.05,bins=360)
dev.off()

png(filename="fish_control_angleplot.png",width=500,height=500, units="px")
plot.circular(fish.env4$ang.circ.d[fish.env4$treatment=="blank"&fish.env4$vlength>0],main="Control",cex=2,stack=TRUE,sep=.05,bins=360)
dev.off()

ggsave(filename = "fish_vectorlengthplot.jpg",fish.vlengths,width=5,height=4,dpi=300)
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
write.csv(inv.summary,"working_data/inv_compvector_summary.csv",row.names = FALSE)

#visualize vector lengths
(inv.vlengths<-ggplot(data=inv.env4)+
    #geom_boxplot(aes(x=treatment,y=vlength))+
    geom_jitter(aes(x=treatment,y=vlength),size=4,alpha=.3,width=.15)+
    ylab("Vector length")+
    xlab("")+
    scale_x_discrete(labels=c("Control","Structure","Live sponge"))+
    theme(axis.text = element_text(size=12),
          axis.title.y=element_text(size=14)))

# joint plot
inv.1<-inv.vlengths|~plot.circular(inv.env4$ang.circ.d[inv.env4$treatment=="real"&inv.env4$vlength>0],main="live sponge")
inv.1<-inv.1+~plot.circular(inv.env4$ang.circ.d[inv.env4$treatment=="fake"&inv.env4$vlength>0],main="structure")
inv.1+~plot.circular(inv.env4$ang.circ.d[inv.env4$treatment=="blank"&inv.env4$vlength>0],main="control")

# save individual plots to make multipanel figure in gimp
png(filename="inv_livesponge_angleplot.png",width=500,height=500, units="px")
plot.circular(inv.env4$ang.circ.d[inv.env4$treatment=="real"&inv.env4$vlength>0],main="Live sponge",cex=2,stack=TRUE,sep=.05,bins=360)
dev.off()

png(filename="inv_structure_angleplot.png",width=500,height=500, units="px")
plot.circular(inv.env4$ang.circ.d[inv.env4$treatment=="fake"&inv.env4$vlength>0],main="Structure",cex=2,stack=TRUE,sep=.05,bins=120)
dev.off()

png(filename="inv_control_angleplot.png",width=500,height=500, units="px")
plot.circular(inv.env4$ang.circ.d[inv.env4$treatment=="blank"&inv.env4$vlength>0],main="Control",cex=2,stack=TRUE,sep=.05,bins=360)
dev.off()

ggsave(filename = "inv_vectorlengthplot.jpg",inv.vlengths,width=5,height=4,dpi=300)


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

# doesn't seem like its actually possible to detect differences with this test and our sample size
# making a trial data set
trialx<-circular(c(0,0,0,0,0),units="degrees")
trialy<-circular(c(180,180,180,180,180),units="degrees")

watson.two.test(trialx,trialy)

# according to the tables I've found, nope not really possible.

# try a manova 

inv.ang<-inv.env4%>%
  filter(vlength>0)

inv.man<-manova(cbind(sin(angle),cos(angle))~treatment,data=inv.ang)

summary(inv.man)# not quite significant but if you look at effect sizes it is large
eta_squared(inv.man) #effect size - from what I've read 0.14 is considered large, but the confidence interval reaches 0


# try a manova where we compare structure plots to live sponge - from looking at 
# the data control plots to do not move in a consistent direction (sd is more than twice structure and live sponge treatments)

inv.ang2<-filter(inv.ang,treatment!="blank")
inv.man2<-manova(cbind(sin(angle),cos(angle))~treatment,data=inv.ang2)

summary(inv.man2)# significant
eta_squared(inv.man2, ci=.95) #effect size - from what I've read 0.14 is considered large
# for this comparison the effect size is large and the confidence intervals are 
# above the "large" cutoff
