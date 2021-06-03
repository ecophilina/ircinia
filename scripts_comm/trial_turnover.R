# install.packages("codyn")
#remotes::install_github("sckott/rphylopic")

library(codyn)
library(rphylopic)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggpubr)

source("scripts_comm/02_community_data_org.R")

#invertebrates

inv.turn<-inv.com.full%>%
  mutate(dummy=1)%>%
  pivot_longer(-treatment:-a.abund.c,names_to="taxa",values_to="abundance")%>%
  filter(!sampling %in% c(5,17))

inv.appear <- turnover(df=inv.turn,
                      time.var = "sampling",
                      species.var = "taxa",
                      abundance.var = "abundance",
                      replicate.var = "plot",
                      metric="appearance")%>%
  mutate(plot=as.numeric(plot))%>%
  left_join(inv.env)


inv.dis <- turnover(df=inv.turn,
                       time.var = "sampling",
                       species.var = "taxa",
                       abundance.var = "abundance",
                       replicate.var = "plot",
                       metric="disappearance")%>%
  mutate(plot=as.numeric(plot))%>%
  left_join(inv.appear)

inv.turnover <- turnover(df=inv.turn,
                         time.var = "sampling",
                         species.var = "taxa",
                         abundance.var = "abundance",
                         replicate.var = "plot")%>%
  mutate(plot=as.numeric(plot))%>%
  left_join(inv.dis)%>%
  left_join(inv.uni)


ggplot(inv.turnover)+
#  geom_line(aes(x=sampling, y=appearance,color=treatment,group=plot))
#  geom_line(aes(x=sampling, y=disappearance,color=treatment,group=plot))
#  geom_line(aes(x=sampling, y=total,color=treatment,group=plot))
  geom_point(aes(y=appearance,x=disappearance,color=treatment,size=spr))+
  facet_wrap(~sampling)


# fish

fish.turn<-fish.com.full%>%
  mutate(dummy=1)%>%
  pivot_longer(-treatment:-a.abund.c,names_to="taxa",values_to="abundance")%>%
  filter(!sampling %in% c(5,17))

fish.appear <- turnover(df=fish.turn,
                       time.var = "sampling",
                       species.var = "taxa",
                       abundance.var = "abundance",
                       replicate.var = "plot",
                       metric="appearance")%>%
  mutate(plot=as.numeric(plot))%>%
  left_join(fish.env)


fish.dis <- turnover(df=fish.turn,
                    time.var = "sampling",
                    species.var = "taxa",
                    abundance.var = "abundance",
                    replicate.var = "plot",
                    metric="disappearance")%>%
  mutate(plot=as.numeric(plot))%>%
  left_join(fish.appear)

fish.turnover <- turnover(df=fish.turn,
                         time.var = "sampling",
                         species.var = "taxa",
                         abundance.var = "abundance",
                         replicate.var = "plot")%>%
  mutate(plot=as.numeric(plot))%>%
  left_join(fish.dis)%>%
  left_join(fish.uni)

#colonial inverts
col.inv.turn<-col.inv.com.full%>%
  mutate(dummy=1)%>%
  pivot_longer(-treatment:-a.abund.c,names_to="taxa",values_to="abundance")%>%
  filter(!sampling %in% c(5,17))

col.inv.appear <- turnover(df=col.inv.turn,
                       time.var = "sampling",
                       species.var = "taxa",
                       abundance.var = "abundance",
                       replicate.var = "plot",
                       metric="appearance")%>%
  mutate(plot=as.numeric(plot))%>%
  left_join(col.inv.env)


col.inv.dis <- turnover(df=col.inv.turn,
                    time.var = "sampling",
                    species.var = "taxa",
                    abundance.var = "abundance",
                    replicate.var = "plot",
                    metric="disappearance")%>%
  mutate(plot=as.numeric(plot))%>%
  left_join(col.inv.appear)

col.inv.turnover <- turnover(df=col.inv.turn,
                         time.var = "sampling",
                         species.var = "taxa",
                         abundance.var = "abundance",
                         replicate.var = "plot")%>%
  mutate(plot=as.numeric(plot))%>%
  left_join(col.inv.dis)%>%
  left_join(col.inv.uni)


# exploratory plots for non-clonal inverts
ggplot(inv.turnover)+
  geom_line(aes(x=sampling, y=appearance,color=treatment,group=plot), alpha = 0.5)
ggplot(inv.turnover)+
  geom_line(aes(x=sampling, y=disappearance,color=treatment,group=plot), alpha = 0.5)
# geom_line(aes(x=sampling, y=total,color=treatment,group=plot))

ggplot(inv.turnover)+
  geom_jitter(aes(y=appearance,x=disappearance,color=treatment,size=spr))+
  facet_wrap(~sampling)

# exploratory plots for clonal inverts
ggplot(col.inv.turnover)+
  geom_line(aes(x=sampling, y=appearance,color=treatment,group=plot), alpha = 0.5)
ggplot(col.inv.turnover)+
  geom_line(aes(x=sampling, y=disappearance,color=treatment,group=plot), alpha = 0.5)
  # geom_line(aes(x=sampling, y=total,color=treatment,group=plot))

ggplot(col.inv.turnover)+
  geom_jitter(aes(y=appearance,x=disappearance,color=treatment,size=spr))+
  facet_wrap(~sampling)


# find animal images
# http://phylopic.org/image/browse/
# try Hyperprosopon argenteum for fish (http://phylopic.org/image/0b9cdf1f-ccbc-4922-8cf6-60f90d07107e/) and blue crab for inverts
fishpng <- image_data("0b9cdf1f-ccbc-4922-8cf6-60f90d07107e", size = 128)[[1]]
crabpng <- image_data("9958579e-5e63-4b7c-8e76-9b1a92d7f7ca", size = 256)[[1]]

# tunicate options
# Tunicata: http://phylopic.org/image/dbbf1325-10e5-4880-a27b-2d9afb5dc55c/
# Ciona intestinalis: http://phylopic.org/name/534a6a2a-1bb7-4163-ab02-e2e69b6d045a
# Ciona savignyi: http://phylopic.org/name/5288b025-0974-4bbd-8dd1-394748322559

# # alternate fish: Prognathodes sp
# fishpng <- image_data("8051e46b-c0b3-4e48-8640-3c84f105f107", size = 128)[[1]]

clonalpng <- image_data("dbbf1325-10e5-4880-a27b-2d9afb5dc55c", size = 256)[[1]]

# or use an image saved to our project folder 
# library(png)
# fishpng <- readPNG("scripts_comm/Hyperprosopon_argenteum.png")

fish.turnover$group <- "fish"
inv.turnover$group <- "invert"
col.inv.turnover$group <- "clonal"  

fish.turnover$spr
inv.turnover$spr
col.inv.turnover$spr

plot_png <- function(dat, 
                     png_list, # a list() with same num of images as groups, in same order as unique(dat$group)
                     png_group = "group", # variable that determines which image
                     x = "disappearance", 
                     y = "appearance",
                     xlim = c(-0.03, 0.8),
                     ylim = c(-0.03, 0.8),
                     alpha = 0.8,
                     scal_fac = c(100) # vector containing scaling factors for each image 
  ) {
  dat$x <- dat[[x]]
  dat$y <- dat[[y]]
  dat$group <- dat[[png_group]]

  p <- ggplot(dat, aes(x, y)) +
    geom_blank() + 
    coord_fixed(xlim = xlim, ylim = ylim) +
    xlab(x) + ylab(y) +
    ggsidekick::theme_sleek(base_size = 18)
  
# loop over groups to get different png images
for (i in 1:length(unique(dat$group))) {
  d <- filter(dat, treatment == "real" & group == unique(dat$group)[i]) 
  for (j in 1:nrow(d)) {
    p <- p + add_phylopic(png_list[[i]], alpha,
      d$x[j]+runif(1, -0.02, 0.02), d$y[j]+runif(1, -0.02, 0.02),
      col = "#E95562FF",
      ysize = (d$spr[j] + 2) / scal_fac[i]
    )
  }
  d <- filter(dat, treatment == "fake" & group == unique(dat$group)[i])
  for (j in 1:nrow(d)) {
    p <- p + add_phylopic(png_list[[i]], alpha,
      d$x[j]+runif(1, -0.02, 0.02), d$y[j]+runif(1, -0.02, 0.02),
      col = "#5A167EFF",
      ysize = (d$spr[j] + 2) / scal_fac[i]
    )
  }
  d <- filter(dat, treatment == "blank" & group == unique(dat$group)[i])
  for (j in 1:nrow(d)) {
    p <- p + add_phylopic(png_list[[i]], alpha,
      d$x[j]+runif(1, -0.02, 0.02), d$y[j]+runif(1, -0.02, 0.02),
      col = "black",
      ysize = (d$spr[j] + 2) / scal_fac[i]
    )
  }
}
  p 
}


# I tried adding tunicates but it just confused things
turnoverdat <- bind_rows(inv.turnover, fish.turnover) 
unique(turnoverdat$group)
png_list <- list(crabpng, fishpng)


# each jitter is random, so just rerun plot code if points are falling on edge of plot area
(p1 <- plot_png(filter(turnoverdat, sampling == 1), png_list, scal_fac = c(70, 100)) +
    ylab("Proportion of Species Gained") + 
    xlab("Proportion of Species Lost") +
    theme(plot.title=element_text(hjust=0.5))+
    ggtitle("1 month"))

(p12 <- plot_png(filter(turnoverdat, sampling == 12), png_list, scal_fac = c(70, 100)) + 
    ylab("") + 
    xlab("Proportion of Species Lost") +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank()) + 
    theme(plot.title=element_text(hjust=0.5))+
    ggtitle("12 months"))

p1 + p12 + plot_layout(widths=c(1,1))

ggsave("turnover-plots.png", width = 10, height = 6)


# since tunicates show a different pattern, I'm putting them on different panels
(p1c <- plot_png(filter(col.inv.turnover, sampling == 1), png_list = list(clonalpng), scal_fac = c(80)) +
    ylab("Proportion of Species Gained") + xlab("Proportion of Species Lost")+ 
    ggtitle("1 month into experiment"))

(p12c <- plot_png(filter(col.inv.turnover, sampling == 12), png_list = list(clonalpng), scal_fac = c(80)) + 
    ylab("") + xlab("Proportion of Species Lost") + 
    ggtitle("12 months into experiment") + 
    theme(axis.text.y = element_blank()))

p1c + p12c + plot_layout(widths=c(1,1))
ggsave("turnover-plots-tunicates.png", width = 10, height = 6)

# (p1 + theme(axis.text.x = element_blank(), axis.title.x = element_blank())) + 
# (p12 + theme(axis.text.x = element_blank(), axis.title.x = element_blank())) + 
#   p1c + p12c + plot_layout(widths=c(1,1))
# 
# ggsave("turnover-plots-2x2.png", width = 10, height = 10)


# create a legend that could be added at some point
(p0 <- ggplot(turnoverdat, 
  aes(x = "disappearance", 
    y = "appearance", colour = treatment)) +
    geom_point() + 
    scale_color_viridis_d(option="A", begin=0, end=0.6)+
    ggsidekick::theme_sleek(base_size = 18)
)

leg <- ggpubr::get_legend(p0 + theme(legend.position = c(0.9,0.9)))

# make the species richness v time 

fspr.sum<-fish.uni%>%
  group_by(treatment,sampling)%>%
  summarize(spr.m=mean(spr),
            spr.sd=sd(spr))%>%
  filter(sampling %in% c(0,1,12))

ispr.sum<-inv.uni%>%
  group_by(treatment,sampling)%>%
  summarize(spr.m=mean(spr),
            spr.sd=sd(spr))%>%
  filter(sampling %in% c(0,1,12))

(fspr<-ggplot()+
  geom_point(data=fish.uni%>%filter(sampling %in% c(0,1,12)),aes(x=as.factor(sampling),y=spr,color=treatment),alpha=.3,position=position_dodge(0.3))+
  geom_point(data=fspr.sum,aes(x=as.factor(sampling),y=spr.m,color=treatment),alpha=.7,size=5,position=position_dodge(0.3))+
  geom_errorbar(data=fspr.sum,aes(x=as.factor(sampling),ymin=spr.m-spr.sd,ymax=spr.m+spr.sd,color=treatment),alpha=.7,width=.3,position=position_dodge(0.3))+ 
  scale_color_viridis_d(option="A", begin=0, end=0.6,name="Treatment",labels=c("Control","Structure","Sponge"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.margin=margin(t=5,r=1,b=1,l=5),
        axis.text.x = element_blank())+
  add_phylopic(fishpng,x=.5,y=4,ysize = 2,alpha=1)+
#  ggtitle("Fish")+
  ylab("")+
  xlab(""))

(ispr<-ggplot()+
    geom_point(data=inv.uni%>%filter(sampling %in% c(0,1,12)),aes(x=as.factor(sampling),y=spr,color=treatment),alpha=.3,position=position_dodge(0.3))+
    geom_point(data=ispr.sum,aes(x=as.factor(sampling),y=spr.m,color=treatment),alpha=.7,size=5,position=position_dodge(0.3))+
    geom_errorbar(data=ispr.sum,aes(x=as.factor(sampling),ymin=spr.m-spr.sd,ymax=spr.m+spr.sd,color=treatment),alpha=.7,width=.3,position=position_dodge(0.3))+ 
    scale_color_viridis_d(option="A", begin=0, end=0.6,name="Treatment",labels=c("Control","Structure","Sponge"))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          plot.margin=margin(t=1,r=1,b=5,l=5))+
    add_phylopic(crabpng,x=.5,y=5,ysize = 3,alpha = 1)+
#    ggtitle("Invertebrates")+
    ylab("")+
    xlab("Months into the Experiment"))

layout1<-c(
  area(t=1,b=2,l=1,r=1),
#  area(t=1,b=1,l=10,r=10),
  area(t=1,b=1,l=2,r=9),
  area(t=2,b=2,l=2,r=9),
  area(t=3,b=6,l=1,r=5),
  area(t=3,b=6,l=6,r=10))

plot(layout1)

wrap_elements(grid::textGrob("Species Richness",rot=90,gp=gpar(fontsize=14)))+fspr+ispr+
  p1+p12+plot_layout(design=layout1,guides="collect")

ggsave("figures/Species_Richness_Turnover.jpg",dpi=300,width=10,height=6)
