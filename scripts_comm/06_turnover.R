# install.packages("codyn")
#remotes::install_github("sckott/rphylopic")

library(codyn)
library(rphylopic)
library(ggplot2)
library(patchwork)
library(grid)
library(ggpubr)
library(tidyverse)

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


#### Make multipanel figure for paper ####
# find animal images
# http://phylopic.org/image/browse/
# try Hyperprosopon argenteum for fish (http://phylopic.org/image/0b9cdf1f-ccbc-4922-8cf6-60f90d07107e/) and blue crab for inverts
fishpng <- image_data("0b9cdf1f-ccbc-4922-8cf6-60f90d07107e", size = 256)[[1]]
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

#### plot species richness vs time ####

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
    geom_point(data=fish.uni%>%filter(sampling %in% c(0,1,12)),
      aes(x=as.factor(sampling),y=spr,color=treatment),alpha=.3,position=position_dodge(0.3))+
    geom_point(data=fspr.sum,
      aes(x=as.factor(sampling),y=spr.m,color=treatment),alpha=.7,size=5,position=position_dodge(0.3))+
    geom_errorbar(data=fspr.sum,
      aes(x=as.factor(sampling),ymin=spr.m-spr.sd,ymax=spr.m+spr.sd,color=treatment),alpha=.7,width=.3,position=position_dodge(0.3))+ 
    scale_color_viridis_d(option="A", begin=0, end=0.65,name="Treatment",labels=c("Control","Structure","Sponge"))+
    # theme_bw()+
    ggsidekick::theme_sleek(base_size = 16) +
    theme(
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      legend.position = "none",
      # panel.grid = element_blank(),
      # plot.margin=margin(t=5,r=1,b=1,l=5),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank())+
    add_phylopic(fishpng,x=.75,y=5,ysize = 1.5,alpha=1)+
    ggtitle("Months into the Experiment"))

(ispr<-ggplot()+
    geom_point(data=inv.uni%>%filter(sampling %in% c(0,1,12)),
      aes(x=as.factor(sampling),y=spr,color=treatment),alpha=.3,position=position_dodge(0.3))+
    geom_point(data=ispr.sum,
      aes(x=as.factor(sampling),y=spr.m,color=treatment),alpha=.7,size=5,position=position_dodge(0.3))+
    geom_errorbar(data=ispr.sum,
      aes(x=as.factor(sampling),ymin=spr.m-spr.sd,ymax=spr.m+spr.sd,color=treatment),alpha=.7,width=.3,position=position_dodge(0.3))+ 
    scale_color_viridis_d(option="A", begin=0, end=0.65,name="Treatment",labels=c("Control","Structure","Sponge"))+
    # theme_bw()+
    ggsidekick::theme_sleek(base_size = 16) +
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
      legend.position = "none",
      # panel.grid = element_blank(),
      # plot.margin=margin(t=1,r=1,b=5,l=5)
    )+
    add_phylopic(crabpng,x=0.75,y=5,ysize = 3,alpha = 1)+
    #    ggtitle("Invertebrates")+
    ylab("")+
    xlab("Months into the Experiment"))

cspr.sum<-col.inv.uni%>%
  group_by(treatment,sampling)%>%
  summarize(spr.m=mean(spr),
    spr.sd=sd(spr))%>%
  filter(sampling %in% c(0,1,12))

(cspr<-ggplot()+
    geom_point(data=col.inv.uni%>%filter(sampling %in% c(0,1,12)),
      aes(x=as.factor(sampling),y=spr,color=treatment),alpha=.3,position=position_dodge(0.3))+
    geom_point(data=cspr.sum,
      aes(x=as.factor(sampling),y=spr.m,color=treatment),alpha=.7,size=5,position=position_dodge(0.3))+
    geom_errorbar(data=cspr.sum,
      aes(x=as.factor(sampling),ymin=spr.m-spr.sd,ymax=spr.m+spr.sd,color=treatment),alpha=.7,width=.3,position=position_dodge(0.3))+ 
    scale_color_viridis_d(option="A", begin=0, end=0.65,name="Treatment",labels=c("Control","Structure","Sponge"))+
    ggsidekick::theme_sleek(base_size = 16) +
    theme(
      axis.title.x = element_blank(),
      # axis.title.y = element_blank(),
      legend.position = "none"
    )+
    add_phylopic(crabpng,x=0.75,y=5,ysize = 3,alpha = 1)+
    ylab("Species Richness")+
    ggtitle("Months into the Experiment"))


# make separate legend for placement using patchwork
l1 <- ggpubr::get_legend(ispr + theme(legend.position = c(0.9,0.9)))


#### make turnover plots with animal shapes ####

fish.turnover$group <- "fish"
inv.turnover$group <- "invert"
col.inv.turnover$group <- "clonal"  

fish.turnover$spr
inv.turnover$spr
col.inv.turnover$spr

# specialized function for plotting different animal shapes ("groups") with or without colouring by treatments ( with categories of "real","fake","control")
plot_png <- function(dat, 
                     png_list, # list() with # of images = groups, in order of unique(dat$group)
                     png_group = "group", # variable that determines which image
                     plot_treatments = T, 
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
    ggsidekick::theme_sleek(base_size = 16)
  
# loop over groups to get different png images
for (i in 1:length(unique(dat$group))) {
  d <- filter(dat, group == unique(dat$group)[i])  # generalizes function to work without treatment column
  if(plot_treatments){
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
  }
  for (j in 1:nrow(d)) {
    if(plot_treatments){
    p <- p + add_phylopic(png_list[[i]], alpha,
      d$x[j]+runif(1, -0.02, 0.02), d$y[j]+runif(1, -0.02, 0.02),
      col = "black",
      ysize = (d$spr[j] + 2) / scal_fac[i]
    )
    } else{
      p <- p + add_phylopic(png_list[[i]], alpha,
        d$x[j], d$y[j],
        col = "black",
        ysize = (d$spr[j] + 2) / scal_fac[i]
      )
    }
  }
}
  p 
}


# I tried adding tunicates but it just confused things
turnoverdat <- bind_rows(inv.turnover, fish.turnover) 
unique(turnoverdat$group)
png_list <- list(crabpng, fishpng)


# each jitter is random, so just rerun plot code if points are falling on edge of plot area
(p1 <- plot_png(filter(turnoverdat, sampling == 1), png_list, scal_fac = c(80, 100)) +
    ylab("Proportion of Species Gained") + xlab("Proportion of Species Lost") +
    # theme(axis.title.x = element_blank(),axis.title.y = element_blank()) + # turn off if not add global axes
    theme(plot.title=element_text(hjust=0.5))+
    ggtitle("1 month"))

(p12 <- plot_png(filter(turnoverdat, sampling == 12), png_list, scal_fac = c(80, 100)) + 
    ylab("Proportion of Species Gained") +  xlab("Proportion of Species Lost") +
    # theme(axis.title.x = element_blank()) + # turn off if not add global axes
    theme(plot.title=element_text(hjust=0.5), axis.title.y = element_blank(), 
      axis.text.y = element_blank(), axis.ticks.y = element_blank())+
    ggtitle("12 months"))

# create legend of sizes of shapes based on species richness 
ldatmin <- turnoverdat %>% 
  select(group, disappearance, appearance, spr) %>%   
  mutate(disappearance = ifelse(group=="fish", 0.15, 0.44), appearance = 0.56, 
    spr = 0) %>%
  distinct()
ldatmid <- turnoverdat %>% 
  select(group, disappearance, appearance, spr) %>%   
  mutate(disappearance = ifelse(group=="fish", 0.15, 0.44), appearance = 0.65, 
    spr = round(max(spr)/2)-1) %>%
  distinct()
ldatmax <- turnoverdat %>% 
  select(group, disappearance, appearance, spr) %>%   
  mutate(disappearance = ifelse(group=="fish", 0.15, 0.44), appearance = 0.75, 
    spr = round(quantile(spr, 0.95))) %>%
  distinct()

ldat <- bind_rows(ldatmin, ldatmid, ldatmax)

(l2 <- plot_png(ldat, png_list, plot_treatments = F, scal_fac = c(70, 100), xlim = c(0.05, 0.5)) + 
    geom_text( aes( x=0.3, y=ldatmax$appearance[1], label= ldatmax$spr[1])) +
    geom_text( aes( x=0.3, y=ldatmid$appearance[1], label= ldatmid$spr[1])) +
    geom_text( aes( x=0.3, y=ldatmin$appearance[1], label= ldatmin$spr[1])) +
    theme_void()+ theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("Species Richness"))

# ## if wanting just turnover by themselves... 
layout1 <- c(
  area(t=1, b=14, l=1, r=6),
  area(t=1, b=14, l=7, r=12),
  area(t=2, b=15, l=13, r=15)
)
p1 + p12 + l2 +
  plot_layout(design=layout1)

ggsave("figures/turnover-plots-legend.png", width = 12.5, height = 6)


# since tunicates show a different pattern, I'm putting them on different panels
(p1c <- plot_png(filter(col.inv.turnover, sampling == 1), png_list = list(clonalpng), scal_fac = c(80)) +
    ylab("Proportion of Species Gained") + xlab("Proportion of Species Lost")+
    # theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
    ggtitle("1 month into experiment"))

(p12c <- plot_png(filter(col.inv.turnover, sampling == 12), png_list = list(clonalpng), scal_fac = c(80)) + 
    ylab("") + xlab("Proportion of Species Lost") +
    # theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
    ggtitle("12 months into experiment") + 
    theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()))

# create legend for tunicates
ldat3 <- filter(ldat, group== "fish") # manually convert legend to scale for tunicates
ldat3$group <- "clonal"
ldat3$spr <- c(0, round(max(col.inv.turnover$spr)/2)-1, round(quantile(col.inv.turnover$spr, 0.95)))
ldat3$disappearance <- c(0.1, 0.1, 0.1)
ldat3$appearance <- c(0.66, 0.72, 0.79)

(l3 <- plot_png(ldat3, png_list = list(clonalpng), plot_treatments = F, 
  scal_fac = c(80), xlim = c(0.05, 0.25)) +
    geom_text( aes( x=0.2, y=0.79, label= round(quantile(col.inv.turnover$spr, 0.95)))) +
    geom_text( aes( x=0.2, y=0.72, label= round(max(col.inv.turnover$spr)/2)-1)) +
    geom_text( aes( x=0.2, y=0.66, label= "0")) +
    theme_void()+ theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("Species Richness"))

# layout2 <- c(
#   area(t=1, b=14, l=1, r=6),
#   area(t=1, b=14, l=7, r=12),
#   area(t=2, b=14, l=4, r=6)
# )
# plot(layout2)
# p1c + p12c + l3 + plot_layout(design=layout2)
# 
# ggsave("figures/turnover-plots-tunicates.png", width = 11, height = 6)

layout2b <- c(
  area(t=1, b=7, l=1, r=9),
  area(t=2, b=3, l=8, r=11),
  area(t=8, b=21, l=1, r=6),
  area(t=8, b=21, l=7, r=12),
  area(t=4, b=18, l=10, r=12)
)

cspr + l1 + p1c + p12c + l3 + plot_layout(design=layout2b)
ggsave("figures/turnover-plots-tunicates2.jpg", width = 10, height = 9)


#### make combined plot for fish and inverts only ####

# remove axes so I can add global ones
(p1 <- p1 + theme(axis.title.x = element_blank(),axis.title.y = element_blank()))
(p12 <- p12 + theme(axis.title.x = element_blank()))

layout<-c(
  area(t=2,b=5,l=1,r=1),# A - y axis label
  area(t=7,b=11,l=1,r=1),# B - y axis label
  area(t=12,b=12,l=2,r=7), # B - y axis label
  area(t=1,b=3,l=2,r=5),# fish richness plot
  area(t=4,b=6,l=2,r=5),# invert richness plot
  area(t=7,b=11,l=2,r=4), # 1 month turnover
  area(t=7,b=11,l=5,r=7), # 12 months turnover
  area(t=2,b=3,l=6,r=6), # treatment legend
  area(t=5,b=7,l=6,r=7), # animal richness legend
  # area(t=0,b=1,l=1,r=1),# A. tag (but positioning didn't work so use title option in plot_annotation to add this)
  area(t=6,b=7,l=1,r=1)# B. tag
  )

plot(layout)

wrap_elements(grid::textGrob("Species Richness",rot=90,vjust =2,gp=gpar(fontsize=16))) +
  wrap_elements(grid::textGrob("Proportion Gained",rot=90,vjust =2,gp=gpar(fontsize=16))) +
  wrap_elements(grid::textGrob("Proportion Lost",vjust =0,gp=gpar(fontsize=16))) +
  fspr + ispr + # richness plots
  p1 + p12 + # turnover plots
  l1 + l2 + # legends
  # wrap_elements(grid::textGrob("A.",vjust =-1,gp=gpar(fontsize=16))) +
  wrap_elements(grid::textGrob("B.",vjust =0.5,hjust=1.1,gp=gpar(fontsize=16))) +
  plot_layout(design=layout) + 
  plot_annotation( # using to add tag that wouldn't cooperate otherwise
    title = '    A.', theme = theme(plot.title = element_text(size = 16, vjust = -3))
  )

ggsave("figures/Species_Richness_Turnover_A_B.jpg",dpi=300,width=8,height=9.5)

#Figure for abundance of fish and non clonal inverts----
#Plot abundance vs time 

fab.sum<-fish.uni%>%
  group_by(treatment,sampling)%>%
  summarize(abund.m=mean(f.abund),
            abund.sd=sd(f.abund))%>%
  filter(sampling %in% c(0,1,12))

iab.sum<-inv.uni%>%
  group_by(treatment,sampling)%>%
  summarize(abund.m=mean(i.abund),
            abund.sd=sd(i.abund))%>%
  filter(sampling %in% c(0,1,12))



(fab<-ggplot()+
    geom_point(data=fish.uni%>%filter(sampling %in% c(0,1,12)),
               aes(x=as.factor(sampling),y=f.abund,color=treatment),alpha=.3,position=position_dodge(0.3))+
    geom_point(data=fab.sum,
               aes(x=as.factor(sampling),y=abund.m,color=treatment),alpha=.7,size=5,position=position_dodge(0.3))+
    geom_errorbar(data=fab.sum,
                  aes(x=as.factor(sampling),ymin=abund.m-abund.sd,ymax=abund.m+abund.sd,color=treatment),alpha=.7,width=.3,position=position_dodge(0.3))+ 
    scale_color_viridis_d(option="A", begin=0, end=0.65,name="Treatment",labels=c("Control","Structure","Sponge"))+
    # theme_bw()+
    ggsidekick::theme_sleek(base_size = 16) +
    theme(
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      legend.position = "none",
      # panel.grid = element_blank(),
      # plot.margin=margin(t=5,r=1,b=1,l=5),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank())+
    add_phylopic(fishpng,x=.75,y=6,ysize = 1.5,alpha=1))


(iab<-ggplot()+
    geom_point(data=inv.uni%>%filter(sampling %in% c(0,1,12)),
               aes(x=as.factor(sampling),y=i.abund,color=treatment),alpha=.3,position=position_dodge(0.3))+
    geom_point(data=iab.sum,
               aes(x=as.factor(sampling),y=abund.m,color=treatment),alpha=.7,size=5,position=position_dodge(0.3))+
    geom_errorbar(data=iab.sum,
                  aes(x=as.factor(sampling),ymin=abund.m-abund.sd,ymax=abund.m+abund.sd,color=treatment),alpha=.7,width=.3,position=position_dodge(0.3))+ 
    scale_color_viridis_d(option="A", begin=0, end=0.65,name="Treatment",labels=c("Control","Structure","Sponge"))+
    # theme_bw()+
    ggsidekick::theme_sleek(base_size = 16) +
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
          legend.position = "none",
          # panel.grid = element_blank(),
          # plot.margin=margin(t=1,r=1,b=5,l=5)
    )+
    add_phylopic(crabpng,x=0.75,y=11,ysize = 4,alpha = 1)+
    #    ggtitle("Invertebrates")+
    ylab(""))



# make separate legend for placement using patchwork
l1 <- ggpubr::get_legend(iab + theme(legend.position = c(0.5,0.5)))




## layout
layouta <- c(
  area(t=2,b=13,l=1,r=1),# y label
  area(t=1, b=7, l=2, r=8), #fish abundance
  area(t=8, b=14, l=2, r=8), #non-clonal invert abundance
  area(t=15,b=15,l=2,r=8), # x label
  area(t=7, b=7, l=9, r=10) #legend
)


plot(layouta)

wrap_elements(grid::textGrob("Abundance",rot=90,gp=gpar(fontsize=16))) +
  fab+
  iab+
  wrap_elements(grid::textGrob("Months into the Experiment",vjust =0,gp=gpar(fontsize=16))) +
  l1+
  plot_layout(design=layouta)

ggsave("figures/fish_ncinvert_abundance figure.png", width = 12.5, height = 6)
