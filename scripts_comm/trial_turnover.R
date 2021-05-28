# install.packages("codyn")
# remotes::install_github("sckott/rphylopic")
library(codyn)
library(rphylopic)
library(ggplot2)
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


# plots

ggplot(fish.turnover)+
  #  geom_line(aes(x=sampling, y=appearance,color=treatment,group=plot))
  #  geom_line(aes(x=sampling, y=disappearance,color=treatment,group=plot))
  #  geom_line(aes(x=sampling, y=total,color=treatment,group=plot))
  geom_point(aes(y=appearance,x=disappearance,color=treatment,size=spr))+
  facet_wrap(~sampling)


# find animal images
# http://phylopic.org/image/browse/
# try Hyperprosopon argenteum for fish (http://phylopic.org/image/0b9cdf1f-ccbc-4922-8cf6-60f90d07107e/) and blue crab for inverts
fishpng <- image_data("0b9cdf1f-ccbc-4922-8cf6-60f90d07107e", size = 128)[[1]]
crabpng <- image_data("9958579e-5e63-4b7c-8e76-9b1a92d7f7ca", size = 256)[[1]]

# or use an image saved to our project folder 
# library(png)
# fishpng <- readPNG("scripts_comm/Hyperprosopon_argenteum.png")

fish.turnover$group <- "fish"
inv.turnover$group <- "invert"
turnoverdat <- bind_rows(inv.turnover, fish.turnover) #%>% rename(lost = disappearance, gained = appearance)

unique(turnoverdat$group)
png_list <- list(crabpng,fishpng)

fish.turnover$spr
inv.turnover$spr

plot_png <- function(dat, png, 
                     png_group = "group",
                     x = "disappearance", 
                     y = "appearance",
                     scal_fac = c(100)
  ) {
  dat$x <- dat[[x]]
  dat$y <- dat[[y]]
  dat$group <- dat[[png_group]]

  p <- ggplot(dat, aes(x, y)) +
    geom_blank() +
    coord_fixed(xlim = c(-0.03, 0.8)
      , ylim = c(0, 0.8)
      ) +
    xlab(x) + ylab(y) +
    ggsidekick::theme_sleek(base_size = 18)

for (i in 1:length(unique(dat$group))) {
  # browser()
  d <- filter(dat, treatment == "real" & group == unique(dat$group)[i]) 
  for (j in 1:nrow(d)) {
    p <- p + add_phylopic(png[[i]], 0.9,
      d$x[j], d$y[j],
      col = "#E95562FF",
      ysize = (d$spr[j] + 2) / scal_fac[i]
    )
  }
  d <- filter(dat, treatment == "fake" & group == unique(dat$group)[i])
  for (j in 1:nrow(d)) {
    p <- p + add_phylopic(png[[i]], 0.9,
      d$x[j]+runif(1, -0.05, 0.05), d$y[j]+runif(1, -0.05, 0.05),
      col = "#5A167EFF",
      ysize = (d$spr[j] + 2) / scal_fac[i]
    )
  }
  d <- filter(dat, treatment == "blank"& group == unique(dat$group)[i])
  for (j in 1:nrow(d)) {
    p <- p + add_phylopic(png[[i]], 0.9,
      d$x[j]+runif(1, -0.05, 0.05), d$y[j]+runif(1, -0.01, 0.05),
      col = "black",
      ysize = (d$spr[j] + 2) / scal_fac[i]
    )
  }
}
  
  
  p
}

(p1 <- plot_png(filter(turnoverdat, sampling == 1), png = png_list, scal_fac = c(70, 100))+ ylab("Proportion of Species Gained") + xlab("Proportion of Species Lost") )

(p12 <- plot_png(filter(turnoverdat, sampling == 12), png = png_list, scal_fac = c(70, 100)) + ylab("") + xlab("Proportion of Species Lost")+ theme(axis.text.y = element_blank()))

library(patchwork)

p1 + p12 + plot_layout()
