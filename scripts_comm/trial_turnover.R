install.packages("codyn")
library(codyn)
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


