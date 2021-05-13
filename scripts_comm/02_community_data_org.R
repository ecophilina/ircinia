# organize data for community analysis
library(tidyverse)

# bring in data that has been prepared in orginal script series
source("scripts/03_reimport.R")  

# prep primary producer data
sg<-sg_shoot%>%
  mutate(sg.sd.global=mean(SD))%>%
  group_by(treatment,plot,sampling,sg.sd.global)%>%
  summarize(sg.sd=mean(SD),t.sd=mean(T.SD))%>%
  mutate(
    #sg.sd=sg.sd-sg.sd.global, # comment this line out if not centering
    sampling=case_when(
      sampling==1~0,
      sampling==2~1,
      sampling==3~5,
      sampling==4~12,
      sampling==5~17)) 

sggrow<-sg_grow%>%
  filter(dist %in% c(0,0.5))%>%
  mutate(grow.global=mean(total.growth.mm2/days))%>%
  group_by(treatment,plot,sampling,grow.global)%>%
  summarize(grow=mean(total.growth.mm2/days))%>%
  mutate(
 #   grow=grow-grow.global,
    sampling=case_when(
      sampling==1~0,
      sampling==2~1,
      sampling==3~5,
      sampling==4~12,
      sampling==5~17))


# Algae. No need to filter for has algae vs doesn't have algae because no plots had a 0 count for algae :)
alg<-algae%>%
  filter(!taxa %in% c("brown.cyanobacteria","green.cyanobacteria","dictyota"))%>%
  dplyr::select(-taxa)%>%
  group_by(treatment,plot,sampling)%>%
  summarize(abund=sum(abundance)/3)%>%
  ungroup()%>%
  distinct()%>%
  mutate(abund.global=mean(abund))%>%
  mutate(
#    abund=abund-abund.global,
    sampling=case_when(
      sampling==1~0,
      sampling==2~1,
      sampling==3~5,
      sampling==4~12,
      sampling==5~17))

productivity<-left_join(sg,sggrow)%>%
  left_join(alg)%>%
  ungroup()%>%
  mutate(sg.prod=grow*t.sd,
         pp.struct = sg.sd + abund,
         sg.prod.global = mean(sg.prod),
         pp.struct.global = mean(pp.struct),
         sg.prod.c = sg.prod - sg.prod.global,
         pp.struct.c = pp.struct - pp.struct.global,
         sg.sd.c = sg.sd - sg.sd.global,
         a.abund.c = abund - abund.global) %>%
  rename(a.abund = abund, a.abund.global = abund.global)

# create community datasets

alg.com.full<-algae%>%
  filter(!taxa %in% c("brown.cyanobacteria","green.cyanobacteria","dictyota"))%>%
  mutate(sampling=case_when(
    sampling==1~0,
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17),
    season=case_when(
      sampling==0~"summer",
      sampling==1~"summer",
      sampling==5~"winter",
      sampling==12~"summer",
      sampling==17~"winter"),
    yr=case_when(
      sampling==0~0,
      sampling==1~1,
      sampling==5~1,
      sampling==12~2,
      sampling==17~2),
    taxa=ifelse(taxa=="cladocephalus","udotea",taxa))%>%
  left_join(productivity)%>%
  pivot_wider(names_from=taxa,values_from=abundance,values_fill=0)

col.inverts<-c("white tunicate","white and green tunicate","pinkish brown tunicate",
               "black and orange tunicate","red and white tunicate","Yellow colonial tunicate",
               "black and white tunicate","black tunicate")
inv.com.full<-inverts%>%
  filter(!taxa %in% col.inverts)%>%
  mutate(sampling=case_when(
    sampling==1~0,
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17),
    season=case_when(
      sampling==0~"summer",
      sampling==1~"summer",
      sampling==5~"winter",
      sampling==12~"summer",
      sampling==17~"winter"),
    yr=case_when(
      sampling==0~0,
      sampling==1~1,
      sampling==5~1,
      sampling==12~2,
      sampling==17~2))%>%
  left_join(productivity)%>%
  pivot_wider(names_from=taxa,values_from=abundance,values_fill=0)

col.inv.com.full<-inverts%>%
  filter(taxa %in% col.inverts)%>%
  mutate(sampling=case_when(
    sampling==1~0,
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17),
    season=case_when(
      sampling==0~"summer",
      sampling==1~"summer",
      sampling==5~"winter",
      sampling==12~"summer",
      sampling==17~"winter"),
    yr=case_when(
      sampling==0~0,
      sampling==1~1,
      sampling==5~1,
      sampling==12~2,
      sampling==17~2))%>%
  left_join(productivity)%>%
  pivot_wider(names_from=taxa,values_from=abundance,values_fill=0)

fish.com.full<-fish%>%
  mutate(sampling=case_when(
    sampling==1~0,
    sampling==2~1,
    sampling==3~5,
    sampling==4~12,
    sampling==5~17),
    season=case_when(
      sampling==0~"summer",
      sampling==1~"summer",
      sampling==5~"winter",
      sampling==12~"summer",
      sampling==17~"winter"),
    yr=case_when(
      sampling==0~0,
      sampling==1~1,
      sampling==5~1,
      sampling==12~2,
      sampling==17~2))%>%
  left_join(productivity)%>%
  pivot_wider(names_from=taxa,values_from=abundance,values_fill=0)

keep.env<-c(colnames(productivity),"season","yr")

alg.env<-alg.com.full %>%
  dplyr::select(all_of(keep.env))

alg.com<-alg.com.full %>%
  dplyr::select(!all_of(keep.env))

alg.uni<-alg.env%>%
  mutate(spr=vegan::specnumber(alg.com),
         div=vegan::diversity(alg.com,index = "shannon"),
         j=div/log(spr),
         j=ifelse(is.na(j),0,j),
         a.abund = rowSums(alg.com))

inv.env<-inv.com.full %>%
  dplyr::select(all_of(keep.env))

inv.com<-inv.com.full %>%
  dplyr::select(!all_of(keep.env))

inv.uni<-inv.env%>%
  mutate(spr=vegan::specnumber(inv.com),
         div=vegan::diversity(inv.com,index = "shannon"),
         j=div/log(spr),
         j=ifelse(is.na(j),0,j),
         i.abund=rowSums(inv.com))

col.inv.env<-col.inv.com.full %>%
  dplyr::select(all_of(keep.env))

col.inv.com<-col.inv.com.full %>%
  dplyr::select(!all_of(keep.env))

col.inv.uni<-col.inv.env%>%
  mutate(spr=vegan::specnumber(inv.com),
         div=vegan::diversity(inv.com,index = "shannon"),
         j=div/log(spr),
         j=ifelse(is.na(j),0,j),
         i.abund=rowSums(inv.com))

fish.env<-fish.com.full %>%
  dplyr::select(all_of(keep.env))

fish.com<-fish.com.full %>%
  dplyr::select(!all_of(keep.env))

fish.uni<-fish.env%>%
  mutate(spr=vegan::specnumber(fish.com),
         div=vegan::diversity(fish.com,index = "shannon"),
         j=div/log(spr),
         j=ifelse(is.na(j),0,j),
         f.abund=rowSums(fish.com))

# make a function to look at residuals
glmm.resids<-function(model){
 t1 <- simulateResiduals(model)
  print(testDispersion(t1))
  plot(t1)
}

