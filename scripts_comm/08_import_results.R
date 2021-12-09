# script to bring in results for paper
library(tidyverse)

#Colonial inverts

#species richness
col.aic<-read_rds("working_data/ColInvAIC_Results.rds")
cispr.struct.sum<-read_rds("working_data/cispr_struct_sum.rds")
cispr.prod.sum<-read_rds("working_data/cispr_prod_sum.rds")
cispr.alg.sum<-read_rds("working_data/cispr_alg_sum.rds")

#Abundance
#No abundance for colonial inverts

#turnover

#compositional vectors
col.inv.vl.aov.sum<-read_rds("working_data/ColInvVLSum.rds")
col.inv.angle.aov.sum<-read_rds("working_data/ColInvAngleSum.rds")

#Non-clonal Inverts

# species richness
inv.aic<-read_rds("working_data/InvSprAIC_Results.rds")
ispr.treat.sum<-read_rds("working_data/InvSprTreatSum.rds")
ispr.alg.sum<-read_rds("working_data/InvSprAlgSum.rds")

#Abundance
inv.abund.aic<-read_rds("working_data/InvAbundAIC_Results.rds")
ia.treat.sum<-read_rds("working_data/InvAbundTreatSum.rds")
ia.treat.prod.sum<-read_rds("working_data/InvAbundTreatProdSum.rds")
ia.treat.struct.sum<-read_rds("working_data/InvAbundTreatStructSum.rds")

#turnover
invgain1<-read_rds("working_data/invgain1.rds")
invgain12<-read_rds("working_data/invgain12.rds")
invloss1<-read_rds("working_data/invloss1.rds")
invloss12<-read_rds("working_data/invloss12.rds")
invgain12.tuk<-read_rds("working_data/invgain12_tukey.rds")


#compositional vectors
inv.vl.aov.sum<-read_rds("working_data/InvVLSum.rds")
inv.angle.aov.sum<-read_rds("working_data/InvAngleSum.rds")

#Fish

#Species richness
fspr.aic<-read_rds("working_data/FishSprAIC.rds")
fspr.treat.sum<-read_rds("working_data/FishSprTreatSum.rds")

#Abundance
fabund.aic<-read_rds("working_data/FishAbundAIC.rds")
fa.treat.sum<-read_rds("working_data/FishAbundTreatSum.rds")

# turnover
fishgain1<-read_rds("working_data/fishgain1.rds")
fishgain12<-read_rds("working_data/fishgain12.rds")
fishloss1<-read_rds("working_data/fishloss1.rds")
fishloss12<-read_rds("working_data/fishloss12.rds")
fishgain12.tuk<-read_rds("working_data/fishgain12_tukey.rds")
fishloss12.tuk<-read_rds("working_data/fishloss12_tukey.rds")

#compositional vectors
fish.vl.aov.sum<-read_rds("working_data/FishVLSum.rds")
fish.vl.aov.tuk<-read_rds("working_data/FishVLTuk.rds")
fish.angle.aov.sum<-read_rds("working_data/FishAngleSum.rds")
fish.angle.aov.tuk<-read_rds("working_data/FishAngleTuk.rds")

ls(pattern=".aic")
# fish models
fabund<-fabund.aic%>%
  separate(Modnames,into=c("cr","model"),sep=3)%>%
  mutate(Community="Fish",
         Response="Abundance",
         Model=case_when(
           model=="treat"~"Sponge Presence",
           model=="treat.prod"~"Sponge Presence + Seagrass Productivity",
           model=="treat.struct"~"Sponge Presence + Seagrass Structure",
           model=="treat.alg"~"Sponge Presence + Macroalgal Structure",
           model=="alg"~"Macroalgal Structure",
           model=="treat.prod.struct"~"Sponge Presence + Seagrass Productivity & Structure",
           model=="treat.prod.alg"~"Sponge Presence + Seagrass Productivity + Macroalgal Structure",
           model=="struct.alg"~"Seagrass & Macroalgal Structure",
           model=="treat.struct.alg"~"Sponge Presence + Seagrass & Macroalgal Structure",
           model=="prod.alg"~"Seagrass Productivity + Macroalgal Structure",
           model=="struct"~"Seagrass Structure",
           model=="prod"~"Seagrass Productivity",
           model=="prod.struct"~"Seagrass Productivity & Structure",
           model=="full"~"Sponge Presence + Seagrass Productivity & Structure + Macroalgal Structure"),
         across(AICc:Cum.Wt,round,2))%>%
  select(Community,Response,Model,K:Cum.Wt)

fspr<-fspr.aic%>%
  separate(Modnames,into=c("cr","model"),sep=5)%>%
  mutate(Community="Fish",
         Response="Species Richness",
         Model=case_when(
           model=="treat"~"Sponge Presence",
           model=="treat.prod"~"Sponge Presence + Seagrass Productivity",
           model=="treat.struct"~"Sponge Presence + Seagrass Structure",
           model=="treat.alg"~"Sponge Presence + Macroalgal Structure",
           model=="alg"~"Macroalgal Structure",
           model=="treat.prod.struct"~"Sponge Presence + Seagrass Productivity & Structure",
           model=="treat.prod.alg"~"Sponge Presence + Seagrass Productivity + Macroalgal Structure",
           model=="struct.alg"~"Seagrass & Macroalgal Structure",
           model=="treat.struct.alg"~"Sponge Presence + Seagrass & Macroalgal Structure",
           model=="prod.alg"~"Seagrass Productivity + Macroalgal Structure",
           model=="struct"~"Seagrass Structure",
           model=="prod"~"Seagrass Productivity",
           model=="prod.struct"~"Seagrass Productivity & Structure",
           model=="full"~"Sponge Presence + Seagrass Productivity & Structure + Macroalgal Structure"),
         across(AICc:Cum.Wt,round,2))%>%
  select(Community,Response,Model,K:Cum.Wt)

fish.models<-bind_rows(fabund,fspr)

# non-clonal invert models
ia<-inv.abund.aic%>%
  separate(Modnames,into=c("cr","model"),sep=3)%>%
  mutate(Community="Non-clonal Invertebrates",
         Response="Abundance",
         Model=case_when(
           model=="treat"~"Sponge Presence",
           model=="treat.prod"~"Sponge Presence + Seagrass Productivity",
           model=="treat.struct"~"Sponge Presence + Seagrass Structure",
           model=="treat.alg"~"Sponge Presence + Macroalgal Structure",
           model=="alg"~"Macroalgal Structure",
           model=="treat.prod.struct"~"Sponge Presence + Seagrass Productivity & Structure",
           model=="treat.prod.alg"~"Sponge Presence + Seagrass Productivity + Macroalgal Structure",
           model=="struct.alg"~"Seagrass & Macroalgal Structure",
           model=="treat.struct.alg"~"Sponge Presence + Seagrass & Macroalgal Structure",
           model=="prod.alg"~"Seagrass Productivity + Macroalgal Structure",
           model=="struct"~"Seagrass Structure",
           model=="prod"~"Seagrass Productivity",
           model=="prod.struct"~"Seagrass Productivity & Structure",
           model=="full"~"Sponge Presence + Seagrass Productivity & Structure + Macroalgal Structure"),
         across(AICc:Cum.Wt,round,2))%>%
  select(Community,Response,Model,K:Cum.Wt)

ispr<-inv.aic%>%
  separate(Modnames,into=c("cr","model"),sep=5)%>%
  mutate(Community="Non-clonal Invertebrates",
         Response="Species Richness",
         Model=case_when(
           model=="treat"~"Sponge Presence",
           model=="treat.prod"~"Sponge Presence + Seagrass Productivity",
           model=="treat.struct"~"Sponge Presence + Seagrass Structure",
           model=="treat.alg"~"Sponge Presence + Macroalgal Structure",
           model=="alg"~"Macroalgal Structure",
           model=="treat.prod.struct"~"Sponge Presence + Seagrass Productivity & Structure",
           model=="treat.prod.alg"~"Sponge Presence + Seagrass Productivity + Macroalgal Structure",
           model=="struct.alg"~"Seagrass & Macroalgal Structure",
           model=="treat.struct.alg"~"Sponge Presence + Seagrass & Macroalgal Structure",
           model=="prod.alg"~"Seagrass Productivity + Macroalgal Structure",
           model=="struct"~"Seagrass Structure",
           model=="prod"~"Seagrass Productivity",
           model=="prod.struct"~"Seagrass Productivity & Structure",
           model=="full"~"Sponge Presence + Seagrass Productivity & Structure + Macroalgal Structure"),
         across(AICc:Cum.Wt,round,2))%>%
  select(Community,Response,Model,K:Cum.Wt)


invert.models<-bind_rows(ia,ispr)
# clonal invertebrates
cispr<-col.aic%>%
  separate(Modnames,into=c("cr","model"),sep=6)%>%
  mutate(Community="Clonal Invertebrates",
         Response="Species Richness",
         Model=case_when(
           model=="treat"~"Sponge Presence",
           model=="treat.prod"~"Sponge Presence + Seagrass Productivity",
           model=="treat.struct"~"Sponge Presence + Seagrass Structure",
           model=="treat.alg"~"Sponge Presence + Macroalgal Structure",
           model=="alg"~"Macroalgal Structure",
           model=="treat.prod.struct"~"Sponge Presence + Seagrass Productivity & Structure",
           model=="treat.prod.alg"~"Sponge Presence + Seagrass Productivity + Macroalgal Structure",
           model=="struct.alg"~"Seagrass & Macroalgal Structure",
           model=="treat.struct.alg"~"Sponge Presence + Seagrass & Macroalgal Structure",
           model=="prod.alg"~"Seagrass Productivity + Macroalgal Structure",
           model=="struct"~"Seagrass Structure",
           model=="prod"~"Seagrass Productivity",
           model=="prod.struct"~"Seagrass Productivity & Structure",
           model=="full"~"Sponge Presence + Seagrass Productivity & Structure + Macroalgal Structure"),
         across(AICc:Cum.Wt,round,2))%>%
  select(Community,Response,Model,K:Cum.Wt)
# table of best models for the results

best.models<-bind_rows(fish.models,invert.models,cispr)%>%
  filter(Delta_AICc<=2)%>%
  arrange(Response)


