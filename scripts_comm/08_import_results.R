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

#compositional vectors
fish.vl.aov.sum<-read_rds("working_data/FishVLSum.rds")
fish.vl.aov.tuk<-read_rds("working_data/FishVLTuk.rds")
fish.angle.aov.sum<-read_rds("working_data/FishAngleSum.rds")
fish.angle.aov.tuk<-read_rds("working_data/FishAngleTuk.rds")

