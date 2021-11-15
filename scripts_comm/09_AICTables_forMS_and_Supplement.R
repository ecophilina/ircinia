# script to make a table of best fit models

source("scripts_comm/08_import_results.R")
if(!require(kableExtra))install.packages("kableExtra");library(kableExtra)

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

# tables of full AIC results for supplement

kbl(fish.models[,c(-1:-2,-7,-10)],booktabs=T,escape = FALSE,
    col.names = c("","$k$","$AIC_c$","$\\Delta AIC_c$",
                  "$Akaike Weight$","$LL$"),
    caption="Table SX. All models tested evaluating competing hypotheses of the drivers of fish abundance and species richness.")%>%
  kable_styling(latex_options = "scale_down",font_size=11)%>%
  pack_rows("Abundance",1,14)%>%
  pack_rows("Species Richness",15,28)


kbl(invert.models[,c(-1:-2,-7,-10)],booktabs=T,escape = FALSE,
    col.names = c("","$k$","$AIC_c$","$\\Delta AIC_c$",
                  "$Akaike Weight$","$LL$"),
    caption="Table SX. All models tested evaluating competing hypotheses of the drivers of non-clonal invertebrate abundance and species richness.")%>%
  kable_styling(latex_options = "scale_down",font_size=11)%>%
  pack_rows("Abundance",1,14)%>%
  pack_rows("Species Richness",15,28)

kbl(cispr[,c(-1:-2,-7,-10)],booktabs=T,escape = FALSE,
    col.names = c("","$k$","$AIC_c$","$\\Delta AIC_c$",
                  "$Akaike Weight$","$LL$"),
    caption="Table SX. All models tested evaluating competing hypotheses of the drivers of clonal invertebrate species richness.")%>%
  kable_styling(latex_options = "scale_down",font_size=11)

# table of best models for the results

best.models<-bind_rows(fish.models,invert.models,cispr)%>%
  filter(Delta_AICc<=2)%>%
  arrange(Response)

kbl(best.models[,c(-1:-2,-7,-10)],booktabs=T,escape = FALSE,
    col.names = c("","$k$","$AIC_c$","$\\Delta AIC_c$",
                  "$\\it Akaike Weight$","$\\it LL$"),
    caption="Table 1. Best models, i.e., those within two 
    $\\Delta AIC_c$ of the top model, explaining changes in fish
    and non-clonal invertebrate abundance and species richness as well 
    as clonal invertebrate species richness.")%>%
  kable_styling(latex_options = "scale_down",font_size=11)%>%
  pack_rows("Abundance",1,4,bold=T)%>%
  pack_rows("Fish",1,1,bold=F,italic=TRUE)%>%
  pack_rows("Non-clonal Invertebrates",2,4,bold=F,italic=TRUE)%>%
  pack_rows("Species Richness",5,10,bold=T)%>%
  pack_rows("Fish",5,5,bold=F,italic=TRUE)%>%
  pack_rows("Non-clonal Invertebrate",6,7,bold=F,italic=TRUE)%>%
  pack_rows("Clonal Invertebrates",8,10,bold=F,italic=TRUE)

  
