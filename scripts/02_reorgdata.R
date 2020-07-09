## this is a project to learn best practices of data management in R ##
## @knitr loaddata
# load original data ####
if(!require(tidyverse))install.packages('tidyverse');library(tidyverse)
if(!require(readxl))install.packages('readxl');library(readxl)

sg_grow<-read_xlsx(here("Original_data","ForFinella_Transplant_data.xlsx"),sheet="Sg_growth")
algae<-read_xlsx(here("Original_data","ForFinella_Transplant_data.xlsx"),sheet="Algae")
sg_shoot<-read_xlsx(here("Original_data","ForFinella_Transplant_data.xlsx"),sheet="sg_shootdensity")
fish<-read_xlsx(here("Original_data","ForFinella_Transplant_data.xlsx"),sheet="fish")
sg_nuts<-read_xlsx(here("Original_data","ForFinella_Transplant_data.xlsx"),sheet="sg_nutrients")
inverts<-read_xlsx(here("Original_data","ForFinella_Transplant_data.xlsx"),sheet="inverts")

## this script is to organize data in a better format ##
## @knitr reorgdata
## reorganize data ####

algae2<-algae %>% # this tells R that I'd like to work with this data set
  pivot_longer(c(-Treatment,-plot,-sampling),names_to = "taxa",values_to = "abundance")%>% 
  # this tells R I'd like to take all the columns other than those indicated with the "-" sign, and turn them into rows
  # the names of the columns will be put into the column "taxa" and the values into "abundance"
  subset(abundance!=0)%>%
  #this says get rid of rows where abundance is equal to 0 (alternative only keep rows where abundance doesn't equal 0)
  rename(treatment=Treatment)
  # this renames the Treatment column so it matches the other data sets

fish2<-fish%>%
  select(c(-treat,-time,-group))%>%
  #this gets rid of the columns with those names. If you want to select columns to keep omit the "-"
  pivot_longer(c(-treatment,-plot,-sampling),names_to = "taxa",values_to = "abundance")%>%
  subset(abundance!=0)

inverts2<-inverts %>%
  select(c(-time,-ID,-date))%>%
  pivot_longer(c(-Treatment,-plot,-sampling),names_to = "taxa",values_to = "abundance")%>%
  subset(abundance!=0)%>%
  rename(treatment=Treatment)

sg_grow2<-sg_grow%>%
  select(c(-treat,-ID2,-ID))%>%
  rename(sampling=time)


sg_nuts2<-sg_nuts %>%
  select(c(-ID,-CNRatio,-NP,-CP))%>%
  rename(mg='Weight(mg)',sampling=samp,treatment=treat)%>%
  mutate(treatment=case_when(
    treatment=="B"~"blank",
    treatment=="F"~"fake",
    treatment=="R"~"real"),
    plot=case_when(
      treatment=="blank"~plot,
      treatment=="fake"~plot+5,
      treatment=="real"~plot+10))
#Mutate allows you to create new variables or redefine old ones. case_when allows you to say what
# the new variable will be given the value of another variable or, as in this case, what the variable was before
# this code basically says when treatment is equal to "B" I'd like to change it to "blank".
# again I did this to keep the data sets consistent

sg_shoot2<-sg_shoot%>%
  select(-treat)%>%
  rename(sampling=time)%>%
  mutate(treatment=case_when(
    treatment=="Blank"~"blank",
    treatment=="Fake"~"fake",
    treatment=="Real"~"real"),
    T.SD=ifelse(T.SD==111,84,T.SD))

# write reorganized data ####
# @knitr writereorg
# now we're going to save all the reorganized data sets the code below shows you how to include the
# date in the name of a file you are saving. 
write.csv(algae2,here("working_data",paste0("algae",".csv")))
write.csv(fish2,here("working_data",paste0("fish",".csv")))
write.csv(inverts2,here("working_data",paste0("inverts",".csv")))
write.csv(sg_grow2,here("working_data",paste0("sg_grow",".csv")))
write.csv(sg_nuts2,here("working_data",paste0("sg_nuts",".csv")))
write.csv(sg_shoot2,here("working_data",paste0("sg_shoot",".csv")))
