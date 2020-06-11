# ---- algaeplot ----
a2<-algae %>%
  pivot_wider(names_from = taxa,values_from = abundance,values_fill = list(abundance=0))%>%
  mutate(total=rowSums(select(.,-treatment,-plot,-sampling)))%>% # whenever you are using a non-tidyverse function within a 
  #tidyverse pipeline you have to use the "." to tell it to use the dataset you had been working with.
  pivot_longer(4:11,names_to = "taxa",values_to = "abundance")%>%
  mutate(prop=abundance/total,
         season=case_when(
           sampling==1~"summer",
           sampling==2~"summer",
           sampling==4~"summer",
           sampling==3~"winter",
           sampling==5~"winter"),
         plot=paste(treatment,plot))
# now I'm going to make my plot
(c1<-ggplot(data=a2)+
    geom_bar(aes(x=sampling,y=prop,fill=taxa),position="stack",stat="identity")+
    facet_wrap(~plot)+
    geom_rect(aes(xmin=2.5,xmax=3.5,ymin=-.01,ymax=1.01),alpha=0.009,size=1,color="black")+
    geom_rect(aes(xmin=4.5,xmax=5.5,ymin=-.01,ymax=1.01),alpha=0.009,size=1,color="black"))

# ---- algaediv.p1----
a3<-algae %>%
    pivot_wider(names_from = taxa,values_from = abundance,values_fill = list(abundance=0))%>%
    mutate(div=diversity(select(.,-treatment,-plot,-sampling)),
           spr=specnumber(select(.,-treatment,-plot,-sampling)),
           season=case_when(
             sampling==1~"summer",
             sampling==2~"summer",
             sampling==4~"summer",
             sampling==3~"winter",
             sampling==5~"winter"))%>%
    select(treatment,plot,sampling,season,div,spr)

bp<-ggplot(data=a3) 
colrs<-c(rep(c("black","red","blue","orange","gray50"),3))
(div.sm<-bp+geom_line(aes(x=sampling,y=div,color=as.factor(plot)))+
   scale_color_manual(values=colrs)+
   facet_wrap(season~treatment,scales="free"))

# ---- algaediv.p2----
(div.smt<-bp+geom_boxplot(aes(x=as.factor(sampling),y=div,fill=season))+
     facet_wrap(~treatment))