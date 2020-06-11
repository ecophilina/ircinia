## this is a project to learn best practices of data management in R ##

## this script will import the reorganized data and do initial data explorations ##

# @knitr reimport

## import reworked data ####
files <- list.files(here("working_data"), pattern = "2020-05-08")
# this bit of code finds all the files in the working data directory that match the date I
# exported the reorganized data you can replace the date with "*.csv" if you want to find all the csv
# files in a folder- or replace the date with any text pattern really.

fd <- data.frame(files = files, dname = files) %>%
  # this line creates a data frame that I'll use to import the files in a couple lines
  separate(dname, into = c("dname", "extra"), sep = "2020")
# this bit separates out of the name of the dataset from the date and the file extension- note
# we still have one column with the full file name

for (i in 1:nrow(fd))assign(fd[i, 2], read.csv(here("working_data", fd[i, 1]))[, -1])
# this bit of code says for every row in the fd data frame assign the file in column 1
# to the name in column 2. You should now have all the reorganized data sets imported.





# div.smt<-bp+geom_boxplot(aes(x=as.factor(sampling),y=div,fill=season))+
#   facet_wrap(~treatment)
# # the patterns I was seeing at the plot level are popping out here too. only plots with real sponges have similar 
# # algal diversity at the end of the experiment as they did at the beginning (sampling 1)
# 
# 
# # Now lets do the same thing with species richness
# # @knitr algaespr ####
# 
# bp+geom_line(aes(x=sampling,y=spr,color=as.factor(plot)))+
#   scale_color_manual(values=colrs)+
#   facet_wrap(~treatment)
# 
# # this makes my brain hurt, but nothing obviously wrong here.
# # now I'm going to look at summer and winter separately
# 
# spr.sm<-bp+geom_line(aes(x=sampling,y=spr,color=as.factor(plot)))+
#   scale_color_manual(values=colrs)+
#   facet_wrap(season~treatment,scales="free")
# # looks like there's less interesting stuff going on here than with diversity- mostly because plots are kind of all
# # over the place.
# 
# spr.smt<-bp+geom_boxplot(aes(x=as.factor(sampling),y=spr,fill=season))+
#   facet_wrap(~treatment)+
#   theme_bw()+
#   theme(panel.grid = element_blank())+
#   ylab("Richness")
# # its uglier and there are more outliers, but it generally looks like blank and real treatments are doing the same things
# # in the summer while plots with real sponges maintain more species in the winter by the end of the experiment.