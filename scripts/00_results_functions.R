decimalplaces <- function(x) {
  if ((x - round(x)) != 0) {
    strs <- strsplit(as.character(format(x, scientific = F)), "\\.")
    n <- nchar(strs[[1]][2])
  } else {
    n <- 0
  }
  return(n) 
}

# summary tables
aov.res<-function(rtable,row){
  tv<-round(rtable[row,4],2)
  pv<-ifelse(rtable[row,5]< 0.001, "p < 0.001", paste("p =", signif(rtable[row,5],1)))
  return(paste("=",tv,",",pv))
}
aov.df<-function(rtable,row){
  return(round(rtable[row,3],0))
}

efsize<-function(rtable,row){
  betas<-signif(rtable$coefficients[row,1], 2)
  betas<-ifelse(abs(betas)<10,
                formatC(signif(betas,2), digits=2, format="fg", flag="#"), 
                round(betas))
  betas2 <- as.numeric(betas)
  decicount <- ifelse(betas2>1 & betas2 <10, 1, decimalplaces(betas2))
  se<-rtable$coefficients[row,2]
  cil<-rtable$coefficients[row,1]-1.96*se
  cil<-format(round(cil, digits=decicount), scientific=F)
  cih<-rtable$coefficients[row,1]+1.96*se
  cih<-format(round(cih, digits=decicount), scientific=F)
  return(paste0(betas,", ",cil," to ",cih))
}

efsize.tmb<-function(rtable,row){
  betas<-signif(rtable$coefficients$cond[row,1], 2)
  betas<-ifelse(abs(betas)<10,
                formatC(signif(betas,2), digits=2, format="fg", flag="#"), 
                round(betas))
  betas2 <- as.numeric(betas)
  decicount <- ifelse(betas2>1 & betas2 <10, 1, decimalplaces(betas2))
  se<-rtable$coefficients$cond[row,2]
  cil<-rtable$coefficients$cond[row,1]-1.96*se
  cil<-format(round(cil, digits=decicount), scientific=F)
  cih<-rtable$coefficients$cond[row,1]+1.96*se
  cih<-format(round(cih, digits=decicount), scientific=F)
  return(paste0(betas,", ",cil," to ",cih))
}

# overall anovas
aov.res2<-function(rtable,row){
  tv<-round(rtable[row,5],2)
  pv<-ifelse(rtable[row,6]< 0.001, "p < 0.001", paste("p =", round(rtable[row,6],2)))
  return(paste("=",tv,",",pv))
}
aov.df2<-function(rtable,row){
  return(paste0(round(rtable[row,3],0),",",round(rtable[row,4],0)))
}

# tmb models
# overall anovas
tmb.res<-function(rtable,row){
  tv<-round(rtable$coefficients$cond[row,3],2)
  pv<-ifelse(rtable$coefficients$cond[row,4]< 0.001, "p < 0.001", paste("p =", round(rtable$coefficients$cond[row,4],2)))
  return(paste("=",tv,",",pv))
}

# degrees of freedom??
