### R script supplied with Pirana
### by Ron Keizer, 2010
###
### Required: - NM table file with DV and TIME on first $TABLE record
###           - Lattice library
###
### Description: This R-script create a plot of DV+PRED versus TIME,
### split by individuals, for multiple selected models
###

library(lattice)
setwd('C:/Users/yfzhang/Dropbox/Programming/Pirana/Metabolism of Arctigenin')

runnum<-'133'
path<-paste('C:/Users/yfzhang/Dropbox/Programming/Pirana/Metabolism of Arctigenin/sdtab',runnum,'.tab',sep='')

tab<- read.table (path, skip=1, header=T) # NONMEM table with ONEHEADER option
tab1<-tab[tab$DV>0,]
tabo<-tab1[order(tab1$TIME),]

xyplot (DV+PRED~TIME | factor(CMT,labels=c('AG','AA'))+factor(DL,labels=c('2.4 mg/kg','4.8 mg/kg','12 mg/kg')), data=tabo,type=c("p","l"),
        main = paste (runnum,": ", "AR Metabolism", sep=""),
        auto.key=T, as.table=T, pch=19,distribute.type=T,
        xlab="Time", ylab="Dependent variable / Pop. prediction"
)

