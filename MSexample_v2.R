## load packages
library("readxl")

## load CPMS function
source("CPMS_v2.R")

## load CPMS_band function
source("CPMS_band.R")


## load data
## LC_MRM data example

#dat=read.csv("open_dat.csv")

## data pre-processing: log transform & remove ID column
dat[,-1]=log(dat[,-1]+1)
dat<-dat[,-1]

p<-ncol(dat)

sum(is.na(dat))
rawdat<-dat<-na.omit(dat)
dim(rawdat)

# generate virtual abnormal test samples
set.seed(123)
ts.idx<-sample(1:nrow(dat),50,replace=F)
test<-dat[ts.idx,]
ts.idx2<-matrix(0,50,10)
for(i in 1:50){
	set.seed(i)
	ts.idx2[i,]<-sample(1:ncol(dat),10,replace=F)
	test[i,ts.idx2[i,]]<-test[i,ts.idx2[i,]]+rnorm(10,1,0.5)
}
dim(test)


# random split the data (train, calibration, normal test sets)
# sample size: 160 / 160/ 57

set.seed(600)
samp<-sample(1:dim(dat)[1],dim(dat)[1], replace=F)
tr.idx<-samp[1:160]
cal.idx<-samp[161:320]
cov.idx<-samp[321:nrow(dat)]

tr.dat<-dat[tr.idx,]
cal.dat<-dat[cal.idx,]
cov.dat<-dat[cov.idx,]
dim(tr.dat);dim(cal.dat);dim(cov.dat)

## outlier sample detection
set.seed(202408)
res<-CPMS(tr.dat=tr.dat,cal.dat=cal.dat,test=rbind(cov.dat,test),alpha=0.05,G=3,samp.plot=T)
summary(res)

## outlier peak detection
res.peak<-CPMS_band(CPMS.obj=res,test=rbind(cov.dat,test),id=107,peak.plot=T)
res.peak$out.peakID

