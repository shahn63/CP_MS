## load packages
library("readxl")

## load data
## LC_MRM data example

info=read_excel("sample_info.xlsx")
dat=read_excel("MRMdata.xlsx")

## data pre-processing

# log transform & remove ID column
dat[,-1]=log(dat[,-1]+1)
dat<-dat[,-1]

# use only non-zero mz-values
idx4<-which(apply(dat,2,sum)==0)
p<-12000-length(idx4) # dimension
dat<-dat[,-idx4]

# split the data
# abnormal test dataset
idx=which(info$Type=="Processed Blanks")
test<-dat[idx,]

# normal (training, calibration, test) dataset
dat<-dat[-idx,]

# random split of normal data
samp<-sample(1:dim(dat)[1],dim(dat)[1], replace=F)
tr.idx<-samp[1:1800] # training set index
cal.idx<-samp[1801:3600] # calibration set index
cov.idx<-samp[3601:4655] # normal test set index

tr.dat<-dat[tr.idx,]
cal.dat<-dat[cal.idx,]
cov.dat<-dat[cov.idx,]

## load CPMS function
source("CPMS.R")

## load CPMS_band function
source("CPMS_band.R")

res<-CPMS(tr.dat=tr.dat,cal.dat=cal.dat,test=cov.dat,alpha=0.05,G=3,samp.plot=T)
summary(res)

res.peak<-CPMS_band(CPMS.obj=res,test=cov.dat,id=207,peak.plot=T)
res.peak$out.peakID

