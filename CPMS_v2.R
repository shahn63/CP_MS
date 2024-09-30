############################################
# Conformal Prediction 
# for outlier detection in MS data
############################################

CPMS<-function(tr.dat=train_data, cal.dat=calibration_data, test=test_data, alpha=0.05, G=3, samp.plot=T){

# mixture model parameter estimation using training data by K-means

# K-means clustering at a given G
ff<-kmeans(tr.dat,centers=G,iter.max=10000)
while(sum(table(ff$cluster)==1)!=0){
	ff<-kmeans(tr.dat,centers=(length(unique(ff$cluster))-sum(table(ff$cluster)==1)),iter.max=10000)
}

ff.idx<-1:length(ff$size) #G
phi<-ff$size/sum(ff$size)
mu<-t(ff$centers[ff.idx,])
if(length(ff.idx)==1){mu<-t(mu)}
sd<-array(0,dim=c(ncol(tr.dat),ncol(tr.dat),length(ff.idx)))

if(length(ff.idx)==1){
	sd[,,1]<- var(tr.dat)
}else{ 
	for(i in 1:length(ff.idx)){
		sd[,,i]<-var(tr.dat[which(ff$cl==i),])
	}
}

# positivie-definite for high-dimensional data
for(i in 1:length(ff.idx)){
	sd[,,i]<-(sd[,,i]+t(sd[,,i]))/2+diag(ncol(tr.dat))*0.001
}

# calculate MH distance/ confomaloty score using calibration data
mh<-matrix(0,nrow=length(ff.idx),ncol=dim(cal.dat)[1])
for(j in 1:length(ff.idx)){
	mh[j,]<-mahalanobis(cal.dat, mu[,j], sd[,,j])
}

conf<-NULL
for(i in 1:length(ff.idx)){
	conf<-cbind(conf,mh[i,]-2*log(phi[i])+determinant(sd[,,i])$modulus)
}
conf.score<-apply(conf,1,min)
ord.mh<-order(conf.score,decreasing=T)

idx.a<-floor((nrow(cal.dat)+1)*alpha)-1
lambda<-conf.score[ord.mh][idx.a] # threshold for conformality score

# outlier detection in the test data

mh.pow<-matrix(0,nrow=length(ff.idx),ncol=dim(test)[1])
for(j in 1:length(ff.idx)){
	mh.pow[j,]<-mahalanobis(test, mu[,j], sd[,,j])
}

Tn<-list()
TT<-NULL
for(i in 1:length(ff.idx)){
	Tn[[i]]<-which(mh.pow[i,] <= lambda+2*log(phi[i])-determinant(sd[,,i])$modulus)
	TT<-union(TT,Tn[[i]])
}
if(length(TT)!=0){
	sampID<-(1:nrow(test))[-TT]
}else{
	sampID<-(1:nrow(test))
}

## conformal score/MH distance plot for outlier sample detection

if(samp.plot==T){
	par(mfrow=c((length(ff.idx) + 1),1), mar=c(0,4,2,1), oma=c(0,0,1,0))

	for(i in 1:length(ff.idx)){
		cutoff<-lambda+2*log(phi[i])-determinant(sd[,,i])$modulus
		plot(1:dim(test)[1],log(mh.pow[i,]),pch=1,cex=0.8, ylim=c(0,max(log(mh.pow),log(cutoff))),
		xlab="Sample ID", ylab=paste("log(MH)"," of cluster",i))
		abline(h=log(cutoff),col=2,lwd=2)
		pow.out<-sampID #c(1:dim(test)[1])[-sort(TT)]	
		abline(v=pow.out,col=1,lty=2)
	}

	pow.out<-c(0,sampID,dim(test)[1])
	## Make y-value=0
	pow.out<-data.frame(pow.out,0)
	## Plotting without box or axis with dot, representing data points                  
	plot(pow.out,bty='n',xaxt='n',yaxt='n',ylab='',xlab='',pch=16,cex=c(0,rep(1,nrow(pow.out)-2),0))
	#axis(side=1,seq(1,dim(test)[1],1),pos=0,labels=F)  ## Using y-value as position 
	for (i in seq(1, dim(test)[1], 1)) {
	lines(c(1, dim(test)[1]), c(0, 0), lwd = 0.1)
	}

	## Placing x-values & x-axis label onto plot
	text(pow.out[-c(1,nrow(pow.out)),],labels=pow.out[-c(1,nrow(pow.out)),1],pos=3,offset=1,font=0.5,cex=0.6,srt=90)
	title("Outlier Sample ID",line=-3)
}


# marginal prediction intervals

uni.conf<-matrix(0,nrow=dim(cal.dat)[1],ncol=dim(cal.dat)[2])
for(i in 1:dim(cal.dat)[2]){
	mag.mh<-NULL
	for(j in 1:length(ff.idx)){
		mh.ij<-mahalanobis(as.matrix(cal.dat[,i]), mu[i,j], sd[i,i,j])
		mag.mh<-cbind(mag.mh, mh.ij-2*log(phi[j])+log(sd[i,i,j]))
	}
	uni.conf[,i]<-apply(mag.mh,1,min)
}

idx.a<-floor((nrow(cal.dat)+1)*alpha)-1
uni.lambda<-NULL
for(i in 1:dim(cal.dat)[2]){
	ord.mh<-order(uni.conf[,i],decreasing=T)
	uni.lambda[i]<-uni.conf[ord.mh,i][idx.a]
}

Lt<-Ut<-matrix(0,nrow=length(ff.idx),ncol=ncol(cal.dat))
for(j in 1:length(ff.idx)){
	for(i in 1:ncol(cal.dat)){
		if((2*log(phi[j])+uni.lambda[i]-log(sd[i,i,j]))<0){
			Ut[j,i]<-mu[i,j]+qnorm(1-alpha/2)*sqrt(sd[i,i,j]/nrow(cal.dat))
			Lt[j,i]<-mu[i,j]-qnorm(1-alpha/2)*sqrt(sd[i,i,j]/nrow(cal.dat))
		}else{
		Ut[j,i]<-mu[i,j]+sqrt((2*log(phi[j])+uni.lambda[i]-log(sd[i,i,j]))*sd[i,i,j])
		Lt[j,i]<-mu[i,j]-sqrt((2*log(phi[j])+uni.lambda[i]-log(sd[i,i,j]))*sd[i,i,j])
		}
	}
}

return(list(numG=G, mu=mu, sd=sd, phi=phi, lmabda=lambda, out.samp=sampID, mag.lower=Lt, mag.upper=Ut))
}


