
CPMS_band<-function(CPMS.obj=CPMS_result, test=test_dat, id=detected_sample_ID, peak.plot=T){
	G<-CPMS.obj$numG
	ff.idx<-1:G
	Lt<-CPMS.obj$mag.lower
	Ut<-CPMS.obj$mag.upper

	pow.out<-pow.outk<-NULL
	for(i in 1:ncol(test)){
		pow.outk<-0
		for(j in 1:length(ff.idx)){
			if(test[id,i] < Lt[j,i] | test[id,i] > Ut[j,i]){
				pow.outk<-pow.outk+1
			}
		}
		if(pow.outk==length(ff.idx)){
			pow.out<-union(pow.out,i)
		}
	}

	if(peak.plot==T){
	linet<-c(1:length(ff.idx)+1,1)
	layout(matrix(1:(length(ff.idx)+1), ncol = 1, byrow = TRUE), heights = c(2,2,2,1.7))

	colidx<-(1:G)+2
	for(i in 1:length(ff.idx)){
		plot(1:ncol(test),test[id,],type="p",pch=16,cex=0.9,col=1,ylim=c(0,max(test,Ut)),xlim=c(0,ncol(test)),lty=linet[length(linet)],
		xlab=" ", ylab=" ");par(new=T)
		plot(1:ncol(test),Lt[i,],type="l",ylim=c(0,max(test,Ut)),col=colidx[i],lwd=1,lty=1,xlim=c(0,ncol(test)),xlab="mz values", ylab="log(intensity+1)");par(new=T)
		lines(1:ncol(test),Ut[i,],type="l",col=colidx[i],lwd=1,xlim=c(0,ncol(test)),lty=1)
		#abline(v=pow.out,col=1,lty=2)
		title(paste("Conformal prediction band Tn",sep="",i," for simulated data (alpha=",alpha,")"))
		legend("topright",c("Outlier sample",paste("Conformal prediction band Tn",sep="",i)),
		col=c(1,colidx[i]),pch=c(16,NA),lty=c(NA,1),lwd=c(0,3))
	}
	
	## Make y-value=0
	pow.out<-data.frame(c(0,pow.out,dim(test)[2]),0)
	## Plotting without box or axis with dot, representing data points                  
	plot(pow.out,bty='n',xaxt='n',yaxt='n',ylab='',xlab='',pch=16,cex=c(0,rep(1,nrow(pow.out)-2),0))
	for (i in seq(1, dim(pow.out)[1], 1)) {
		lines(c(1, dim(test)[2]), c(0, 0), lwd = 0.1)
	}

	## Placing x-values & x-axis label onto plot
	text(pow.out[-c(1,nrow(pow.out)),],labels=pow.out[-c(1,nrow(pow.out)),1],pos=3,offset=1,font=0.5,cex=0.6,srt=90)
	title("Outlier Sample ID",line=-3)
	}

	return(list(out.peakID=pow.out[-1,1]))
}


