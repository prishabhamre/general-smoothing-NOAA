NOAA1<-NOAAnew

plot(NOAA1[,3],NOAA1[,2],xlab="temperature rise",ylab="rate of billion dollar
weather disasters")

dumb<-bin.mean(NOAA1[,3],NOAA1[,2],6)
dumg<-gauss.mean(NOAA1[,3],NOAA1[,2],.063)
gauss.reg(NOAA1[,3],NOAA1[,2],.078,do.plot=T)
gauss.mean.trunc(NOAA1[,3],NOAA1[,2],.063,20,do.plot=T)
gauss.reg.trunc(NOAA1[,3],NOAA1[,2],.08,17,do.plot=T)
lines(lowess(NOAA1[,3],NOAA1[,2]),col=7)
lines(smooth.spline(NOAA1[,3],NOAA1[,2]),col=8)
(smooth.spline(NOAA1[,3],NOAA1[,2]))$df

bin<-greedybin(bin.mean,NOAA1,3,3,6)
mean<-greedy(gauss.mean,NOAA1,3,2,.063)
reg<-greedy(gauss.reg,NOAA1,3,2,.078)
meantrunc<-greedytrunc(gauss.mean.trunc,NOAA1,3,2,.063,20)
regtrunc<-greedytrunc(gauss.reg.trunc,NOAA1,3,2,.08,17)

newb<-bin.mean(NOAA1[,3],NOAA1[,2],bin$nbin)
newg<-gauss.mean(NOAA1[,3],NOAA1[,2],mean$theta)
gauss.reg(NOAA1[,3],NOAA1[,2],reg$theta,do.plot=T)
gauss.mean.trunc(NOAA1[,3],NOAA1[,2],meantrunc$theta,meantrunc$nnn,do.plot=T)
gauss.reg.trunc(NOAA1[,3],NOAA1[,2],regtrunc$theta,regtrunc$nnn,do.plot=T)


