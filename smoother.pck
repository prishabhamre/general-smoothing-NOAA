smoother.pck <-
c("smoother.pck", "bin.mean", "gauss.mean", "gauss.reg", "gauss.mean.trunc", 
"gauss.reg.trunc", "my.hat.w")
bin.mean <-
function(x,y,nbin,xcol=2)
{
o1<-order(x)
x1<-x[o1]
y1<-y[o1]
r1<-range(x)
inc<-(r1[2]-r1[1])/nbin
yvec<-NULL
smat<-NULL
for(i in 1:nbin){
        bin.low<-r1[1]+(i-1)*inc
        bin.high<-r1[1]+i*inc
        
        I1<-x1>=bin.low
if(i<nbin){
        I2<-x1<bin.high
}else{
        I2<-x1<=(bin.high+200)
}
        I3<-as.logical(I1*I2)
        yval<-mean(y1[I3])
        n1<-sum(I3)
        matdum<-NULL
        for(i in 1:n1){
        matdum<-rbind(matdum,I3*1/n1)
        }
        smat<-rbind(smat,matdum)
        yvec<-c(yvec,rep(yval,n1))
}
n99<-length(x1)
dferror<-length(x1)-sum(diag(2*smat-smat%*%(t(smat))))
delta1<-sum(diag(t(diag(n99)-smat)%*%(diag(n99)-smat)))
R<-t(diag(n99)-smat)%*%(diag(n99)-smat)
delta2<-2*sum(diag(R%*%R))
lines(x1,yvec,col=xcol)
ypred<-y1
ypred<-smat%*%y1
resid<-y-ypred
PRESS<-sum((resid/(1-diag(smat)))^2)
list(smat=smat,df=sum(diag(smat)),dferror=dferror,delta1=delta1,delta2=delta2,resid=resid,pred=ypred,x=x,press=PRESS)
                
}
gauss.mean <-
function(x,y,lambda,xcol=3,do.plot=T)
{
o1<-order(x)
x1<-x[o1]
y1<-y[o1]
r1<-range(x)
smat<-NULL
n1<-length(x1)
for(i in 1:n1){
        v1<-dnorm(x1,x1[i],lambda)
        v1<-v1/sum(v1)
        smat<-rbind(smat,v1)
}
yhat<-smat%*%y1
if(do.plot){
lines(x1,yhat,col=xcol)
}
n99<-length(x1)
dferror<-length(x1)-sum(diag(2*smat-smat%*%(t(smat))))
delta1<-sum(diag(t(diag(n99)-smat)%*%(diag(n99)-smat)))
R<-t(diag(n99)-smat)%*%(diag(n99)-smat)
delta2<-2*sum(diag(R%*%R))
resid<-y1-smat%*%y1
ypred<-y1
ypred[o1]<-smat%*%y1
PRESS<-sum((resid/(1-diag(smat)))^2)
list(smat=smat,df=sum(diag(smat)),dferror=dferror,delta1=delta1,delta2=delta2,resid=resid,pred=ypred,press=PRESS)
        
}
gauss.reg <-
function(x,y,lambda,xcol=4,do.plot=T)
{
o1<-order(x)
x1<-x[o1]
y1<-y[o1]
r1<-range(x)
smat<-NULL
n1<-length(x1)
for(i in 1:n1){
        v1<-dnorm(x1,x1[i],lambda)
        v1<-v1/sum(v1)
        H1<-my.hat.w(x1,v1)
        smat<-rbind(smat,H1[i,])
}
yhat<-smat%*%y1
if(do.plot){
lines(x1,yhat,col=xcol)
}
n99<-length(x1)
dferror<-length(x1)-sum(diag(2*smat-smat%*%(t(smat))))
delta1<-sum(diag(t(diag(n99)-smat)%*%(diag(n99)-smat)))
R<-t(diag(n99)-smat)%*%(diag(n99)-smat)
delta2<-2*sum(diag(R%*%R))
resid<-y1-smat%*%y1
ypred<-y1
ypred[o1]<-smat%*%y1
PRESS<-sum((resid/(1-diag(smat)))^2)
list(smat=smat,df=sum(diag(smat)),dferror=dferror,delta1=delta1,delta2=delta2,resid=resid,pred=ypred,press=PRESS)
}
gauss.mean.trunc <-
function(x,y,lambda,nnn,xcol=5,do.plot=T)
{
o1<-order(x)
x1<-x[o1]
y1<-y[o1]
r1<-range(x)
smat<-NULL
n1<-length(x1)
trunc.val<-n1-nnn
for(i in 1:n1){
        v1<-dnorm(x1,x1[i],lambda)
        o2<-order(v1)
        thresh<-v1[o2[trunc.val]]
        v1<-v1*(v1>thresh)
        v1<-v1/sum(v1)
        smat<-rbind(smat,v1)
}
yhat<-smat%*%y1
if(do.plot){
lines(x1,yhat,col=xcol)
}
n99<-length(x1)
dferror<-length(x1)-sum(diag(2*smat-smat%*%(t(smat))))
delta1<-sum(diag(t(diag(n99)-smat)%*%(diag(n99)-smat)))
R<-t(diag(n99)-smat)%*%(diag(n99)-smat)
delta2<-2*sum(diag(R%*%R))
resid<-y1-smat%*%y1
ypred<-y1
ypred[o1]<-smat%*%y1
PRESS<-sum((resid/(1-diag(smat)))^2)
list(smat=smat,df=sum(diag(smat)),dferror=dferror,delta1=delta1,delta2=delta2,resid=resid,pred=ypred,press=PRESS)
                
}
gauss.reg.trunc <-
function(x,y,lambda,nnn,xcol=6,do.plot=T)
{
o1<-order(x)
x1<-x[o1]
y1<-y[o1]
r1<-range(x)
smat<-NULL
n1<-length(x1)
trunc.val<-n1-nnn
for(i in 1:n1){
        v1<-dnorm(x1,x1[i],lambda)
        o1<-order(v1)
        thresh<-v1[o1[trunc.val]]
        v1<-v1*(v1>thresh)
        v1<-v1/sum(v1)
        H1<-my.hat.w(x1,v1)
        smat<-rbind(smat,H1[i,])
}
yhat<-smat%*%y1
if(do.plot){
lines(x1,yhat,col=xcol)
}
n99<-length(x1)
dferror<-length(x1)-sum(diag(2*smat-smat%*%(t(smat))))
delta1<-sum(diag(t(diag(n99)-smat)%*%(diag(n99)-smat)))
R<-t(diag(n99)-smat)%*%(diag(n99)-smat)
delta2<-2*sum(diag(R%*%R))
resid<-y1-smat%*%y1
ypred<-y1
ypred[o1]<-smat%*%y1
PRESS<-sum((resid/(1-diag(smat)))^2)
list(smat=smat,df=sum(diag(smat)),dferror=dferror,delta1=delta1,delta2=delta2,resid=resid,pred=ypred,press=PRESS)

        
}
my.hat.w <-
function(x,wt){
x1<-cbind(1,x)
x1%*%solve(t(x1)%*%diag(wt)%*%x1)%*%t(x1)%*%(diag(wt))
}

greedy<-function(func, data, xcol, ycol, theta){
  #press with original paramaters
  press0<- func(data[,xcol],data[,ycol],theta,F)$press
  #stores original press
  press00<-press0
  inc<- 0
  #stores original parameter
  theta0<-theta
  #loop runs up to 100 times with lowest current press val
  while(inc<100){
    epsilon<- rnorm(1,0,.01)
    #random value to the parameter lambda
    theta1 <- (theta+epsilon)
    #recalculate press statistic and store
    press1<- func(data[,xcol],data[,ycol],theta1)$press
    if (press1<press0){
      press0=press1
      inc=0
      theta=theta1
    }
    else{
      inc<-inc+1
    }
  }
  list(theta = theta, theta0=theta0, press=press0, press0=press00)
  print(list(theta = theta, theta0=theta0, press=press0, press0=press00))
}

greedytrunc<- function(func, data, xcol, ycol, theta,nnn){
  #press with original paramaters
  press0<- func(data[,xcol],data[,ycol],theta,nnn,do.plot=F)$press
  #stores original press
  press00<-press0
  #stores original nnn value
  nnn0<-nnn
  inc<- 0
  #stores original parameter
  theta0<-theta
  #loop runs up to 100 times with lowest current press val
  while(inc<100){
    epsilon<- rnorm(1,0,.01)
    epsilon1<- sample(c(-1,0,1),1)
    #add random value epsilon to the parameter lambda
    theta1 <- (theta+epsilon)
    #add random value epsilon1 to parameter nnn
    nnn1<- nnn+epsilon1
    #recalculate press statistic and store
    press1<- func(data[,xcol],data[,ycol],theta1,nnn1,do.plot=F)$press
    if (press1<press0){
      press0=press1
      inc=0
      theta=theta1
      nnn=nnn1
    }
    else{
      inc<-inc+1
    }
  }
  list(theta = theta, theta0=theta0,nnn0=nnn0,nnn=nnn, press=press0, press0=press00)
  print(list(theta = theta, theta0=theta0,nnn0=nnn0,nnn=nnn, press=press0, press0=press00))
}

greedybin<- function(func, data, xcol, ycol,nbin){
  #press with original paramaters
  press0<- func(data[,xcol],data[,ycol],nbin)$press
  #stores original press
  press00<-press0
  #stores original nbin value
  nbin0<-nbin
  nbin00<-nbin
  inc<- 0
  #loop runs up to 100 times with lowest current press val
  while(inc<100){
    epsilon<- sample(-2:2,1)
    #add random value epsilon to parameter nbin
    nbin1<- nbin0+epsilon
    while(nbin1>9||nbin1<1){
      epsilon<-sample(-2:2,1)
      nbin1<- nbin0+epsilon
    }
    #recalculate press statistic and store
    press1<- func(data[,xcol],data[,ycol],nbin1)$press
    if (press1<press0){
      press0=press1
      nbin0=nbin1
    }
    else{
      inc<-inc+1
    }
  }
  list(nbin0=nbin00, nbin=nbin0, press=press0,press0=press00)
  print(list(nbin0=nbin00, nbin=nbin0, press=press0,press0=press00))
}

