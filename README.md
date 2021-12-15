# nonparam-predict-mpower2_61-yahoo.com

rm(list=ls())
library(limSolve)
library(forecast)
library(stable)
m=1000;n1=100;h<-5;H1=.9;G1=1;alpha=1.5
c1=matrix(0,nrow=n1,ncol=2*m)
for (t in 1:n1){
for (j1 in 1:m){
#c1[n1-t+1,j1]=cos(n1*pi*j1/m)*cos(pi*j1*(t-1)/m)+sin(n1*pi*j1/m)*sin(pi*j1*(t-1)/m)
#c1[n1-t+1,j1+m-1]=sin(n1*pi*j1/m)*cos((t-1)*pi*j1/m)-sin((t-1)*pi*j1/m)*cos(n1*pi*j1/m)
c1[t,j1]=2*cos(t/H1*pi*j1/m)/G1
c1[t,j1+m]=-2*sin(t/H1*pi*j1/m)/G1
}
}
det(c1)
#v1=seq(50,200,50);v2=length(v1)
#c11=matrix(0,nrow=n1-v2,ncol=2*m-v2)
#c11=c1[-v1,-v1]
n=n1+h;it=1000
pr1=matrix(0,nrow=it,ncol=h);pr2=matrix(0,nrow=it,ncol=h)
pr3=matrix(0,nrow=it,ncol=h);pr11=matrix(0,nrow=it,ncol=h)
pr22=matrix(0,nrow=it,ncol=h);pr111=matrix(0,nrow=it,ncol=h)
pr222=matrix(0,nrow=it,ncol=h);pr33=matrix(0,nrow=it,ncol=h)
pr333=matrix(0,nrow=it,ncol=h);pr4=matrix(0,nrow=it,ncol=h)
pr44=matrix(0,nrow=it,ncol=h);pr444=matrix(0,nrow=it,ncol=h)
pr5=matrix(0,nrow=it,ncol=h);pr55=matrix(0,nrow=it,ncol=h)
pr555=matrix(0,nrow=it,ncol=h)
phi0=c(0.7, -0.3);theta=c(-0.2, 0.1)
for (k in 1:it){	
#Gaussian
#x0=arima.sim(n=n1+h,list(ar = c(0.7, -0.3), ma = c(-0.2, 0.1)),
#innov=rnorm(n1+h,0,1))
#stable
x0=0;x0[1]=0;x0[2]=0;z1=rstable(n1+h,alpha,0,1,0,1)
for (i0 in 3:(n1+h)){
x0[i0]=phi0[1]*x0[i0-1]+phi0[2]*x0[i0-2]+z1[i0]+theta[1]*z1[i0-1]+theta[2]*z1[i0-2]	
}
################
x1=x0[c(1:n1)]
A=Solve(c1,x1);w2=0
for (i1 in 1:h){
w1=0;t1=n1+i1
for (j in 1:m){
#w1[j]=(cos(n1*pi*j/m)*cos(pi*j*i1/m)-sin(n1*pi*j/m)*sin(pi*j*i1/m))*A[j]
w1[j]=2*cos(t1/H1*pi*j/m)/G1*A[j]
#w1[j+m]=(sin(n1*pi*j/m)*cos(i1*pi*j/m)+sin(i1*pi*j/m)*cos(n1*pi*j/m))*A[j+m]
w1[j+m]=-2*sin(t1/H1*pi*j/m)/G1*A[j+m]
}
w2[i1]=sum(w1)
pr1[k,i1]=sum(w1)-x0[n1+i1]
pr11[k,i1]=100*abs(pr1[k,i1]/x0[n1+i1])
pr111[k,i1]=pr1[k,i1]/mean(abs(diff(x1)))
}
#################
#Gaussian
#fit <- arima(x1, order = c(2,0,2),method=c("CSS"))
#q1=predict(fit,n.ahead=h)
#q2=q1$pred
#pr2[k,]=q2-x0[c((n1+1):(n1+h))]
#####################
#stable
phi=matrix(0,nrow=4,ncol=4)
phi[1,]=c(phi0,theta)
phi[2,1]=1;phi[4,3]=1
yt=matrix(0,nrow=4,ncol=1)
yt[,1]=c(x0[n1],x0[n1-1],z1[n1],z1[n1-1])
e11=matrix(c(1,0,0,0),nrow=4,ncol=1)
for (h1 in 1:h){
pr2[k,h1]=t(e11)%*%(phi%^%h1)%*%yt-x0[n1+h1]
}
#####################
fit1 <- nnetar(x1, size=2)
f1=forecast(fit1, h=5)
pr3[k,]=f1$mean-x0[c((n1+1):(n1+h))]

#s=ssa(x1)
#f2 <- forecast(s, groups = list(1:6), method = "recurrent", bootstrap = TRUE, len = 10, R = 10)
#pr4[k,]=f2$mean-x0[c((n1+1):(n1+h))]
#pr44[k,i1]=100*pr4[k,i1]/x0[n1+i1]
#pr444[k,i1]=pr4[k,i1]/mean(abs(diff(x1)))
###########################################
###########################################
M=floor(2*n1/3)
rmseh=matrix(0,nrow=M/2-1,ncol=M/2-1)
for (l in 2:(M/2)){
for (r in 1:(l-1)){
k1=n1-l+1
xx=matrix(0,nrow=l,ncol=k1)
for (i in 1:k1){
xx[,i]=x0[i:(l+i-1)]
}
ei=eigen(xx%*%t(xx))
ev=ei$values
eve=ei$vectors
pim=matrix(0,nrow=l-1,ncol=r)
j=1:r
pim[,j]=eve[c(1:(l-1)),j]
pim1=0
j1=1:r
pim1[j1]=eve[l,j1]
nu2=sum(pim1^2)
R1=rowSums(nu2*pim)/(1-nu2)
R2=0
i2=1:(l-1)
R2=R1[l-1-i2+1];rmse=0;x2=0
for (i3 in 1:(n1-M-h+1)){
#####################
#Reconstructed time series
#for (i4 in i3:(M+i3-1)){
#xx0=matrix(0,nrow=l,ncol=k1);xx1=matrix(0,nrow=l,ncol=k1)
#for (y1 in 1:l){
#for (y2 in 1:k1){
#if (y1+y2==i4+1) xx0[y1,y2]=xx[y1,y2] else xx0[y1,y2]=0
#if (y1+y2==i4+1) xx1[y1,y2]=1 else xx1[y1,y2]=0
#}
#}
#x2[i4]=sum(xx0)/sum(xx1)
#}
##################
#Initial time series
x2[i3:(M+i3-1)]=x0[i3:(M+i3-1)]
for (h1 in 1:h){
e1=0
j3=1:(l-1)
e1[j3]=R2[j3]*x2[M+i3+h1-1-j3]
x2[M+i3+h1-1]=sum(e1)
}
rmse[i3]=sqrt(mean((x2[(M+i3):(M+h+i3-1)]-x0[(M+i3):(M+h+i3-1)])^2))
}
rmseh[l-1,r]=mean(rmse)
}
}

rmseh1=matrix(0,nrow=M/2-1,ncol=M/2-1)
for (l in 2:(M/2)){
for (r in 1:(M/2-1)){
if (rmseh[l-1,r]==0) rmseh1[l-1,r]=max(rmseh) else rmseh1[l-1,r]=rmseh[l-1,r]
}} 
d2=which(rmseh1 == min(rmseh1), arr.ind = TRUE)
l=d2[1,1]+1
r=d2[1,2]
#############################################
xx=matrix(0,nrow=l,ncol=k1)
for (i in 1:k1){
xx[,i]=x0[i:(l+i-1)]
}
ei=eigen(xx%*%t(xx))
ev=ei$values
eve=ei$vectors
pim=matrix(0,nrow=l-1,ncol=r)
j=1:r
pim[,j]=eve[c(1:(l-1)),j]
pim1=0
j1=1:r
pim1[j1]=eve[l,j1]
nu2=sum(pim1^2)
R1=rowSums(nu2*pim)/(1-nu2)
R2=0
i2=1:(l-1)
R2=R1[l-1-i2+1]
x2[1:n1]=x0[1:n1]
for (h1 in 1:h){
e1=0
j3=1:(l-1)
e1[j3]=R2[j3]*x2[n1-j3+1+h1-1]
x2[n1+h1]=sum(e1)
}
pr4[k,]=x2[(n1+1):(n1+h)]-x0[(n1+1):(n1+h)]

for (i11 in 1:h){
pr22[k,i11]=100*abs(pr2[k,i11]/x0[n1+i11])
pr222[k,i11]=pr2[k,i11]/mean(abs(diff(x1)))
pr33[k,i11]=100*abs(pr3[k,i11]/x0[n1+i11])
pr333[k,i11]=pr3[k,i11]/mean(abs(diff(x1)))
pr44[k,i11]=100*abs(pr4[k,i11]/x0[n1+i11])
pr444[k,i11]=pr4[k,i11]/mean(abs(diff(x1)))
}
###########################################
###########################################
}
