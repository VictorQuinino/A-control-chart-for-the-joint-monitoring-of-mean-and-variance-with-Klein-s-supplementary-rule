library(pracma)
clear()
tic()
# Definition of parameters
#Target ARL0=250. 
#In-control
u0=0
sd0=1
#out-of-control - Table 3. 
u1=0
sd1=1
lambda1=0.1
lambda2=0.1
#lambda1=0.8
#lambda2=0.8
#Chen et. al. (2001)
k1=2.81
k2=2.86
#k1=3.09
#k2=3.88

n=5
LCL1=u0-k1*((lambda1/(2-lambda1))^0.5)*(sd0/(n^0.5))
UCL1=u0+k1*((lambda1/(2-lambda1))^0.5)*(sd0/(n^0.5))
c=log(sd0^2)-1/(n-1)-1/(3*(n-1)^2)+2/(15*(n-1)^4)
d=2/(n-1)+2/((n-1)^2)+4/(3*(n-1)^3)-16/(15*(n-1)^5)
LCL2=c-k2*(d*lambda2/(2-lambda2))^0.5
UCL2=c+k2*(d*lambda2/(2-lambda2))^0.5

runs=1000000
Result<-matrix(0, runs, 1)

for(i in 1:runs){
  ua=0 #starting value
  sa=c #starting value
  s=0
  D<-c()
  s<-0
  while(s<1){
    RR<-rnorm(n,u1,sd1)
    R=mean(RR)
    V=var(RR)
    T=lambda1*R+(1-lambda1)*ua
    W=lambda2*log(V)+(1-lambda2)*sa
    
    if (T>UCL1 | T<LCL1 |W>UCL2 |W<LCL2){
      D<-rbind(D,1)
      s=s+1
    }else{
      D<-rbind(D,1)
      s=0
    }
    ua=T
    sa=W
  }
  Result[i,1]=sum(D)
  
  
}

ARL<-mean(Result[,1])
cat("ARL1=",ARL,"\n")
SDRL<-(var(Result[,1]))^0.5
cat("SDRL1=",SDRL,"\n")

toc()

