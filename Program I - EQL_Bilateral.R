library(cubature)
library(pracma)
tic()
clear ()
#In-control
u0 = 0
s0 = 1
sd0=s0
n = 5
#Klein Limits. Change as needed.
LCLxb = -0.916886 
UCLxb = 0.916886 
LCs2a=9.97299
LCs2b = 0.634607
#Xb-S2 Limits
LSCxb=1.37924
LICxb=-1.37924
LSCqui=18.5093
LICqui= 0.089926 
runs=10000 #Simulation EWMA
corridas=10000 #Simulation EQL
#EWMA 
lambda1=0.1
lambda2=0.1
lambda1=0.8
lambda2=0.8
#Chen et. al. (2001)
k1=2.81
k2=2.86
k1=3.09
k2=3.88
#EWMA Limits 
LCL1=u0-k1*((lambda1/(2-lambda1))^0.5)*(sd0/(n^0.5))
UCL1=u0+k1*((lambda1/(2-lambda1))^0.5)*(sd0/(n^0.5))
c=log(sd0^2)-1/(n-1)-1/(3*(n-1)^2)+2/(15*(n-1)^4)
d=2/(n-1)+2/((n-1)^2)+4/(3*(n-1)^3)-16/(15*(n-1)^5)
LCL2=c-k2*(d*lambda2/(2-lambda2))^0.5
UCL2=c+k2*(d*lambda2/(2-lambda2))^0.5


#Limits for Double Integrals
lowerLimit <- c(0, 0.25)  # Lower bound for u and s
upperLimit <- c(2, 2)  # Upper bound for u and s
#Uniform Distribution
R1=(1/(upperLimit[1]-lowerLimit[1]))
R2=(1/(upperLimit[2]-lowerLimit[2]))

# Function that calculates the EQL Klein
Inte_double <- function(params) {
  u=params[1]
  s=params[2]
  # Probability Calculation
  pxl<-pnorm(LCLxb,u,s/(n^0.5))
  pxu=1-pnorm(UCLxb,u,s/(n^0.5))
  pxc=1-pxu-pxl
  
  psu=1-pchisq(LCs2a*(s0^2/s^2),(n-1))
  psl=pchisq(LCs2b*(s0^2/s^2),(n-1))
  psc=1-psu-psl
  
  MarkovChain <- matrix(0, 25, 25)
  
  c1<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
  c2<- c(pxc*psc,0,pxc*psl,pxc*psu,0,pxu*psc,0,pxu*psl,pxu*psu,0,pxl*psc,0,pxl*psl,pxl*psu,0,0,0,0,0,0,0,0,0,0,0)
  c3<- c(pxc*psc,pxc*psu,0,0,pxc*psl,pxu*psc,pxu*psu,0,0,pxu*psl,pxl*psc,pxl*psu,0,0,pxl*psl,0,0,0,0,0,0,0,0,0,0)
  c4<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
  c5<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
  c6<- c(pxc*psc,pxc*psu,pxc*psl,0,0,0,0,0,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,0,0,0,0,0)
  c7<- c(pxc*psc,0,pxc*psl,pxc*psu,0,0,0,0,0,0,pxl*psc,0,pxl*psl,pxl*psu,0,pxu*psc,0,pxu*psl,pxu*psu,0,0,0,0,0,0)
  c8<- c(pxc*psc,pxc*psu,0,0,pxc*psl,0,0,0,0,0,pxl*psc,pxl*psu,0,0,pxl*psl,pxu*psc,pxu*psu,0,0,pxu*psl,0,0,0,0,0)
  c9<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
  c10<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
  c11<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,0,0,0,0,0,0,0,0,0,0,pxl*psc,pxl*psu,pxl*psl,0,0)
  c12<- c(pxc*psc,0,pxc*psl,pxc*psu,0,pxu*psc,0,pxu*psl,pxu*psu,0,0,0,0,0,0,0,0,0,0,0,pxl*psc,0,pxl*psl,pxl*psu,0)
  c13<- c(pxc*psc,pxc*psu,0,0,pxc*psl,pxu*psc,pxu*psu,0,0,pxu*psl,0,0,0,0,0,0,0,0,0,0,pxl*psc,pxl*psu,0,0,pxl*psl)
  c14<-c1
  c15<-c1
  c16<-c1
  c17<-c1
  c17<-c1
  c18<-c1
  c19<-c1
  c20<-c1
  c21<-c1
  c22<-c1
  c23<-c1
  c24<-c1
  c25<-c1
  
  MarkovChain <- rbind(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, 
                       c14, c15, c16, c17, c18, c19, c20, c21, c22, c23, c24,c25)
  
  A <- t(MarkovChain) - diag(25)
  A[25,] <- rep(1, 25)
  B <- rep(0, 25)
  B[25] <- 1
  
  Solved_markov_chain <- solve(A) %*% B
  ARL1 <- 1 / sum(Solved_markov_chain[c(4,5,9,10,14,15,16,17,18,19,20,21,22,
                                        23,24,25)])
  
  P2=(u^2+abs(s^2-1))*ARL1
  
  return(P2)
}
# Function that calculates the EQL Xb-S2
Inte <- function(params1){
  u=params1[1]
  s=params1[2]
  pc=pnorm(LICxb,u,s/(n^0.5))
  pa=1-pnorm(LSCxb,u,s/(n^0.5))
  pb=1-pa-pc
  
  pas=1-pchisq(LSCqui*(s0^2/s^2),(n-1))
  pcs=pchisq(LICqui*(s0^2/s^2),(n-1))
  pbs=1-pas-pcs
  ARL1<-1/(1-pb*pbs)
  P1=(u^2+abs(s^2-1))*ARL1
  
  return(P1)
}

resultado1 <- adaptIntegrate(Inte, lowerLimit = lowerLimit, upperLimit = upperLimit)
EQLA=resultado1$integral*(R1*R2)

resultado2 <- adaptIntegrate(Inte_double, lowerLimit = lowerLimit, upperLimit = upperLimit)
EQLB=resultado2$integral*(R1*R2)

# Function that calculates the EQL EWMA. Change as needed.
Inte_EWMA <- function(params2) {
  u=params2[1]
  ss=params2[2]
  Result<-matrix(0, runs, 1)
  
  for(i in 1:runs){
    ua=0 #starting value
    sa=c #starting value
    s=0
    D<-c()
    s<-0
    while(s<1){
      RR<-rnorm(n,u,ss)
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
  P3=(u^2+abs(ss^2-1))*ARL
  
  return(P3)
}


RR1=runif(corridas,lowerLimit[1],upperLimit[1])
RR2=runif(corridas,lowerLimit[2],upperLimit[2])
Saida=matrix(0,corridas,3)
for (i in 1:corridas){
 a1=RR1[i]
 a2=RR2[i]
 v=c(a1,a2)
 P1=Inte(v)
 P2=Inte_double(v)
 P3=Inte_EWMA(v)
  Saida[i,1]=P1
  Saida[i,2]=P2
  Saida[i,3]=P3
}

EQLA_MC=mean(Saida[,1])
EQLB_MC=(mean(Saida[,2]))
EQLC_MC=(mean(Saida[,3]))

cat("Numerical Integration","\n")
cat("EQL-Xb-S2=",EQLA,"\n")
cat("EQL-Xb-S2 Klein=",EQLB,"\n")

cat("Monte Carlo Simulation","\n")
cat("EQL-Xb-S2 MC=",EQLA_MC,"\n")
cat("EQL-Xb-S2 Klein MC=",EQLB_MC,"\n") 
cat("EQL-EWMA Xb-S2 MC=",EQLC_MC,"\n") 
toc()