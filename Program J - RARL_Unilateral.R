library(cubature)
library(pracma)
tic()
clear ()
#In-control 
u0 = 0
s0 = 1
n =7
corridas=1000000 #EQL by Monte Carlo simulation
#Klein Limits. Change as needed.
LCLxba = -0.7422006
UCLxba = 0.7422006
LCs2a=13.22683 
#Xb-S2 Limits
LSCxb=1.209286
LICxb=-1.209286
LSCqui=21.7843

#Limits for Double Integrals. Change as needed.
lowerLimit <- c(0, 1)  # Lower bound for u and s
upperLimit <- c(1.5, 1.5)  # Limite superior para u e s
#Uniform Distribution
R1=(1/(upperLimit[1]-lowerLimit[1]))
R2=(1/(upperLimit[2]-lowerLimit[2]))

# Function that calculates the ARL Klein.Change as needed. 
Inte_double <- function(params) {
  u=params[1]
  s=params[2]
  
  # Probability Calculation
  pxi = pnorm(LCLxba, u, s / sqrt(n))
  pxs = 1 - pnorm(UCLxba, u, s / sqrt(n))
  pxc = 1 - pxs - pxi
  
  psf = 1 - pchisq(LCs2a * (s0^2 / s^2), (n - 1))
  psc = 1 - psf
  
  MarkovChain <- matrix(0, 15, 15)
  c1 <- c(pxc * psc, pxs * psc, 0, pxi * psc, 0, pxc * psf, pxs * psf, 0, pxi * psf, 0, 0, 0, 0, 0, 0) 
  c2 <- c(pxc * psc, 0, pxs * psc, pxi * psc, 0, pxc * psf, 0, pxs * psf, pxi * psf, 0, 0, 0, 0, 0, 0) 
  c3 <- c(pxc * psc, pxs * psc, 0, pxi * psc, 0, pxc * psf, pxs * psf, 0, pxi * psf, 0, 0, 0, 0, 0, 0) 
  c4 <- c(pxc * psc, pxs * psc, 0, 0, pxi * psc, pxc * psf, pxs * psf, 0, 0, pxi * psf, 0, 0, 0, 0, 0) 
  c5 <- c(pxc * psc, pxs * psc, 0, pxi * psc, 0, pxc * psf, pxs * psf, 0, pxi * psf, 0, 0, 0, 0, 0, 0) 
  c6 <- c(pxc * psc, pxs * psc, 0, pxi * psc, 0, 0, 0, 0, 0, 0, pxc * psf, pxs * psf, 0, pxi * psf, 0) 
  c7 <- c(pxc * psc, 0, pxs * psc, pxi * psc, 0, 0, 0, 0, 0, 0, pxc * psf, 0, pxs * psf, pxi * psf, 0)   
  c8 <- c(pxc * psc, pxs * psc, 0, pxi * psc, 0, pxc * psf, pxs * psf, 0, pxi * psf, 0, 0, 0, 0, 0, 0)   
  c9 <- c(pxc * psc, pxs * psc, 0, 0, pxi * psc, 0, 0, 0, 0, 0, pxc * psf, pxs * psf, 0, 0, pxi * psf)   
  c10 <- c(pxc * psc, pxs * psc, 0, pxi * psc, 0, pxc * psf, pxs * psf, 0, pxi * psf, 0, 0, 0, 0, 0, 0)   
  c11 <- c(pxc * psc, pxs * psc, 0, pxi * psc, 0, pxc * psf, pxs * psf, 0, pxi * psf, 0, 0, 0, 0, 0, 0) 
  c12 <- c(pxc * psc, pxs * psc, 0, pxi * psc, 0, pxc * psf, pxs * psf, 0, pxi * psf, 0, 0, 0, 0, 0, 0) 
  c13 <- c(pxc * psc, pxs * psc, 0, pxi * psc, 0, pxc * psf, pxs * psf, 0, pxi * psf, 0, 0, 0, 0, 0, 0) 
  c14 <- c(pxc * psc, pxs * psc, 0, pxi * psc, 0, pxc * psf, pxs * psf, 0, pxi * psf, 0, 0, 0, 0, 0, 0) 
  c15 <- c(pxc * psc, pxs * psc, 0, pxi * psc, 0, pxc * psf, pxs * psf, 0, pxi * psf, 0, 0, 0, 0, 0, 0) 
  MarkovChain <- rbind(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15)
  
  A <- t(MarkovChain) - diag(15)
  A[15,] <- rep(1, 15)
  B <- rep(0, 15)
  B[15] <- 1
  
  Solved_markov_chain <- solve(A) %*% B
  ARL1 <- 1 / sum(Solved_markov_chain[c(3, 5, 8, 10, 11, 12, 13, 14, 15)])
  ARL1 <- ARL1  
  P2=ARL1
  return(P2)
}
# Function that calculates the ARL Xb-S2. Change as needed.
Inte <- function(params1){
  u=params1[1]
  s=params1[2]
  pc=pnorm(LICxb,u,s/(n^0.5))
  pa=1-pnorm(LSCxb,u,s/(n^0.5))
  pb=1-pa-pc
  
  pas=1-pchisq(LSCqui*(s0^2/s^2),(n-1))
  pbs=1-pas
  ARL1<-1/(1-pb*pbs)
  P1=ARL1
  return(P1)
}
#Function that calculates the RARL.Change as needed.
Raz <- function(params){
u=params[1]
s=params[2]
R1=Inte_double(params)/Inte(params)
return(R1)
}
resultado1 <- adaptIntegrate(Raz, lowerLimit = lowerLimit, upperLimit 
                             = upperLimit)
RARL=(resultado1$integral*(R1*R2))


#Simulation Monte Carlo

RR1=runif(corridas,lowerLimit[1],upperLimit[1])
RR2=runif(corridas,lowerLimit[2],upperLimit[2])
Saida=matrix(0,corridas,2)
for (i in 1:corridas){
 a1=RR1[i]
 a2=RR2[i]
 params=c(a1,a2)
 params1=c(a1,a2)
 P1=Inte(params1)
 P2=Inte_double(params)
  Saida[i,1]=P1/P1 #Change as needed.
  Saida[i,2]=P2/P1 #Change as needed.
  
  
}

#Change as needed.
RARLA_MC=mean(Saida[,1])
RARLB_MC=(mean(Saida[,2]))
cat("Integration Numeric","\n")
cat("RARL-Xb-S2 Klein=",RARL,"\n")

cat("Monte Carlo Simulation","\n")
cat("RARLA_MC=",RARLA_MC,"\n")
cat("RARLB_MC Klein=",RARLB_MC,"\n") 

toc()