###### load RCPP based functions ####
source("semicompAFT.R")

n<-1000   #number of observations
sigma<-0.5 #true sigma
c01<-2  #true h0_01
c02<-3  #true h0_02
c12<-2  #true h0_12
m<-15  #parameter that dictates censoring percentage

## probabilities for Bernoulli variables
ber_p1<-0.5
ber_p2<-0.5
ber_p3<-0.25
ber_p4<-0.6

## True coefficients
beta01<-c(1,1/2,0,0)
beta02<-c(0,1,1,0)
beta12<-c(1/2,1/2,0,1)
true_beta01<-beta01[beta01!=0]
true_beta02<-beta02[beta02!=0]
true_beta12<-beta12[beta12!=0]

## Variables indexes
var_index01<-c(1,2)
var_index02<-c(2,3)
var_index12<-c(1,2,4)

set.seed(3)

#data sampling function
#change gamma with of without frailty
sample_data<-function(n,ber_p1,ber_p2,ber_p3,ber_p4,sigma,beta01,beta02,beta12,c01,c02,c12,m){
  ###Bernoulli variables
  X1 <- runif(n,-1,1)
  X2 <- rbinom(n,1,ber_p2)
  X3 <- runif(n,-1,1)
  X4 <- runif(n,-1,1)
  X  <- cbind(X1,X2,X3,X4)
  colnames(X) <- c("X1","X2","X3","X4")
  
  ###sample frailty values
  
  #gamma <- rep(1,n) #without frailty
  gamma <-rgamma(n,shape=1/sigma,scale=sigma) #with frailty
  
  ###compute the failure times
  unif01 <- runif(n,min=0,max=1)
  unif02 <- runif(n,min=0,max=1)
  
  ###compute T1 and T2
  T1 <- sqrt(-2*log(unif01)*exp(2*X%*%beta01)/(c01*gamma))
  colnames(T1)<-c("T1")
  T2 <- sqrt(-2*log(unif02)*exp(2*X%*%beta02)/(c02*gamma))
  
  ##if T1<T2 then compute T12
  unif12       <- runif(n,min=0,max=1)
  T2_0         <- sqrt((-2*log(unif12)/(c12*gamma))*exp(2*X%*%beta12)+T1^2)
  T2[T2>T1]    <- T2_0[T2>T1]
  colnames(T2) <- c("T2")
  
  ###sample censoring values
  Cens <- runif(n=n,min=0,max=m)
  
  ###compute the ovserved data
  V      <- pmin(T1,T2,Cens); colnames(V) <- c("V")
  delta1 <- ifelse(T1<=pmin(T2,Cens),1,0); colnames(delta1) <- c("delta1")
  W      <- delta1*pmin(T2,Cens); colnames(W) <- c("W")
  delta2 <- ifelse(T2<=pmin(T1,Cens),1,0); colnames(delta2)<-c("delta2")
  delta3 <- delta1*(T2<=Cens); colnames(delta3)<-c("delta3")
  
  ###data constraction
  data <- cbind(X,V,W,delta1,delta2,delta3,true_gamma=gamma)
  #data <-as.data.frame(data)
  return(data)
}

data<-sample_data(n=n,ber_p1=ber_p1,ber_p2=ber_p2,ber_p3=ber_p3,ber_p4=ber_p4,sigma=sigma,beta01=beta01,beta02=beta02,beta12=beta12,c01=c01,c02=c02,c12=c12,m=m)

# examples for estimation with or without frailty

var_index01<-c(1,2)
var_index02<-c(2,3)
var_index12<-c(1,2,4)
X01<-as.matrix(data[,var_index01])
X02<-as.matrix(data[,var_index02])
X12<-as.matrix(data[,var_index12])
delta1<-data[,"delta1"]
delta2<-data[,"delta2"]
delta3<-data[,"delta3"]
V<-data[,"V"]
W<-delta1*data[,"W"]


# with frailty
results_with_frailty<-estimation_with_frailty(X01=X01,X02=X02,X12=X12,V=V,W=W,delta1=delta1,delta2=delta2,delta3=delta3,B=0,print=T)

# without frailty
results_without_frailty<-estimation_without_frailty(X01=X01,X02=X02,X12=X12,V=V,W=W,delta1=delta1,delta2=delta2,delta3=delta3,B=0,print=T)
