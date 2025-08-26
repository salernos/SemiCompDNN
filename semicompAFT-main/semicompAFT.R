
###### Load required packages ######
suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(RcppArmadillo))
suppressPackageStartupMessages(library(aftgee))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(inline))
suppressPackageStartupMessages(library(lbfgs))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(tidyr))

###### load RCPP based functions ####
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
source(here("R", "semicompAFT-main", "rcpp_based_functions.R"))

###### R Functions ###########

# function for computing bandwidths
bandwidths<-function(V,W,delta1,delta2,delta3,ah01,ah02,ah12){
  n<-length(V)
  delta1_0 <- delta1
  delta2_0 <- delta2
  W        <- W[delta1_0==1]
  delta3_0 <- delta3[delta1_0==1]
  sd1_ev <- sd(log(V[delta1_0==1]))
  sd1    <- sd(log(V))
  sd2_ev <- sd(log(V[delta2_0==1]))
  sd2    <- sd(log(V))
  sd12_ev<- sd(log(W[delta3_0==1]))
  sd12   <- sd(log(W))
  #Find the bandwidths
  a_n1   <- ah01*c((8*sqrt(2)/3)^(1/5)*sd1_ev *(n)^(-1/5),4^(1/3)*sd1*n^(-1/3))
  a_n2   <- ah02*c((8*sqrt(2)/3)^(1/5)*sd2_ev *(n)^(-1/5),4^(1/3)*sd2*n^(-1/3))
  a_n12  <- ah12*c((8*sqrt(2)/3)^(1/5)*sd12_ev*(sum(delta1==1))^(-1/5),4^(1/3)*sd12*sum(delta1==1)^(-1/3))
  return(rbind(a_n1=a_n1,a_n2=a_n2,a_n12=a_n12))
}

# function for estimating h0_01 or h0_02
h0_hat_01_02_f<-function(t,beta,gamma,a_n,X,V,delta){
  n   <- dim(X)[1]
  R   <- log(V)-X%*%beta
  RRn <- matrix(rep(R,length(t))-rep(log(t),each=n),nrow=n)

  kernel_density  <- dnorm(RRn[delta==1,,drop=F]/a_n[1])
  kernel_integral <- pnorm(RRn/a_n[2])
  s3 <- colSums(kernel_density)/(a_n[1]*t)
  s4 <- colSums(gamma*kernel_integral)
  h0 <- s3/s4
  h  <-ifelse(s4==0,0,h0)
  return(h)
}

# function for estimating bootstrap h0_01 or h0_02
h0_hat_01_02_perturbed_f<-function(t,beta,gamma,a_n,X,V,delta,G){
  n   <- dim(X)[1]
  R   <- log(V)-X%*%beta
  RRn <- matrix(rep(R,length(t))-rep(log(t),each=n),nrow=n)

  kernel_density  <- dnorm(RRn[delta==1,,drop=F]/a_n[1])
  kernel_integral <- pnorm(RRn/a_n[2])
  s3 <- (t(G[delta==1])%*%kernel_density)/(a_n[1]*t)
  s4 <- colSums(G*gamma*kernel_integral)

  h <- s3/c(s4)
  return(h)
}

# compute integral of h0_01 or h0_02 from x to y
H0_hat_01_02_f<-function(x,y,beta,gamma,a_n,X,V,delta){
  int  <- integrate(h0_hat_01_02_f,gamma=gamma,beta = beta,a_n=a_n,X=X,V=V,delta=delta,lower=x,upper=y,subdivisions = 1000)$val
   return(int)
}

# compute integral of bootstrap h0_01 or h0_02 from x to y
H0_hat_01_02_perturbed_f<-function(x,y,beta,gamma,a_n,X,V,delta,G){
  int <-  integrate(h0_hat_01_02_perturbed_f,gamma=gamma,beta = beta,a_n=a_n,X=X,V=V,delta=delta,G=G,lower=x,upper=y,subdivisions = 10000,rel.tol = 1e-8)$val
  return(int)
}

# compute H0_01 or H0_02 at observed times
H0_hat_01_02_observed_f<-function(beta,gamma,a_n,X,V,delta){
  n    <- dim(X)[1]
  e_R1 <- V*exp(-X%*%beta)
  e_R1 <- cbind(e_R1,order=1:n)
  colnames(e_R1) <- c("e_R1","order")
  e_R1_ord       <- e_R1[order(e_R1[,1]),]
  e_R1_ord2col   <- cbind(c(0,e_R1_ord[-n,1]),e_R1_ord)

  cl <- parallel::makeCluster(4)
  doParallel::registerDoParallel(cl) # register the cluster
  v1 <- foreach(i = 1:n
                ,.combine = "rbind"
                ,.export = c("h0_hat_01_02_f")
  ) %dopar% {
    int <- integrate(h0_hat_01_02_f,gamma=gamma,beta = beta,a_n=a_n,X=X,V=V,delta=delta,lower=e_R1_ord2col[i,1],upper=e_R1_ord2col[i,2],subdivisions = 10000)$val
    return(int)
  }
  parallel::stopCluster(cl)

  if (sum(is.na(v1))>0){
    v2            <- v1
    v2[is.na(v1)] <- 0
    ind_na        <- which(is.na(v1))
    len_na        <-length(ind_na)
    v2            <-cumsum(v2)
    for (k in 1:len_na){
      m       <- ind_na[k]
      int     <- integrate(h0_hat_01_02_f,gamma=gamma,beta = beta,a_n=a_n,X=X,V=V,delta=delta,lower=0,upper=e_R1_ord2col[m,2],subdivisions = 1000)$val-sum(v2[1:(m-1)])
      v2[m:n] <- v2[m:n]+int-v2[m]
    }
    v2           <- cbind(v2,e_R1_ord2col[,2],e_R1_ord2col[,3])
    colnames(v2) <- c("integral_ord","e_R1_ord","order")
    H01_m        <- merge(e_R1,v2,by="order",sort=F)[,"integral_ord"]

  } else {
    v1           <- cbind(cumsum(v1),e_R1_ord2col[,2],e_R1_ord2col[,3])
    colnames(v1) <- c("integral_ord","e_R1_ord","order")
    H01_m        <- merge(e_R1,v1,by="order",sort=F)[,"integral_ord"]
  }
  return(H01_m)
}

# compute perturbed H0_01 or H0_02 at observed times
H0_hat_01_02_observed_perturbed_f<-function(beta,gamma,a_n,X,V,delta,G){
  n    <- dim(X)[1]
  e_R1 <- V*exp(-X%*%beta)
  e_R1 <- cbind(e_R1,order=1:n)
  colnames(e_R1) <- c("e_R1","order")
  e_R1_ord       <- e_R1[order(e_R1[,1]),]
  e_R1_ord2col   <- cbind(c(0,e_R1_ord[-n,1]),e_R1_ord)

  cl <- parallel::makeCluster(4)
  doParallel::registerDoParallel(cl) # register the cluster
  v1 <- foreach(i = 1:n
                ,.combine = "rbind"
                ,.export = c("h0_hat_01_02_perturbed_f")
  ) %dopar% {
    t <-  try({
      int <-   integrate(h0_hat_01_02_perturbed_f,gamma=gamma,beta = beta,a_n=a_n,X=X,V=V,delta=delta,G=G,lower=e_R1_ord2col[i,1],upper=e_R1_ord2col[i,2],subdivisions = 1000)$val
    })
    if (inherits(t, "try-error")) {
      int <-  distrEx::GLIntegrate(h0_hat_01_02_perturbed_f,gamma=gamma,beta = beta,a_n=a_n,X=X,V=V,delta=delta,G=G,lower=e_R1_ord2col[i,1],upper=e_R1_ord2col[i,2])
    }
    return(int)
  }
  parallel::stopCluster(cl)

  v1           <- cbind(cumsum(v1),e_R1_ord2col[,2],e_R1_ord2col[,3])
  colnames(v1) <- c("integral_ord","e_R1_ord","order")
  H01_m        <- merge(e_R1,v1,by="order",sort=F)[,"integral_ord"]
  return(H01_m)
}

# function for estimating h0_12
h0_hat_12_f<-function(t,beta,gamma,a_n12,X12,V,W,delta1,delta3){
  X0       <- X12[delta1==1,]
  W0       <- W[delta1==1]
  V0       <- V[delta1==1]
  delta3   <- delta3[delta1==1]
  nR       <- dim(X0)[1]
  lt       <- length(t)
  x_beta   <- X0%*%beta
  R        <- log(W0)-x_beta
  mat_logt <- rep(log(t),each=nR)
  RRn      <- matrix(rep(R,lt)-mat_logt,nrow=nR)
  R2v      <- log(V0)-x_beta
  RRn2v    <- matrix(rep(R2v,lt)-mat_logt,nrow=nR)

  kernel_density    <- dnorm(RRn[delta3==1,,drop=F]/a_n12[1])
  kernel_integral_W <- pnorm(RRn/a_n12[2])
  kernel_integral_V <- pnorm(RRn2v/a_n12[2])
  gamma_sel         <- gamma[delta1==1]

  s3 <- colSums(kernel_density)/(a_n12[1]*t)
  s4 <- colSums(gamma_sel*(kernel_integral_W-kernel_integral_V))
  h0 <- s3/s4
  h  <-ifelse(s4==0,0,h0)
  return(h)
}

# function for estimating perturbed h0_12
h0_hat_12_perturbed_f<-function(t,beta,gamma,a_n12,X12,V,W,delta1,delta3,G){
  X0       <- X12[delta1==1,]
  W0       <- W[delta1==1]
  V0       <- V[delta1==1]
  delta3   <- delta3[delta1==1]
  nR       <- dim(X0)[1]
  lt       <- length(t)
  x_beta   <- X0%*%beta
  R        <- log(W0)-x_beta
  mat_logt <- rep(log(t),each=nR)
  RRn      <- matrix(rep(R,lt)-mat_logt,nrow=nR)
  R2v      <- log(V0)-x_beta
  RRn2v    <- matrix(rep(R2v,lt)-mat_logt,nrow=nR)

  kernel_density    <- dnorm(RRn[delta3==1,,drop=F]/a_n12[1])
  kernel_integral_W <- pnorm(RRn/a_n12[2])
  kernel_integral_V <- pnorm(RRn2v/a_n12[2])
  gamma_sel         <- gamma[delta1==1]

  s3 <- (t(G[delta1==1][delta3==1])%*%kernel_density)/(a_n12[1]*t)
  s4 <- colSums((G[delta1==1])*gamma_sel*(kernel_integral_W-kernel_integral_V))
  h0 <- s3/s4
  h  <- ifelse(s4==0,0,h0)
  return(h)
}

#compute integral of h0_12 from x to y
H0_hat_12_f<-function(x,y,beta,gamma,a_n12,X12,V,W,delta1,delta3){
  if (x==0) {int<-0} else {
    if (y==0) {int<-0} else
            int<-integrate(h0_hat_12_f,gamma=gamma,beta= beta,a_n12 = a_n12,X12=X12,V=V,W=W,delta1=delta1,delta3=delta3,lower=x,upper = y,subdivisions = 1000)$val
   }
  return(int)}

#compute integral of perturbed h0_12 from x to y
H0_hat_12_perturbed_f<-function(x,y,beta,gamma,a_n12,X12,V,W,delta1,delta3,G){
  if (x==0) {int<-0} else {
    if (y==0) {int<-0} else
      int<-integrate(h0_hat_12_perturbed_f,gamma=gamma,beta= beta,a_n12 = a_n12,X12=X12,V=V,W=W,delta1=delta1,delta3=delta3,G=G,lower=x,upper = y,subdivisions = 10000,rel.tol = 1e-8)$val
  }
  return(int)}

# compute H0_12 at observed times
H0_hat_12_observed_f_new<-function(beta,gamma,a_n12,X12,V,W,delta1,delta3){
  n          <- length(V)
  exp_x_beta <- exp(-X12%*%beta)
  e_Rw12     <- W*exp_x_beta
  e_Rw12[delta1==0] <- 0
  e_Rv12     <- V*exp_x_beta
  e_Rv12[delta1==0] <- 0

  e_Rw12_new <- data.frame(times=e_Rw12,type=rep("W",n),order=1:n)
  e_Rv12_new <- data.frame(times=e_Rv12,type=rep("V",n),order=1:n)
  times      <- rbind(e_Rv12_new,e_Rw12_new)
  times_ord  <- times[order(times[,1]),]

  times_ord2col <- cbind(c(0,times_ord[-2*n,1]),times_ord,1:(2*n))

  cl  <-  parallel::makeCluster(4)
  doParallel::registerDoParallel(cl) # register the cluster
  integral_interval <- foreach(i = 1:(2*n)
                               ,.combine = "rbind"
                               ,.export = c("h0_hat_12_f","H0_hat_12_f")
  ) %dopar% {
    int<-
      tryCatch({
        suppressWarnings (
          H0_hat_12_f(times_ord2col[i,1],times_ord2col[i,2],gamma=gamma,beta = beta,a_n12=a_n12,X12=X12,W=W,V=V,delta1=delta1,delta3=delta3)
        )},error=function(e) NA)
    return(int)
  }
  parallel::stopCluster(cl)

  integral_cumsum      <- data.frame(integral=cumsum(integral_interval),type=times_ord2col[,"type"],order=times_ord2col[,"order"],times=times_ord2col[,2])
  integral_cumsum_wide <- spread(integral_cumsum[,-4],type,integral)

  H12 <- integral_cumsum_wide$W-integral_cumsum_wide$V

  ind_na<-which(is.na(H12))
  len_na<-length(ind_na)
  if (len_na>0){
    for (k in 1:len_na){
      m<-ind_na[k]
      H12[m]<-integrate(h0_hat_12_f,gamma=gamma,beta= beta,a_n12 = a_n12,X12=X12,V=V,W=W,delta1=delta1,delta3=delta3,lower=(V*exp(-X12%*%beta))[m],upper = (W*exp(-X12%*%beta))[m],subdivisions = 10000)$val
    }
  }
  return(H12)
}

# compute perturbed H0_12 at observed times
H0_hat_12_observed_perturbed_f_new<-function(beta,gamma,a_n12,X12,V,W,delta1,delta3,G){
  n<-length(V)
  exp_x_beta <- exp(-X12%*%beta)
  e_Rw12     <- W*exp_x_beta
  e_Rw12[delta1==0] <- 0
  e_Rv12     <- V*exp_x_beta
  e_Rv12[delta1==0] <- 0

  e_Rw12_new <- data.frame(times=e_Rw12,type=rep("W",n),order=1:n)
  e_Rv12_new <- data.frame(times=e_Rv12,type=rep("V",n),order=1:n)
  times      <- rbind(e_Rv12_new,e_Rw12_new)
  times_ord  <- times[order(times[,1]),]

  times_ord2col <- cbind(c(0,times_ord[-2*n,1]),times_ord,1:(2*n))

  cl  <-  parallel::makeCluster(4)
  doParallel::registerDoParallel(cl) # register the cluster
  integral_interval <- foreach(i = 1:(2*n)
                               ,.combine = "rbind"
                               ,.export = c("h0_hat_12_perturbed_f")
  ) %dopar% {
    if (times_ord2col[i,1]==0) {int<-0} else {
      if (times_ord2col[i,2]==0) {int<-0} else {
        int<-tryCatch({
          suppressWarnings(integrate(h0_hat_12_perturbed_f,G=G,gamma=gamma,beta= beta,a_n12 = a_n12,X12=X12,V=V,W=W,delta1=delta1,delta3=delta3,lower=times_ord2col[i,1],upper = times_ord2col[i,2],subdivisions = 10000)$val)
        },error=function(e) NA)
      }}
    return(int)
  }
  parallel::stopCluster(cl)

  integral_cumsum      <- data.frame(integral=cumsum(integral_interval),type=times_ord2col[,"type"],order=times_ord2col[,"order"],times=times_ord2col[,2])
  integral_cumsum_wide <- spread(integral_cumsum[,-4],type,integral)

  H12 <- integral_cumsum_wide$W-integral_cumsum_wide$V

  ind_na<-which(is.na(H12))
  len_na<-length(ind_na)
  if (len_na>0){
    for (k in 1:len_na){
      m<-ind_na[k]
      H12[m]<-integrate(h0_hat_12_perturbed_f,G=G,gamma=gamma,beta= beta,a_n12 = a_n12,X12=X12,V=V,W=W,delta1=delta1,delta3=delta3,lower=(V*exp(-X12%*%beta))[m],upper = (W*exp(-X12%*%beta))[m],subdivisions = 10000)$val
    }
  }
  return(H12)}

# compute posterior expectations of gamma and log(gamma)
posterior_expectations<-function(sigma,delta1,delta2,delta3,H0_01_obs,H0_02_obs,H0_12_obs){
  g     <- delta1+delta2+delta3
  k     <- H0_01_obs+H0_02_obs+delta1*H0_12_obs
  a     <- 1/sigma+g
  b     <- 1/sigma+k
  E1    <- a/b
  E2    <- digamma(a)-log(b)
  sumE1 <- sum(E1)
  sumE2 <- sum(E2)
  return(list(E1=E1,E2=E2,sumE1=sumE1,sumE2=sumE2))
}

# estimating sigma
sigma_opt<-function(init_sigma,sumE1,sumE2,n){
  optim(par=init_sigma,fn=function(x) {-( (1/x)*log(1/x)+(1/x-1)*sumE2/n-(1/x)*sumE1/n-log(gamma(1/x)) )}, method="L-BFGS-B",lower=0.01)$par
}

# data sampling function
sample_data<-function(n,ber_p1,ber_p2,ber_p3,ber_p4,sigma,beta01,beta02,beta12,c01,c02,c12,m){

  ###Bernoulli variables
  X1 <- runif(n,-1,1)
  X2 <- rbinom(n,1,ber_p2)
  X3 <- runif(n,-1,1)
  X4 <- runif(n,-1,1)
  X  <- cbind(X1,X2,X3,X4)
  colnames(X) <- c("X1","X2","X3","X4")

  ###sample frailty values
  gamma <- rgamma(n,shape=1/sigma,scale=sigma)

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
  delta2 <- ifelse(T2<=pmin(T1,Cens),1,0); colnames(delta2) <- c("delta2")
  delta3 <- delta1*(T2<=Cens); colnames(delta3) <- c("delta3")

  ###data constraction
  data <- cbind(X,V,W,delta1,delta2,delta3,true_gamma=gamma)
  #data <-as.data.frame(data)
  return(data)
}

# statistics for estimated parameters
stat_f<-function(x,true_val){
  mean <- rowMeans2(x)
  esd  <- rowSds(x)
  cr   <- rowMeans2(x-1.96*esd<=true_val&true_val<=x+1.96*esd)
  stat <- rbind(true_val,mean,esd,cr)
  return(stat)
}

########## Estimation functions

# function for estimation with frailty
estimation_with_frailty<-function(X01,X02,X12,V,W,delta1,delta2,delta3,zeta_beta=50,zeta_h=0.01,initial_sigma=2,stop_iter_num=1000,conv_betas_bound=0.00001,conv_Hs_bound=0.0001,conv_sigma_bound=0.0001,B,print)
{
  n_pars<-1+dim(X01)[2]+dim(X02)[2]+dim(X12)[2]
  n <- length(V)
  D <- delta1+delta2+delta3
  initial_gamma <- rep(1,n)

  ## bandwidths computation
  ab01  <-ab02 <- ab12 <- zeta_beta
  ah01  <-ah02 <- ah12 <- zeta_h
  band  <- bandwidths(V=V,W=W,delta1=delta1,delta2=delta2,delta3=delta3,ah01=ah01,ah02=ah02,ah12=ah12)
  a_n1  <- band[1,]
  a_n2  <- band[2,]
  a_n12 <- band[3,]

  #estimation

  #initiate first iteration
  m=1
  gamma_m <- initial_gamma
  sigma_m <- initial_sigma

  ####set environments for cpp functions####
  #environment for optimization of l01
  env01             <- new.env()
  env01[["X"]]      <- X01
  env01[["X0"]]     <- X01[delta1==1,]
  env01[["delta"]]  <- delta1
  env01[["V"]]      <- V
  env01[["V0"]]     <- V[delta1==1]
  env01[["a"]]      <- a_n1*ab01
  #environment for optimization of l02
  env02             <- new.env()
  env02[["X"]]      <- X02
  env02[["X0"]]     <- X02[delta2==1,]
  env02[["delta"]]  <- delta2
  env02[["V"]]      <- V
  env02[["V0"]]     <- V[delta2==1]
  env02[["a"]]      <- a_n2*ab02
  #environment for optimization of l12
  env12             <- new.env()
  env12[["X"]]      <- X12[delta1==1,]
  env12[["X0"]]     <- X12[delta3==1,]
  env12[["delta3"]] <- delta3[delta1==1]
  env12[["V"]]      <- V[delta1==1]
  env12[["W"]]      <- W[delta1==1]
  env12[["W0"]]     <- W[delta3==1]
  env12[["a"]]      <- a_n12*ab12
  env12[["n"]]      <- n
  ######
  a <- Sys.time()
  # set overall environments values for cpp functions
  env01[["gamma_m"]] <- env02[["gamma_m"]] <- gamma_m
  env12[["gamma_m"]] <- gamma_m[delta1==1]

  # estimating Betas
  beta_hat_01_m <- aftsrr(Surv(V, delta1) ~ X01, se = "ISMB",B=0)$beta
  beta_hat_02_m <- aftsrr(Surv(V, delta2) ~ X02, se = "ISMB",B=0)$beta
  beta_hat_12_m <- aftsrr(Surv(W[delta1==1], delta3[delta1==1]) ~ X12[delta1==1,], se = "ISMB",B=0)$beta


  # estimating H0s
  H0_01_obs_m  <- H0_hat_01_02_observed_f(beta=beta_hat_01_m,gamma=gamma_m,a_n=a_n1,X=X01,V=V,delta=delta1)
  H0_02_obs_m  <- H0_hat_01_02_observed_f(beta=beta_hat_02_m,gamma=gamma_m,a_n=a_n2,X=X02,V=V,delta=delta2)
  H0_12_obs_m  <- H0_hat_12_observed_f_new(beta=beta_hat_12_m,gamma=c(gamma_m),a_n12=a_n12,X12=X12,V=V,W=W,delta1 = delta1,delta3 = delta3)

  # computing E1 and E2
  posterior_expectations_m_plus1 <- posterior_expectations(sigma=sigma_m,H0_01_obs_m,H0_02_obs_m,H0_12_obs_m,delta1=delta1,delta2=delta2,delta3=delta3)
  E1_m_plus1                     <- posterior_expectations_m_plus1[["E1"]]
  env01[["gamma_m"]]             <- env02[["gamma_m"]]<- E1_m_plus1
  env12[["gamma_m"]]             <- E1_m_plus1[delta1==1]

  # estimate sigma
  sigma_m_plus1 <- initial_sigma

  #print results of the first iteration
  if(print==T){
  print("Estimation proccess begun")
  print(paste("iteration m=",m))
  print(paste(paste0(colnames(X01),"_01=",round(beta_hat_01_m,6))))
  print(paste(paste0(colnames(X02),"_02=",round(beta_hat_02_m,6))))
  print(paste(paste0(colnames(X12),"_12=",round(beta_hat_12_m,6))))
  print(paste("s=",round(sigma_m_plus1,6)))
  }

  ##### Next iterations ####

  # define convergence indicators
  conv_beta01 <- conv_beta02 <- conv_beta12 <- conv_H0_01 <- conv_H0_02 <- conv_H0_12 <- conv_sigma <-conv_betas<-conv_Hs<- conv<-F

  while(conv==F){
    # define stopping iteration
    if(m == stop_iter_num) {
      conv <- F
      break
    }

    init <- Sys.time()
    m<-m+1

    ## estimating Betas
    if(conv_beta01==F){
    l01 <- lbfgs(l_s_m_01_02_f.CPP(),grad_beta01_02_f.CPP(),vars=beta_hat_01_m, environment =env01,invisible = 1)
    beta_hat_01_m_plus1 <- l01$par
    } else {beta_hat_01_m_plus1 <- beta_hat_01_m}

    if(conv_beta02==F){
    l02 <- lbfgs(l_s_m_01_02_f.CPP(),grad_beta01_02_f.CPP(),vars=beta_hat_02_m, environment =env02,invisible = 1)
    beta_hat_02_m_plus1<-l02$par
    } else {beta_hat_02_m_plus1<-beta_hat_02_m}

    if (conv_beta01==T&conv_beta02==T){
    l12 <- lbfgs(l_s_m_12_f.CPP(),grad_beta12_f.CPP(),vars=beta_hat_12_m, environment =env12,invisible = 1,linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING")
    beta_hat_12_m_plus1 <- l12$par
    } else {beta_hat_12_m_plus1 <-beta_hat_12_m}

    ##estimating H0s
    start<-Sys.time()
    if ((conv_beta01&conv_Hs)==F){
      H0_01_obs_m_plus1 <-H0_hat_01_02_observed_f(beta=beta_hat_01_m_plus1,gamma=c(E1_m_plus1),a_n=a_n1,X=X01,V=V,delta=delta1)
    } else {H0_01_obs_m_plus1<-H0_01_obs_m}
    end<-Sys.time()
    #print(end-start)

    start<-Sys.time()
    if ((conv_beta02&conv_Hs)==F) {
      H0_02_obs_m_plus1<- H0_hat_01_02_observed_f(beta=beta_hat_02_m_plus1,gamma=c(E1_m_plus1),a_n=a_n2,X=X02,V=V,delta=delta2)
    } else {H0_02_obs_m_plus1<-H0_02_obs_m}
    end<-Sys.time()
    #print(end-start)

    start<-Sys.time()
    if ((conv_beta12&conv_Hs)==F) {
      H0_12_obs_m_plus1<-H0_hat_12_observed_f_new(beta=beta_hat_12_m_plus1,gamma=c(E1_m_plus1),a_n12=a_n12,X12=X12,V=V,W=W,delta1 = delta1,delta3 = delta3)
    } else {H0_12_obs_m_plus1<-H0_12_obs_m}
    end<-Sys.time()
    #print(end-start)

    ## computing E1 and E2
    posterior_expectations_m_plus2 <- posterior_expectations(sigma=sigma_m_plus1,H0_01_obs_m_plus1,H0_02_obs_m_plus1,H0_12_obs_m_plus1,delta1=delta1,delta2=delta2,delta3=delta3)
    E1_m_plus2                     <- posterior_expectations_m_plus2[["E1"]]
    env01[["gamma_m"]]             <- env02[["gamma_m"]] <- E1_m_plus2
    env12[["gamma_m"]]             <- E1_m_plus2[delta1==1]

    sumE1_m_plus2 <- posterior_expectations_m_plus2[["sumE1"]]
    sumE2_m_plus2 <- posterior_expectations_m_plus2[["sumE2"]]

    ## estimating sigma
    sigma_m_plus2 <- sigma_opt(init_sigma = sigma_m_plus1,sumE1 = sumE1_m_plus2,sumE2 = sumE2_m_plus2,n=n)

    ## check for convergence
    #betas' convergence
    if ((sum(abs(beta_hat_01_m_plus1-beta_hat_01_m)<conv_betas_bound))==length(beta_hat_01_m)) {conv_beta01 <- T}
    if ((sum(abs(beta_hat_02_m_plus1-beta_hat_02_m)<conv_betas_bound))==length(beta_hat_02_m)) {conv_beta02 <- T}
    if ((sum(abs(beta_hat_12_m_plus1-beta_hat_12_m)<conv_betas_bound))==length(beta_hat_12_m)) {conv_beta12 <- T}

    #H0s' convergence
    if (mean(abs(H0_01_obs_m_plus1-H0_01_obs_m))<conv_Hs_bound) {conv_H0_01 <- T}
    if (mean(abs(H0_02_obs_m_plus1-H0_02_obs_m))<conv_Hs_bound) {conv_H0_02 <- T}
    if (mean(abs(H0_12_obs_m_plus1-H0_12_obs_m))<conv_Hs_bound) {conv_H0_12 <- T}
    conv_Hs  <- conv_H0_01&conv_H0_02&conv_H0_12

    #sigmas' convergence
    if (abs(sigma_m_plus2-sigma_m_plus1)<conv_sigma_bound) {conv_sigma <- T}

    #overall convergence
    conv <- conv_sigma&conv_Hs&conv_beta01&conv_beta02&conv_beta12

    ## print the results of current iteration
    if (print==T){
    itter_time<-round(Sys.time()-init,2)
    print(paste("iteration m=",m))
    print(paste0(colnames(X01),"_01=",round(beta_hat_01_m_plus1,6)))
    print(paste0(colnames(X02),"_02=",round(beta_hat_02_m_plus1,6)))
    print(paste0(colnames(X12),"_12=",round(beta_hat_12_m_plus1,6)))
    print(paste("s=",round(sigma_m_plus2,6)))
    print(paste("conv_betas=",conv_beta01,conv_beta02,conv_beta12,"   ","conv_H0s=",conv_H0_01,conv_H0_02,conv_H0_12,"    conv_sigma=",conv_sigma,"    sigma_diff=",round(abs(sigma_m_plus2-sigma_m_plus1),6)))
    print(paste("time=",itter_time,"  sumE1=",round(sumE1_m_plus2,3),"  sumE2=",round(sumE2_m_plus2,3)))
    }

    ## save iterations' results
    beta_hat_01_m    <- beta_hat_01_m_plus1
    beta_hat_02_m    <- beta_hat_02_m_plus1
    beta_hat_12_m    <- beta_hat_12_m_plus1
    H0_01_obs_m      <- H0_01_obs_m_plus1
    H0_02_obs_m      <- H0_02_obs_m_plus1
    H0_12_obs_m      <- H0_12_obs_m_plus1
    E1_m_plus1       <- E1_m_plus2
    sigma_m_plus1old <- sigma_m_plus1
    sigma_m_plus1    <- sigma_m_plus2
  }

  #save the converged results
  sigma_conv  <- sigma_m_plus1old
  beta01_conv <- beta_hat_01_m_plus1
  beta02_conv <- beta_hat_02_m_plus1
  beta12_conv <- beta_hat_12_m_plus1
  tau_conv    <- c(sigma_conv=sigma_conv,beta01_conv=beta01_conv,beta02_conv=beta02_conv,beta12_conv=beta12_conv,
                   E1=posterior_expectations_m_plus2$E1,E2=posterior_expectations_m_plus2$E2,
                   H01_conv=H0_01_obs_m_plus1,H02_conv=H0_02_obs_m_plus1,H12_conv=H0_12_obs_m_plus1)

  #### store results for summary ####

  ###### Summary #####
  est_par           <- cbind(sigma_m_plus1,t(beta_hat_01_m),t(beta_hat_02_m),t(beta_hat_12_m),sum(delta1),sum(delta2),sum(delta3),difftime(Sys.time(),a,units = "hours"),m)
  colnames(est_par) <- c("sigma",paste0(colnames(X01),"-01"),paste0(colnames(X02),"-02"),paste0(colnames(X12),"-12"),
                         "illness","death01","death12","duration","iter_num")

  if(print==T){print(est_par)}

  # Bootstrap for variance estimation
  if(B>0){
  if(print==T){print("Bootstrap process begun")}

  pert_all<-sigma_pert<-b01_pert<-b02_pert<-b12_pert<-boot_time<-iter<-h0_01<-h0_02<-h0_12<-H0_01<-H0_02<-H0_12<-NULL
  i=1
  while (length(sigma_pert)<B){
    aa<-Sys.time()

    ####set environments for cpp functions####
    #environment for optimization of l01
    env01             <- new.env()
    env01[["X"]]      <- X01
    env01[["X0"]]     <- X01[delta1==1,]
    env01[["delta"]]  <- delta1
    env01[["V"]]      <- V
    env01[["V0"]]     <- V[delta1==1]
    env01[["a"]]      <- a_n1*ab01
    #environment for optimization of l02
    env02             <- new.env()
    env02[["X"]]      <- X02
    env02[["X0"]]     <- X02[delta2==1,]
    env02[["delta"]]  <- delta2
    env02[["V"]]      <- V
    env02[["V0"]]     <- V[delta2==1]
    env02[["a"]]      <- a_n2*ab02
    #environment for optimization of l12
    env12             <- new.env()
    env12[["X"]]      <- X12[delta1==1,]
    env12[["X0"]]     <- X12[delta3==1,]
    env12[["delta3"]] <- delta3[delta1==1]
    env12[["V"]]      <- V[delta1==1]
    env12[["W"]]      <- W[delta1==1]
    env12[["W0"]]     <- W[delta3==1]
    env12[["a"]]      <- a_n12*ab12
    env12[["n"]]      <- n

    #########
    G         <- rexp(n,1)
    env01$G   <- env02$G <- G
    env01$G0  <- G[delta1==1]
    env02$G0  <- G[delta2==1]
    env12$G   <- G[delta1==1]
    env12$G0  <- G[delta3==1]

    #set first iteration
    m= 1
    gamma_m       <- rep(1,n)
    sigma_m       <- initial_sigma
    env01$gamma_m <- env02$gamma_m <- gamma_m
    env12$gamma_m <- gamma_m[delta1==1]

    # estimating Betas
    beta_hat_01_m <-  lbfgs(l_s_m_01_02pert_f.CPP(),grad_beta01_02pert_f.CPP(),env=env01,vars=as.numeric(tau_conv[ ,grep("beta01_conv",colnames(tau_conv))]),invisible = T)$par
    beta_hat_02_m <-  lbfgs(l_s_m_01_02pert_f.CPP(),grad_beta01_02pert_f.CPP(),env=env02,vars=as.numeric(tau_conv[ ,grep("beta02_conv",colnames(tau_conv))]),invisible = T)$par
    beta_hat_12_m <-  aftsrr(Surv(W[delta1==1], delta3[delta1==1]) ~ X12[delta1==1,], se = "ISMB",B=0)$beta

    # estimating H0s
    H0_01_obs_m   <- H0_hat_01_02_observed_perturbed_f(beta=beta_hat_01_m,gamma=c(gamma_m),a_n=a_n1,X=X01,V=V,delta=delta1,G=G)
    H0_02_obs_m   <- H0_hat_01_02_observed_perturbed_f(beta=beta_hat_02_m,gamma=c(gamma_m),a_n=a_n2,X=X02,V=V,delta=delta2,G=G)
    H0_12_obs_m   <- H0_hat_12_observed_perturbed_f_new(beta=beta_hat_12_m,gamma=c(gamma_m),a_n12=a_n12,X12=X12,V=V,W=W,delta1 = delta1,delta3 = delta3,G=G)

    # computing E1 and E2
    posterior_expectations_m_plus1 <- posterior_expectations(sigma=sigma_m,H0_01_obs_m,H0_02_obs_m,H0_12_obs_m,delta1=delta1,delta2=delta2,delta3=delta3)
    E1_m_plus1                     <- posterior_expectations_m_plus1[["E1"]]
    E2_m_plus1                     <- posterior_expectations_m_plus1[["E2"]]
    env01[["gamma_m"]]             <- env02[["gamma_m"]]<- E1_m_plus1
    env12[["gamma_m"]]             <- E1_m_plus1[delta1==1]

    # estimate sigma
    sigma_m_plus1 <- initial_sigma

    ## print results of the first iteration
    if (print==T){
    print(paste("m=",m))
    print(paste("b01=",round(beta_hat_01_m,6)))
    print(paste("b02=",round(beta_hat_02_m,6)))
    print(paste("b12=",round(beta_hat_12_m,6)))
    print(paste("s=",round(sigma_m_plus1,6)))
    }

    ##### Next iterations ####

    # define convergence indicators
    conv_beta01 <- conv_beta02 <- conv_beta12 <- conv_H0_01 <- conv_H0_02 <- conv_H0_12 <- conv_sigma <-conv_betas<-conv_Hs<- conv<-F

    while(conv==F){
      # define stopping iteration
      if(m == stop_iter_num) {
        conv <- F
        break
      }

      init<-Sys.time()
      m<-m+1

      ## estimating Betas
      if (conv_beta01==F) {
      beta_hat_01_m_plus1 <- lbfgs(l_s_m_01_02pert_f.CPP(),grad_beta01_02pert_f.CPP(),env=env01,vars=beta_hat_01_m,invisible = T)$par
      } else  {beta_hat_01_m_plus1 <-beta_hat_01_m}

      if (conv_beta02==F) {
      beta_hat_02_m_plus1 <- lbfgs(l_s_m_01_02pert_f.CPP(),grad_beta01_02pert_f.CPP(),env=env02,beta_hat_02_m,invisible = T)$par
      } else {beta_hat_02_m_plus1 <-beta_hat_02_m}

      if (conv_beta01==T&conv_beta02==T) {
      l12 <- lbfgs(l_s_m_12pert_f.CPP(),grad_beta12pert_f.CPP(),env=env12,beta_hat_12_m,invisible = T,linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING")
      beta_hat_12_m_plus1 <- l12$par
      } else {beta_hat_12_m_plus1 <-beta_hat_12_m}

      #estimating H0s
      start<-Sys.time()
      if ((conv_beta01&conv_Hs)==F) {
        H0_01_obs_m_plus1 <- H0_hat_01_02_observed_perturbed_f(beta=beta_hat_01_m_plus1,gamma=c(E1_m_plus1),a_n=a_n1,X=X01,V=V,delta=delta1,G=G1)
      } else {H0_01_obs_m_plus1 <-H0_01_obs_m}
      end<-Sys.time()
      #print(end-start)

      start<-Sys.time()
      if ((conv_beta02&conv_Hs)==F) {
        H0_02_obs_m_plus1<- H0_hat_01_02_observed_perturbed_f(beta=beta_hat_02_m_plus1,gamma=c(E1_m_plus1),a_n=a_n2,X=X02,V=V,delta=delta2,G=G2)
      } else {H0_02_obs_m_plus1<-H0_02_obs_m}
      end<-Sys.time()
      #print(end-start)

      start<-Sys.time()
      if ((conv_beta12&conv_Hs)==F) {
        H0_12_obs_m_plus1<-H0_hat_12_observed_perturbed_f_new (beta=beta_hat_12_m_plus1,gamma=c(E1_m_plus1),a_n12=a_n12,X12=X12,V=V,W=W,delta1 = delta1,delta3 = delta3,G=G12)
      } else {H0_12_obs_m_plus1<-H0_12_obs_m}
      end<-Sys.time()

      ## computing E1 and E2
      posterior_expectations_m_plus2 <- posterior_expectations(sigma=sigma_m_plus1,H0_01_obs_m_plus1,H0_02_obs_m_plus1,H0_12_obs_m_plus1,delta1=delta1,delta2=delta2,delta3=delta3)
      E1_m_plus2                     <- posterior_expectations_m_plus2[["E1"]]
      E2_m_plus2                     <- posterior_expectations_m_plus2[["E2"]]
      env01[["gamma_m"]]             <- env02[["gamma_m"]] <- E1_m_plus2
      env12[["gamma_m"]]             <- E1_m_plus2[delta1==1]

      sumE1_m_plus2 <- posterior_expectations_m_plus2[["sumE1"]]
      sumE2_m_plus2 <- posterior_expectations_m_plus2[["sumE2"]]

      ## estimating sigma
      sigma_m_plus2 <- optim(par=sigma_m_plus1,
                             fn=function(x) {-( (1/x-1)*(t(G)%*%E2_m_plus2)/n-(1/x)*(t(E1_m_plus2)%*%G)/n+(-log(gamma(1/x))+(1/x)*log(1/x))*mean(G) )},
                             method="L-BFGS-B",lower=0.0075)$par

      ## check convergence
      #betas' convergence
      if ((sum(abs(beta_hat_01_m_plus1-beta_hat_01_m)<conv_betas_bound))==length(beta_hat_01_m)) {conv_beta01 <- T}
      if ((sum(abs(beta_hat_02_m_plus1-beta_hat_02_m)<conv_betas_bound))==length(beta_hat_02_m)) {conv_beta02 <- T}
      if(conv_beta01==T&conv_beta02==T){
        if ((sum(abs(beta_hat_12_m_plus1-beta_hat_12_m)<conv_betas_bound))==length(beta_hat_12_m)) {conv_beta12 <- T}
      }

      #H0s' convergence
      if (mean(abs(H0_01_obs_m_plus1-H0_01_obs_m))<conv_Hs_bound) {conv_H0_01 <- T}
      if (mean(abs(H0_02_obs_m_plus1-H0_02_obs_m))<conv_Hs_bound) {conv_H0_02 <- T}
      if (mean(abs(H0_12_obs_m_plus1-H0_12_obs_m))<conv_Hs_bound) {conv_H0_12 <- T}
      conv_Hs<-conv_H0_01&conv_H0_02&conv_H0_12

      #sigmas' convergence
      if (abs(sigma_m_plus2-sigma_m_plus1)<conv_sigma_bound) {conv_sigma <- T}

      #define when the simulation finally converged
      conv <- conv_sigma&conv_beta01&conv_beta02&conv_beta12&conv_Hs

      ## print the results of current iteration
      if (print==T){
      itter_time <- round(Sys.time()-init,2)
      print(paste("m=",m))
      print(paste("b01=",round(beta_hat_01_m_plus1,6)))
      print(paste("b02=",round(beta_hat_02_m_plus1,6)))
      print(paste("b12=",round(beta_hat_12_m_plus1,6)))
      print(paste("s=",round(sigma_m_plus2,6)))
      print(paste("conv_betas=",conv_beta01,conv_beta02,conv_beta12,"   ","conv_H0s=",conv_H0_01,conv_H0_02,conv_H0_12,"    conv_sigma=",conv_sigma,"    sigma_diff=",round(abs(sigma_m_plus2-sigma_m_plus1),6)))
      print(paste("time=",itter_time,"  sumE1=",round(sumE1_m_plus2,2),"  sumE2=",round(sumE2_m_plus2,2)))
      }

      ## save iterations' results
      beta_hat_01_m    <- beta_hat_01_m_plus1
      beta_hat_02_m    <- beta_hat_02_m_plus1
      beta_hat_12_m    <- beta_hat_12_m_plus1
      H0_01_obs_m      <- H0_01_obs_m_plus1
      H0_02_obs_m      <- H0_02_obs_m_plus1
      H0_12_obs_m      <- H0_12_obs_m_plus1
      E1_m_plus1       <- E1_m_plus2
      sigma_m_plus1old <- sigma_m_plus1
      sigma_m_plus1    <- sigma_m_plus2
    }

      sigma_pert0      <- sigma_m_plus1
      b01_pert0        <- beta_hat_01_m_plus1
      b02_pert0        <- beta_hat_02_m_plus1
      b12_pert0        <- beta_hat_12_m_plus1
      boot_time0       <- difftime(Sys.time(), aa, units = "hours")
      iter0            <- m

      sigma_pert       <- cbind(sigma_pert,sigma_pert0)
      b01_pert         <- cbind(b01_pert,b01_pert0)
      b02_pert         <- cbind(b02_pert,b02_pert0)
      b12_pert         <- cbind(b12_pert,b12_pert0)

      boot_time0       <- difftime(Sys.time(), aa, units = "hours")
      boot_time        <- cbind(boot_time,boot_time0)
      iter             <- cbind(iter,iter0)

      print(paste(i,length(sigma_pert)))
      pert_all           <- rbind(pert_all,cbind(sigma_pert0,t(b01_pert0),t(b02_pert0),t(b12_pert0),boot_time0,iter0,i))
      colnames(pert_all) <- c("sigma",paste0("beta01-",colnames(X01)),paste0("beta02-",colnames(X02)),paste0("beta12-",colnames(X12)),
                              "boot_time","m","boots_i")
      i=i+1

  }

  se_sigma <- sqrt(diag(cov(t(sigma_pert))))
  se_b01   <- sqrt(diag(cov(t(b01_pert))))
  se_b02   <- sqrt(diag(cov(t(b02_pert))))
  se_b12   <- sqrt(diag(cov(t(b12_pert))))

  pert_se  <- cbind(se_sigma, t(se_b01), t(se_b02),t(se_b12))
  colnames(pert_se) <- c("sigma",paste0(colnames(X01),"-01"),paste0(colnames(X02),"-02"),paste0(colnames(X12),"-12"))
  } else {pert_se=NA}

  est        <- est_par[1:n_pars]
  names(est) <- c("sigma",paste0(colnames(X01),"-01"),paste0(colnames(X02),"-02"),paste0(colnames(X12),"-12"))
  return(list(ests = cbind(coefficient=est,ESE=c(pert_se)),
              H0s = list(H0_01_obs_m_plus1, H0_02_obs_m_plus1, H0_12_obs_m_plus1)))
}

#function for estimation without frailty
estimation_without_frailty<-function(X01,X02,X12,V,W,delta1,delta2,delta3,zeta_beta=50,zeta_h=0.01,B,print)
{
  n_pars<-1+dim(X01)[2]+dim(X02)[2]+dim(X12)[2]
  n <- length(V)
  D <- delta1+delta2+delta3
  initial_gamma <- rep(1,n)

  ## bandwidths computation
  ab01  <-ab02 <- ab12 <- zeta_beta
  ah01  <-ah02 <- ah12 <- zeta_h
  band  <- bandwidths(V=V,W=W,delta1=delta1,delta2=delta2,delta3=delta3,ah01=ah01,ah02=ah02,ah12=ah12)
  a_n1  <- band[1,]
  a_n2  <- band[2,]
  a_n12 <- band[3,]

  #estimation
 #first iteration
m=1
gamma_m <- initial_gamma

####set envitronments for cpp functions####
#environment for optimization of l01
env01             <- new.env()
env01[["X"]]      <- X01
env01[["X0"]]     <- X01[delta1==1,]
env01[["delta"]]  <- delta1
env01[["V"]]      <- V
env01[["V0"]]     <- V[delta1==1]
env01[["a"]]      <- a_n1*ab01
#environment for optimization of l02
env02             <- new.env()
env02[["X"]]      <- X02
env02[["X0"]]     <- X02[delta2==1,]
env02[["delta"]]  <- delta2
env02[["V"]]      <- V
env02[["V0"]]     <- V[delta2==1]
env02[["a"]]      <- a_n2*ab02
#environment for optimization of l12
env12             <- new.env()
env12[["X"]]      <- X12[delta1==1,]
env12[["X0"]]     <- X12[delta3==1,]
env12[["delta3"]] <- delta3[delta1==1]
env12[["V"]]      <- V[delta1==1]
env12[["W"]]      <- W[delta1==1]
env12[["W0"]]     <- W[delta3==1]
env12[["a"]]      <- a_n12*ab12
env12[["n"]]      <- n

env01[["gamma_m"]]             <- env02[["gamma_m"]]<- gamma_m
env12[["gamma_m"]]             <- gamma_m[delta1==1]
#########
# estimating Betas
beta_hat_01_m <- lbfgs(l_s_m_01_02_f.CPP(),grad_beta01_02_f.CPP(),vars=rep(0,dim(X01)[2]), environment =env01,invisible = 1)$par
beta_hat_02_m <- lbfgs(l_s_m_01_02_f.CPP(),grad_beta01_02_f.CPP(),vars=rep(0,dim(X02)[2]), environment =env02,invisible = 1)$par
beta_hat_12_m <- lbfgs(l_s_m_12_f.CPP(),grad_beta12_f.CPP(),vars=rep(0,dim(X12)[2]), environment =env12,invisible = 1,linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING")$par

# estimating H0s
H0_01_obs_m  <- H0_hat_01_02_observed_f(beta=beta_hat_01_m,gamma=gamma_m,a_n=a_n1,X=X01,V=V,delta=delta1)
H0_02_obs_m  <- H0_hat_01_02_observed_f(beta=beta_hat_02_m,gamma=gamma_m,a_n=a_n2,X=X02,V=V,delta=delta2)
H0_12_obs_m  <- H0_hat_12_observed_f_new(beta=beta_hat_12_m,gamma=c(gamma_m),a_n12=a_n12,X12=X12,V=V,W=W,delta1 = delta1,delta3 = delta3)

#print results of the first iteration
print(paste("b01=",round(beta_hat_01_m,10)))
print(paste("b02=",round(beta_hat_02_m,10)))
print(paste("b12=",round(beta_hat_12_m,10)))
#print(paste("s=",round(sigma_m_plus1,10)))

##### Next iterations ####

beta01_conv<-beta_hat_01_m
beta02_conv<-beta_hat_02_m
beta12_conv<-beta_hat_12_m
tau_conv_nf<-list(sigma_conv=NA,beta01_conv=beta01_conv,beta02_conv=beta02_conv,beta12_conv=beta12_conv,
                  E1=gamma_m,E2=rep(0,n),
                  H01_conv=H0_01_obs_m,H02_conv=H0_02_obs_m,H12_conv=H0_12_obs_m)

#function for bootstrap without frailty
if (B>0){
if(print==T){print("Bootstrap process begun")}

pert_all_nf<-sigma_pert<-b01_pert<-b02_pert<-b12_pert<-boot_time<-iter<-h0_01<-h0_02<-h0_12<-H0_01<-H0_02<-H0_12<-NULL
i=1
while (max(0,dim(b01_pert)[2])<B){
  aa<-Sys.time()

  ####set environments for cpp functions####
  #environment for optimization of l01
  env01             <- new.env()
  env01[["X"]]      <- X01
  env01[["X0"]]     <- X01[delta1==1,]
  env01[["delta"]]  <- delta1
  env01[["V"]]      <- V
  env01[["V0"]]     <- V[delta1==1]
  env01[["a"]]      <- a_n1*ab01
  #environment for optimization of l02
  env02             <- new.env()
  env02[["X"]]      <- X02
  env02[["X0"]]     <- X02[delta2==1,]
  env02[["delta"]]  <- delta2
  env02[["V"]]      <- V
  env02[["V0"]]     <- V[delta2==1]
  env02[["a"]]      <- a_n2*ab02
  #environment for optimization of l12
  env12             <- new.env()
  env12[["X"]]      <- X12[delta1==1,]
  env12[["X0"]]     <- X12[delta3==1,]
  env12[["delta3"]] <- delta3[delta1==1]
  env12[["V"]]      <- V[delta1==1]
  env12[["W"]]      <- W[delta1==1]
  env12[["W0"]]     <- W[delta3==1]
  env12[["a"]]      <- a_n12*ab12
  env12[["n"]]      <- n

  #########
  G         <-rexp(n,1)
  env01$G   <-env02$G<-G
  env01$G0  <-G[delta1==1]
  env02$G0  <-G[delta2==1]
  env12$G   <-G[delta1==1]
  env12$G0  <-G[delta3==1]

  m       <- 1
  gamma_m <- rep(1,n)
  sigma_m <- NA
  env01$gamma_m <- env02$gamma_m<-gamma_m
  env12$gamma_m <- gamma_m[delta1==1]

    beta_hat_01_m <-lbfgs(l_s_m_01_02pert_f.CPP(),grad_beta01_02pert_f.CPP(),env=env01,vars=tau_conv_nf$beta01_conv,invisible = T)$par
    beta_hat_02_m <-lbfgs(l_s_m_01_02pert_f.CPP(),grad_beta01_02pert_f.CPP(),env=env02,vars=tau_conv_nf$beta02_conv,invisible = T)$par
    beta_hat_12_m <-lbfgs(l_s_m_12pert_f.CPP(),grad_beta12pert_f.CPP(),vars=tau_conv_nf$beta12_conv, env=env12,invisible = 1,linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING")$par

    # print results of the first iteration
    print(paste("m=",m))
    print(paste("b01=",round(beta_hat_01_m,10)))
    print(paste("b02=",round(beta_hat_02_m,10)))
    print(paste("b12=",round(beta_hat_12_m,10)))

    sigma_pert0 <- NA
    b01_pert0   <- beta_hat_01_m
    b02_pert0   <- beta_hat_02_m
    b12_pert0   <- beta_hat_12_m
    boot_time0  <- difftime(Sys.time(), aa, units = "hours")
    iter0       <- m

    sigma_pert  <- cbind(sigma_pert,sigma_pert0)
    b01_pert    <- cbind(b01_pert,b01_pert0)
    b02_pert    <- cbind(b02_pert,b02_pert0)
    b12_pert    <- cbind(b12_pert,b12_pert0)
    boot_time0  <- difftime(Sys.time(), aa, units = "hours")
    boot_time   <- cbind(boot_time,boot_time0)
    iter        <- cbind(iter,iter0)

    print(iter0)
    print(paste(i,dim(b01_pert)[2]))
    pert_all_nf<-rbind(pert_all_nf,cbind(sigma_pert0,t(b01_pert0),t(b02_pert0),t(b12_pert0),boot_time0,iter0,i))
    colnames(pert_all_nf)<-c("sigma",paste0("beta01-",colnames(X01)),
                             paste0("beta02-",colnames(X02)),
                             paste0("beta12-",colnames(X12)),
                             "boot_time","m","boots_i")
    #write.csv(pert_all,file=paste0(dir,"bootstrap_all_no_frailty",".csv"))
    i=i+1
}

se_sigma  <- sqrt(diag(cov(t(sigma_pert))))
se_b01    <- sqrt(diag(cov(t(b01_pert))))
se_b02    <- sqrt(diag(cov(t(b02_pert))))
se_b12    <- sqrt(diag(cov(t(b12_pert))))

se_all_nf <- cbind(se_sigma, t(se_b01), t(se_b02),t(se_b12))
colnames(se_all_nf) <- c("sigma",paste0("beta01-",colnames(X01)),paste0("beta02-",colnames(X02)),paste0("beta12-",colnames(X12)))

} else {se_all_nf=NA}

result <- cbind(coefficient=c(unlist(tau_conv_nf)[1:n_pars]),ESE=c(se_all_nf))
rownames(result) <- c("sigma",paste0("beta01-",colnames(X01)),paste0("beta02-",colnames(X02)),paste0("beta12-",colnames(X12)))
return( result <-cbind(coefficient=c(unlist(tau_conv_nf)[1:n_pars]),ESE=c(se_all_nf))
)
}
