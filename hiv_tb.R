library(here)
library(tidyr)
library(lhs)
library(odin)


hiv_tb <- odin::odin({
  nhiv <- user(2)
  n    <- nhiv
  #=====mixing matrix==============
  mix_matrix[,] <- 1  
  Beta[,]        <- beta*mix_matrix[i,j]
  lambda_prod[,] <- Beta[i, j] * prop_infec[j]*I[j]/N[j] 
  lambda[]       <- sum(lambda_prod[i, ])            
  output(lambda) <- TRUE                                
  

  N[]  <- S[i] + E[i] + L[i] + I[i]+ T[i] + R[i]
  totN <- sum(N[])
  
  
  ## on/off version
  scalefac <- (1-exp(-max(t-syear,0)/sut))
  rrcdr1 <- if(t > eyear ) 1.0 else (rracf * scalefac + 1-scalefac)
  rrcdr  <- if(t > syear ) rrcdr1 else 1.0
  
  cdrp <- user()                                      #HIV +ve CDR
  cdrn <- user()                                      #HIV -ve CDR
  cdr[1]   <- rrcdr*cdrn*(recov[1] + mutb[1] +mu[1])/(1-cdrn);
  cdr[2]   <- rrcdr*cdrp*(recov[2] + mutb[2] +mu[2])/(1-cdrp);
  
  dim(cdr) <- n
  
  
  # differential progression
  reactn    <-user() 
  fastn     <-user()
  relapsen  <-user()
  mut       <-user()
  mun       <-user()
  # relative ratios for hivp
  rrMortH   <-user()
  rrprog    <-user()  # irr of hiv progression
  rrInfectH <-user()  # rr of infectiousness in hivp
  rrmu      <-user()
  # HIV neg TB params
  react[1]   <-reactn
  fast[1]    <-fastn
  relapse[1] <-relapsen
  mutb[1]    <-mut
  mu[1]      <-mun
  #HIV +ve TB params
  react[2]   <-reactn*rrprog
  fast[2]    <-fastn*rrprog
  mutb[2]    <-mut*rrMortH
  mu[2]      <-mun*rrmu
  relapse[2] <-relapsen*rrprog
  
  #===recovery============
  recov[1] <- recovn
  recov[2] <-0
  psi[1] <-psin
  psi[2] <-0
  #====prop Infectious======
  pInfecious    <- user() # infectious TB hiv_neg
  
  prop_infec[1] <- pInfecious
  prop_infec[2] <- rrInfectH*pInfecious
  
  
  #=====dimensions=============
  
  dim(fast)  <-n
  dim(react) <-n
  dim(mutb)  <-n
  dim(recov) <-n
  dim(psi)   <-n
  dim(mu)    <-n
  
  
  deriv(S[]) <- I[i] * mutb[i]  + N[i] * mu[i] + T[i]*(1-theta)*tau- S[i] * (lambda[i] + mu[i]) 
  deriv(E[]) <- S[i] * lambda[i] + L[i]*(lambda[i]*(1-psi[i]))- (mu[i]+stab)*E[i]-fast[i]*E[i]  
  deriv(L[]) <- E[i] * stab + I[i] * recov[i]   - (mu[i]  +lambda[i]*(1-psi[i]))*L[i]- react[i]*L[i] 
  deriv(I[]) <- fast[i]*E[i]+ react[i]*L[i]  + R[i]*relapse[i]-I[i]*(mu[i]+mutb[i]+recov[i]+ cdr[i])
  deriv(T[]) <- I[i] * cdr[i] - T[i] * (tau + mu[i]) 
  deriv(R[]) <- theta * (tau)*T[i] - R[i] * ( relapse[i] + mu[i]) 
  
  hivPop  <-user(0.01)   # hiv prevalence
  hivnPop <- 1- hivPop

  #initialisation====
  ari0  <- user(0.01)
  Linit <- 1-exp(-ari0/mu[1])
  
  I0[1] <- (ari0/beta)*hivnPop
  I0[2] <- (ari0/beta)*hivPop
  E0[1] <- (Linit*0.05)*hivnPop
  L0[1] <- (Linit*0.95)*hivnPop
  E0[2] <- (Linit*0.05)*hivPop
  L0[2] <- (Linit*0.95)*hivPop
  T0[] <- I0[i] * 2 /3
  R0[] <- theta * T0[i]
  S0[1] <- hivnPop - E0[1] - L0[1] - I0[1] - T0[1] - R0[1]
  S0[2] <- hivPop - E0[2] - L0[2] - I0[2] - T0[2] - R0[2]
  
  
  initial(S[]) <- S0[i]
  initial(E[]) <- E0[i]
  initial(L[]) <- L0[i]
  initial(I[]) <- I0[i]
  initial(T[]) <- T0[i]
  initial(R[]) <- R0[i]
  
  incfast[]  <- fast[i]*E[i]  
  incslow[]  <- react[i]*L[i] 
  notifall[] <- 1e5*cdr[i]*I[i] 
  pr[]       <- 1e2* incfast[i] /(incfast[i] + incslow[i]+relapse[i]*R[i]) 
  prevall[]  <- 1e5*I[i]/N[i]
  incall[]   <- 1e5*(incfast[i] + incslow[i] + relapse[i]*R[i]) 

  output(prev)<-  hivnPop*prevall[1]+ hivPop*prevall[2]
  output(inc) <-  hivnPop*incall[1] + hivPop*incall[2]
  output(notif)<- hivnPop*notifall[1]+ hivPop*notifall[2]
  output(HpnR) <- prevall[1]/prevall[2]
  output(PR)   <- hivnPop*pr[1]+hivPop*pr[2]
  
  dim(incfast) <- n
  dim(incslow) <- n
  dim(relapse) <- n
  dim(incall)  <- n
  dim(notifall)<- n
  dim(prevall) <- n
  dim(pr)      <- n
  
  #====outputs=========
  output(prevall) <- TRUE
  output(incall)  <- TRUE
  output(notifall)<- TRUE
  output(N)       <- TRUE
  output(totN)    <- TRUE
  
  
  beta        <- user()
  recovn      <- user()             
  stab        <- user()
  psin        <- user()             
  tau         <- user()          
  theta       <- user()          
  rracf       <- user()
  syear       <- user()
  eyear       <- user()
  sut         <- user() #scale-up time
  # dimensions
  dim(mix_matrix)   <- c(n,n)
  dim(Beta)         <- c(n,n)
  dim(lambda_prod)  <- c(n,n)
  dim(lambda) <- n
  dim(prop_infec)<-n
  dim(S)      <- n
  dim(E)      <- n
  dim(L)      <- n
  dim(I)      <- n
  dim(T)      <- n
  dim(R)      <- n
  dim(N)      <- n
  
  dim(S0)     <- n
  dim(E0)     <- n
  dim(L0)     <- n
  dim(I0)     <- n
  dim(T0)     <- n
  dim(R0)     <- n
},target="c")


pms<-list(
  ari0=0.02,          
  hivPop=0.2, 
  beta = 5*2,
  reactn = 0.00102, 
  fastn= 0.0936, 
  stab=1.859, 
  recovn=0.1868,
  mut= 0.2067, 
  psin = 0.7920,    
  relapsen =0.01925, 
  mun=0.01612, 
  tau=2,
  theta=0.9,
  rrprog=8,          
  pInfecious= 0.5,    
  rrInfectH= 0.78,   
  rrMortH = 6,       
  cdrn = 0.6,        
  cdrp = 0.6,        
  rrmu = 3.01,  
  rracf=2.0,
  syear=100,
  eyear=105,
  sut=1) 

# check model 
model <- hiv_tb(user=pms)
t <- seq(0, 200, length.out = 1e3)
y <- data.table(model$run(t))
y[,year:=floor(t)]
y[year==250, 'prev']
ggplot(y, aes(t, prev)) +geom_line()


# parameters for prior distributions
priors <- list(
  ari0k=2,ari0scale=2.5e-2, 
  hivPop_a=1, 
  hivPop_b=9, 
  psi_a = 77.9,
  psi_b= 20.7,
  recov_m =-1.677488,
  recov_s = 0.1830687,
  stab_m =0.6200967,
  stab_s =0.06877732, 
  fast_m = -2.368396,
  fast_s = 0.3221463,
  react_m = -6.886404,
  react_s = 0.5751487,
  relapse_m = -3.95, 
  relapse_s=0.27,
  mutb_m = -1.576473,
  mutb_s=0.08773353,
  cdrn_m=65.31527,
  cdrn_s=32.36644,
  cdrp_m=83.07731,
  cdrp_s=11.54504,
  beta_m=1.678, 
  beta_s=0.371,
  rrprog_a=5.38626,
  rrprog_b=0.53917, 
  rrmut_a=1.8, 
  rrmut_b=0.04306934,
  rrmu_m =1.081835,
  rrmu_s = 0.2005234)

priorsquantiles <- list(
  ari0 = function(x) qgamma(x,scale=priors$ari0scale,
                            shape=priors$ari0k), 
  hivPop = function(x) qbeta(x,priors$hivPop_a,priors$hivPop_b), 
  beta = function(x) qlnorm(x,log(2)+priors$beta_m,priors$beta_s),  #TODO check this
  reactn = function(x) qlnorm(x,priors$react_m,priors$react_s),
  fastn = function(x) qlnorm(x,priors$fast_m,priors$fast_s),
  stab = function(x) qlnorm(x,priors$stab_m,priors$stab_s),
  recovn = function(x) qlnorm(x,priors$recov_m,priors$recov_s),
  mut   = function(x) qlnorm(x,priors$mutb_m,priors$mutb_s),
  psin = function(x) qbeta(x,priors$psi_a,priors$psi_b),
  relapsen = function(x) qlnorm(x,priors$relapse_m,priors$relapse_s),
  rrprog=function(x) qgamma(x,priors$rrprog_a, priors$rrprog_b),
  rrMortH=function(x) qlnorm(x,priors$rrmut_a, priors$rrmut_b),
  cdrn=function(x) qbeta(x,priors$cdrn_m, priors$cdrn_s),
  cdrp=function(x) qbeta(x,priors$cdrp_m, priors$cdrp_s),
  rrmu=function(x) qlnorm(x,priors$rrmu_m, priors$rrmu_s),
  rracf= function(x) qgamma(x,1,rate=.25) + 1, 
  sut = function(x) qgamma(x,1,scale=0.75)
)



pz <- pms
for(nm in names(priorsquantiles)) pz[[nm]] <- priorsquantiles[[nm]](0.5)

pz[['hivPop']] <- 0.2 


# ------------ generate data -------------
nmz <- names(priorsquantiles)

NN <- 1e3
t <- seq(0, 200, length.out = 1e3) 
L <- randomLHS(NN,length(priorsquantiles)) 
PL <- AL <- YAL <- YAL0 <- list() 

for(i in 1:NN){
  if(!i%%50) cat('i = ',i,'\n')

    pz <- pms
  for(j in 1:length(nmz))  
    pz[[nmz[j]]] <- priorsquantiles[[nmz[j]]](L[i,j])
  ## run model
  model$set_user(user=pz)
  y <- model$run(t) #with ACF
  pzz <-pz         

  pz$rracf <- 1.0  #without ACF
  model$set_user(user=pz)
  y0 <- model$run(t)
  
  # compute quantities and capture
  Y  <- as.data.table(y); 
  Y0 <- as.data.table(y0)

  Y  <- Y[t>=pz$syear-10 & t < pz$syear+20]; 
  Y0 <- Y0[t>=pz$syear-10 & t < pz$syear+20];
  Y[,year:=floor(t)]; Y0[,year:=floor(t)]
  #============================

  
  Y[,Q:=1+floor(t*4) %% 4]; Y0[,Q:=1+floor(t*4) %% 4]
  Y[,t:=NULL]; Y0[,t:=NULL]; 
  Y  <- Y[,lapply(.SD,mean),by=.(year,Q)]
  Y0 <- Y0[,lapply(.SD,mean),by=.(year,Q)]
  Y[,id:=i]; Y0[,id:=i];
  YAL[[i]]  <- as.data.table(Y)
  YAL0[[i]] <- as.data.table(Y0)
  PL[[i]]   <- as.data.table(pzz)
  
} 

PL   <- rbindlist(PL)    
YAL  <- rbindlist(YAL)   
YAL0 <- rbindlist(YAL0) 
PL[,id:=1:nrow(PL)]


save(PL,file=here('PL.Rdata'))
save(YAL,file=here('YAL.Rdata'))
save(YAL0,file=here('YAL0.Rdata'))








