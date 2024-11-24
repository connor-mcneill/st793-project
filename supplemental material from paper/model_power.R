## This is the code for Appendix D.2 of the Supplementary material

## Parallel computing
#######################################################################
 
args <- commandArgs(trailingOnly=TRUE)

core <- as.numeric(args[1])

#######################################################################
 
## RP test
#### Input:
##  X:  n x p matrix of covariates.
##  Y:  n dimensional vector.
##  rho: projection ratio.
##  D: random projection times.
RP.test<-function(X,Y,rho,D=10){
  
  n = nrow(X)
  p = ncol(X)
  k = ceiling(rho*n)
  
  I.P1 = diag(n) - (1/n)*matrix(1,n,n)

  P.array = array(NA,dim = c(k,p,D))
  P.array = array(rnorm((k*p*D),0,1),dim=c(k,p,D))
  
  P=NULL
  P=(1/sqrt(p))*as.matrix(apply(P.array,c(1,2),mean))
  U=NULL
  U = P%*%t(X)%*% I.P1
  H=NULL
  H = t(U)%*%solve((U%*%t(U)))%*%U
  
  T.multi = ((t(Y)%*%H%*%Y)/k)/((t(Y)%*%( I.P1-H)%*%Y)/(n-k-1))
  deno = sqrt(2/(n*rho*(1-rho)))
  multi = (T.multi-1)/deno
  
  return(multi)
}
 

#### Input:
##  n, p, m : data dimensions.
##  beta 
##  Gamma: required by generating X.
##  data.gen:  for data generation.
##  rho: projection ratio.
simulation<-function(n,p,m,beta,Gamma,data.gen,
                     rho,sign=0.95,index){
  k=ceiling(rho*n)
  count.RP=0
  L=length(index)
  record.RP=rep(NA,L)

  for(i in 1:L){
    set.seed(index[i])
    data=NULL
    data=data.gen(n,m,p,Gamma,beta)
    X=as.matrix(data$X)
    Y=as.vector(data$Y)
    
    RP.value = RP.test(X,Y,rho)
    if (RP.value>qnorm(sign,0,1)){
    count.RP = count.RP + 1
    }
    if ( sum(beta^2)==0 ){
    record.RP[i]=RP.value
    }}
  
  count=list(count.RP=count.RP)
  return(list(count=count,record.RP=record.RP))
  
}

##----------Data generation ---------------------------------

## Generate data following Model 4 ------- 
data.model4<-function(n,m,p,Gamma,beta){
    
    Z1 = NULL;  Z1 = matrix(sqrt(3/5)*rt((n*(m/2)),5,ncp=0), n, (m/2))  
    Z2 = NULL;  Z2 = matrix(((rbinom((n*(m/2)), 1, prob=0.5)-0.5)*2), n, (m/2))
    Z = NULL;   Z=cbind(Z1,Z2)
    X = NULL;   X = Z %*% Gamma 
    e = NULL;   e = sqrt(3/5)*rt(n,5,ncp=0)
    Y = NULL;   Y = cos(X %*% beta)+e     
    
    return(list(X=X,Y=Y)) }

## Generate data following Model 5 ------- 
data.model5<-function(n,m,p,Gamma,beta){
    
    Z1 = NULL;  Z1 = matrix(sqrt(3/5)*rt((n*(m/2)),5,ncp=0), n, (m/2))   
    Z2 = NULL;  Z2 = matrix(((rbinom((n*(m/2)), 1, prob=0.5)-0.5)*2), n, (m/2))
    Z = NULL;   Z=cbind(Z1,Z2)
    X = NULL;   X = Z %*% Gamma 
    e = NULL;   e = sqrt(3/5)*rt(n,5,ncp=0)
    Y = NULL;   Y = (X %*% beta)^2+e
    
    return(list(X=X,Y=Y)) }


 
setting<-function(n,p,no_block,sp_small,omega2,seed){
  
  Gamma=NULL;m=NULL;beta_sp=NULL;beta_ran=NULL
  set.seed(seed)
  
  s=ceiling(n^1)
  M=ceiling(n^0.8)
  diag.t=rep(0,(p-M))
  v=4
  for (i in 1:(p-M)){
    diag.t[i]=i^(-v)
  }
  sum.v=sum(diag.t)
  dia.tempt=((s-M)/sum.v)*diag.t
  dia=sqrt(c(rep(1,M),dia.tempt))
  dia=sqrt(dia)
  
  if(no_block==1){ 
    save=rnorm(p^2,0,1)
    orignal=matrix(save,p,p)
    symtric=orignal+t(orignal)           
    O=as.matrix(eigen(symtric)$vectors)  
    }else{
    p_sub=NULL
    p_sub=p/no_block
    orignal=matrix(0,p,p)
    for(i in 1:no_block){
      orignal[((i-1)*p_sub+1):(i*p_sub), ((i-1)*p_sub+1):(i*p_sub)]=matrix(rnorm(p_sub^2,0,1),p_sub,p_sub)
    }
    symtric=orignal+t(orignal)           
    O=as.matrix(eigen(symtric)$vectors) 
    }
  
  Gamma=NULL
  Gamma=O%*%diag(dia)%*%t(O) 
  Sigma=Gamma%*%t(Gamma)
  m=p
  
  beta_sp_tempt = rep(0,p)
  beta_sp_tempt[sample(1:p,sp_small)] =  rep(1, sp_small)

  beta_ran_tempt=O[,1:100]%*%rnorm(100)
  
  t1=as.numeric(t(beta_sp_tempt)%*%Sigma%*%beta_sp_tempt)
  t2=as.numeric(t(beta_ran_tempt)%*%Sigma%*%beta_ran_tempt)
  
  beta_sp=sqrt(omega2/t1)* beta_sp_tempt
  beta_ran = sqrt(omega2/t2)*beta_ran_tempt
  
  return(list(Gamma=Gamma,m=m,beta_sp=beta_sp,beta_ran=beta_ran))
}

##--------------------------------------------
n=400;p=1000 

setting_tempt=setting(n,p,no_block=5, sp_small=10, omega2=1, seed=20)
## generate different Sigma:
##    Sigma1: no_block=1; Sigma2: no_block=100; Sigma3: no_block=5;

Gamma=as.matrix(setting_tempt$Gamma)
m=as.numeric(setting_tempt$m)
beta_sp_ori=as.vector(setting_tempt$beta_sp)
beta_ran_ori=as.vector(setting_tempt$beta_ran)

index=(20*core-19):(20*core)
omega_seq=seq(from=0,to=2,by=0.05)
true_array=array(NA,dim=c(2,2,length(omega_seq)))

for (i in 1:length(omega_seq)){

omega2=NULL;   omega2=omega_seq[i]
beta_sp=NULL;  beta_sp=sqrt(omega2)*beta_sp_ori
beta_ran=NULL; beta_ran=sqrt(omega2)*beta_ran_ori

simu_sp2=simulation(n,p,m,beta=beta_sp,Gamma,data.gen= data.model5,rho=0.2,sign=0.95,index)
simu_sp4=simulation(n,p,m,beta=beta_sp,Gamma,data.gen= data.model5,rho=0.4,sign=0.95,index)
simu_ran2=simulation(n,p,m,beta=beta_ran,Gamma,data.gen=data.model5,rho=0.2,sign=0.95,index)
simu_ran4=simulation(n,p,m,beta=beta_ran,Gamma,data.gen=data.model5,rho=0.4,sign=0.95,index)

true_array[1, 1, i]=as.numeric(simu_sp2$count$count.RP)
true_array[1, 2, i]=as.numeric(simu_sp4$count$count.RP)
true_array[2, 1, i]=as.numeric(simu_ran2$count$count.RP)
true_array[2, 2, i]=as.numeric(simu_ran4$count$count.RP)
 }

filename <- paste0("true_array",core,".txt")
 
saveRDS(true_array, 
        file = filename)
 