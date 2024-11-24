## This is the code for Appendix D.4 of the Supplementary material

## Parallel computing
#######################################################################

args <- commandArgs(trailingOnly=TRUE)

core <- as.numeric(args[1])

#######################################################################

library(glmnet)

Classical.test<-function(X,Y){
  
  n = nrow(X)
  p = ncol(X)
 
  beta_0=rep(0,p)
  
  esti = glm(Y~ X + 0,  family = binomial(link = "logit"))
  esti_beta = NULL; esti_beta = as.vector(esti$coefficients)
  esti_I_0 = NULL; esti_I_0 = (1/n)*(1/4)*t(X)%*%X
 
  LRT.value = NULL
  LRT.value = esti$null.deviance-esti$deviance
 
  Wald.value = NULL
  Wald.value = t(esti_beta)%*%solve(vcov(esti))%*%esti_beta
  
  Score.value = NULL
  Score.value = (1/n)*t(Y-(1/2))%*%X%*%solve(esti_I_0)%*%t(X)%*%(Y-(1/2))
  
  return(list(LRT.value=LRT.value,Wald.value=Wald.value, Score.value=Score.value))
  
}

 
simulation<-function(n,p,RP.Classical=TRUE,beta,Gamma,data.gen,sign=0.95,index){
  
  count.LRT =0
  count.Wald =0
  count.Score =0
  
  L=length(index)
  record.LRT=rep(NA,L)
  record.Wald =rep(NA,L)
  record.Score=rep(NA,L)
 
  for(i in 1:L){
    set.seed(index[i])
    data=NULL
    data=data.gen(n,p, mu=rep(0,p),Gamma,beta)
    X=as.matrix(data$X)
    Y=as.vector(data$Y)
     
    if (RP.Classical==FALSE){
      count.LRT=NA
      count.Wald=NA
      count.Score=NA
      
    } else {
      Classical.value = NULL
      Classical.value = Classical.test(X,Y)
      LRT.value = as.numeric(Classical.value$LRT.value)
      Wald.value = as.numeric(Classical.value$Wald.value)
      Score.value = as.numeric(Classical.value$Score.value)
      
      if ( LRT.value > qchisq(sign,p) ){
        count.LRT = count.LRT + 1
      }
      if ( Wald.value > qchisq(sign,p) ){
        count.Wald = count.Wald + 1
      }
      if ( Score.value > qchisq(sign,p) ){
        count.Score = count.Score + 1
      }
      
      if ( sum(beta^2)==0 ){
        record.LRT[i]=LRT.value
        record.Wald[i] =Wald.value
        record.Score[i]=Score.value
      }}}
    
  count=list(count.LRT=count.LRT, count.Wald=count.Wald, count.Score=count.Score)
  record=list(record.LRT=record.LRT, record.Wald=record.Wald, record.Score=record.Score)
  
  return(list(count=count, record=record))
  
}

## generate data following the logistic model ------------------------

logi<-function(x){
  y = NULL
  y = exp(x)/(1+exp(x))
  return(y) 
}

logi.deri<-function(x){
  y = NULL
  y = exp(x)/((1+exp(x))^2)
  return(y)
}

data.bino<-function(n,p,mu,Gamma,beta){
  Z = NULL;   Z= matrix(rnorm(n*p), n,p)  
  X = NULL;   X = Z %*% Gamma + tcrossprod(rep(1,n), mu)   
  pro = NULL; pro = logi(X %*% beta)
  Y = NULL;   Y = rbinom(n, 1, prob=pro)
  return(list(X=X,Y=Y)) }
 

##-----------------------------------------------------------------
index=1+(800*core-799):(800*core) 

p_seq=c(1,5,10,50,100,200) ## dimension of p
true_null=array(NA,dim=c(length(p_seq),3,length(index)))

for(i in 1:length(p_seq)){
  n=1000; p=p_seq[i]
  Gamma=diag(1,p)
  simu_null=simulation(n,p,beta=rep(0,p),RP.Classical=TRUE,Gamma,data.gen=data.bino,sign=0.95,index)
  true_null[i,1,]=as.numeric(simu_null$record$record.LRT)
  true_null[i,2,]=as.numeric(simu_null$record$record.Wald)
  true_null[i,3,]=as.numeric(simu_null$record$record.Score)
}

## save data -------------------------------

filename <- paste0("true_null",core,".txt")
saveRDS(true_null, 
        file = filename)
 
 