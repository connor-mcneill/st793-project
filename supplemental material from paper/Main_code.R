## This is the code for Section 5.1, Appendix D.1 and D.2 of the Supplementary material

## Parallel computing 
#######################################################################
 
args <- commandArgs(trailingOnly=TRUE)

core <- as.numeric(args[1])

#######################################################################
 
library(glmnet)

## MCL test: 
# X: n x p matrix of covariates.
# Y: n dimensional vector.
#  nfolds: number of folds used for cross-validation of Lasso/logistic Lasso, default = 5.
#  lambda: tuning parameter for the logistic lasso, default = NA. If specified, nfolds will be overwritten.
#  tune.1, tune.2 are tuning parameters in the nodewise lasso. Should be values between 0 and 1.
#  fdr: false discovery rate for multiple testing.
#  please refer to Ma, Cai, and Li (2020) for details.
logistic.test <- function(X,Y,nfolds=5, lambda = 0, 
                          tune.1 = 1.5, tune.2 = 1, 
                          intercept = F, fdr = 0.05){
  f = function(x){
    exp(x)/(1+exp(x))
  }
  p = dim(X)[2]
  n = dim(X)[1]
  
  if(lambda == 0){
    logistic.cv = cv.glmnet(x = X, y = Y, family = "binomial", alpha = 1,  
                            intercept=intercept, nfolds = nfolds, type.measure = "class")
    lambda = logistic.cv$lambda.min
  }
  
  my.logistic.fit = glmnet(x = X, y = Y, family = "binomial", alpha = 1,  
                           intercept=intercept, lambda = lambda, standardize = F)
  
  b.hat = coef(my.logistic.fit)
   
  W.n1 = c(exp(X%*% as.matrix(b.hat)[-1,])/(1+exp(X%*% as.matrix(b.hat)[-1,]))^2)^(-1)
  zeta.try = matrix(nrow =5,ncol = p)
  tau.try = matrix(nrow =5,ncol = p)
  
  V = matrix(ncol=n, nrow = p)
  tau = c()
  for(i in 1:p){
    nodewise.try = glmnet(x= X[,-i], y = X[,i], family = "gaussian", alpha = 1, intercept = F, nlambda = 5, standardize = F)
    for(lambda.i in 1:5){
      V[i,] = X[,i]-X[,-i]%*%nodewise.try$beta[,lambda.i]
      zeta.try[lambda.i,i] = max(abs(as.vector(V[i,]%*%X[,-i])/sqrt(sum((V[i,])^2*W.n1))))
      tau.try[lambda.i,i] = sqrt(sum((V[i,])^2*W.n1))/(V[i,]%*% X[,i])
    }
    zeta0 = sqrt(2*log(p))
    if(min(zeta.try[,i])>sqrt(2*log(p))) zeta0 = tune.1*min(zeta.try[,i])
    lambda.chosen = order(nodewise.try$lambda[which(zeta.try[,i] <= zeta0)], decreasing=T)[1]
    tau[i] = tau.try[lambda.chosen,i]
    lambda.chosen = order(nodewise.try$lambda[tau.try[,i]<=tune.2*tau[i]],decreasing = F)[1]
    tau[i] = tau.try[lambda.chosen,i]
    V[i,] = X[,i]-X[,-i]%*%nodewise.try$beta[,lambda.chosen]
  }

  V2 = t((t(V)*W.n1))
  b.check = c()
  for(j in 1:p){
    b.check[j] = b.hat[j+1]+(V2[j,]%*%(Y-f(X %*% as.matrix(b.hat)[-1])))/(V[j,] %*% X[,j])
  }
  test.value = max((b.check/tau)^2)
  return(test.value)
}


## GC test:
##  X:  n x p matrix of covariates.
##  Y:  n dimensional vector.
##  g.assume, g.assume.deri:  model assumption for mean and its first derivative.
##  V.assume:  model assumption for variance.
##  please refer to Guo and Chen (2016) for details.
CSX.test<-function(X,Y, g.assume, g.assume.deri, V.assume){
  n = nrow(X)
  p = ncol(X)
  mu.0 = rep(g.assume(0), n)  
  var.ratio=g.assume.deri(0)/ V.assume(0)
  Phi.0 = rep(var.ratio, n)            
  Un.tempt = 0
  tr.Var.tempt = 0
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      tempt.value = NULL
      tempt.value = (Y[i]- mu.0[i]) * (Y[j]- mu.0[j]) * Phi.0[j] * Phi.0[i] * crossprod(X[i,],X[j,])
      Un.tempt = Un.tempt + tempt.value
      tr.Var.tempt = tr.Var.tempt + tempt.value^2
    }
  }
  Un = NULL
  Un = (1/n) * 2 * Un.tempt
  tr.Var = NULL
  tr.Var = (1/(n*(n-1))) * 2 * tr.Var.tempt
  CSX.test = Un/sqrt(2*tr.Var) 
  return (CSX.test)
}


## RP test
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
 

##  n, p, m : data dimensions.
##  beta 
##  Gamma: required by generating X.
##  data.gen:  for data generation.
##  g.assume, g.assume.deri:  model assumption for mean and its first derivative.
##  V.assume: model assumption for variance.
##  rho: projection ratio.
simulation<-function(n,p,m,RP.testornpt=TRUE,CSX.testornot=TRUE,TC.testornot=TRUE,beta,Gamma,data.gen,
                      g.assume,g.assume.deri,V.assume,
                     rho,sign=0.95,index){
  k=ceiling(rho*n)
  count.RP=0
  count.CSX=0
  count.TC=0
  L=length(index)
  record.RP=rep(NA,L)
  criti_Gambel= 2*log(p)-log(log(p))-log(pi)-2*log(log(1/(sign)))
  
  for(i in 1:L){
    set.seed(index[i])
    data=NULL
    data=data.gen(n,m,p,Gamma,beta)
    X=as.matrix(data$X)
    Y=as.vector(data$Y)
    
    if ( RP.testornpt==FALSE ){
      count.RP=NA
    } else {
      RP.value = RP.test(X,Y,rho)
      if ( RP.value > qnorm(sign,0,1) ){
        count.RP = count.RP + 1
      }
      if ( sum(beta^2)==0 ){
        record.RP[i]=RP.value
      }}
    
    if ( CSX.testornot==FALSE ){
      count.CSX=NA
    } else {
      CSX.value = CSX.test(X,Y, g.assume, g.assume.deri, V.assume)
      if ( CSX.value > qnorm(sign,0,1) ){
        count.CSX = count.CSX + 1
      }}
    
    if (TC.testornot == FALSE){
      count.TC = NA
    }else{
      TC.value = logistic.test(X,Y, nfolds=5, lambda = 0, 
                               tune.1 = 1.5, tune.2 = 1.01, 
                               intercept = F, fdr = 0.05)
      if ( TC.value > criti_Gambel){
        count.TC = count.TC + 1}}
  }
  count=list(count.RP=count.RP,count.CSX=count.CSX,count.TC=count.TC)
  return(list(count=count,record.RP=record.RP))
}

##----------------------------------------------------
## Data generation ----------------------------------

## Generate data following the logistic model --------
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

data.bino<-function(n,m,p, Gamma, beta){
  Z = NULL;   Z = matrix(rnorm(n*m,0,1), n, m) 
  ## Distribution of Z:
  ##   Uniform case:    Z = matrix(runif(n*m,-sqrt(3),sqrt(3)), n, m)   
  ##   Mixture type:    
  #               Z1 = matrix(sqrt(3/5)*rt((n*(m/2)),5,ncp=0), n, (m/2))   
  #               Z2 = matrix(((rbinom((n*(m/2)), 1, prob=0.5)-0.5)*2), n, (m/2))
  #               Z=cbind(Z1,Z2)
  X = NULL;   X = Z %*% Gamma 
  pro = NULL; pro = logi(X %*% beta)
  Y = NULL;   Y = rbinom(n, 1, prob=pro)
  
  return(list(X=X,Y=Y)) }



## Generate data following Poisson model ------- 
data.pois<-function(n,m,p,Gamma,beta){
  Z = NULL;   Z = matrix(rnorm(n*m,0,1), n, m) 
  ## Distribution of Z:
  ##   Uniform case:    Z = matrix(runif(n*m,-sqrt(3),sqrt(3)), n, m)   
  ##   Mixture type:    
  #               Z1 = matrix(sqrt(3/5)*rt((n*(m/2)),5,ncp=0), n, (m/2))   
  #               Z2 = matrix(((rbinom((n*(m/2)), 1, prob=0.5)-0.5)*2), n, (m/2))
  #               Z=cbind(Z1,Z2)
  X = NULL;   X = Z %*% Gamma 
  pro = NULL; pro = as.vector(exp(X %*% beta))
  Y = NULL;   Y = rpois(n, lambda=pro)
  
  return(list(X=X,Y=Y)) 
}

## Generate data following Model 1 ------- 
data.model1<-function(n,m,p,Gamma, beta){
  
  Z1 = NULL;  Z1 = matrix(sqrt(3/5)*rt((n*(m/2)),5,ncp=0), n, (m/2)) 
  Z2 = NULL;  Z2 = matrix(((rbinom((n*(m/2)), 1, prob=0.5)-0.5)*2), n, (m/2))
  Z = NULL;   Z=cbind(Z1,Z2)
  X = NULL;   X = Z %*% Gamma 
  e = NULL;   e = sqrt(3/5)*rt(n,5,ncp=0)
  Y = NULL;   Y = sin(X %*% beta)+e  
  
  return(list(X=X,Y=Y)) }



## Generate data following Model 2 ------- 
data.model2<-function(n,m,p,Gamma,beta){
  
  Z1 = NULL;  Z1 = matrix(sqrt(3/5)*rt((n*(m/2)),5,ncp=0), n, (m/2))
  Z2 = NULL;  Z2 = matrix(((rbinom((n*(m/2)), 1, prob=0.5)-0.5)*2), n, (m/2))
  Z = NULL;   Z=cbind(Z1,Z2)
  X = NULL;   X = Z %*% Gamma 
  e = NULL;   e = sqrt(3/5)*rt(n,5,ncp=0)
  Y = NULL;   Y = (X %*% beta)^2+X %*% beta+e 
  
  return(list(X=X,Y=Y)) }


## Generate data following Model 3 ------- 
data.model3<-function(n,m,p,Gamma,beta){
  
  Z1 = NULL;  Z1 = matrix(sqrt(3/5)*rt((n*(m/2)),5,ncp=0), n, (m/2)) 
  Z2 = NULL;  Z2 = matrix(((rbinom((n*(m/2)), 1, prob=0.5)-0.5)*2), n, (m/2))
  Z = NULL;   Z=cbind(Z1,Z2)
  X = NULL;   X = Z %*% Gamma 
  e = NULL;   e = sqrt(3/5)*rt(n,5,ncp=0)
  Y = NULL;   Y = (X %*% beta)^2+sin(X %*% beta)+cos(X %*% beta)+e   
  
  return(list(X=X,Y=Y)) }


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


##-------------------------------------------
n=400;p=1000 # or (n,p)=(600,3000)
setting_tempt=setting(n,p,no_block=1, sp_small=10, omega2=0.4, seed=70)
## To generate different Sigma:
##    Sigma1: no_block=1; Sigma2: no_block=100; Sigma3: no_block=5;
## To generate beta with different b^2:
##    for the logistic model: omega2 = 0.4 or 0.8
##    for Poission model: omega2 = 0.1 or 0.2
##    for models 1-3: omega2 = 0.1 or 0.2

Gamma=as.matrix(setting_tempt$Gamma)
m=as.numeric(setting_tempt$m)
beta_sp=as.vector(setting_tempt$beta_sp)  
beta_ran=as.vector(setting_tempt$beta_ran)
index=(20*core-19):(20*core)  


##----------Task 1: Comparison in the table--------------------

true_m=matrix(NA,3,3)
true_null=rep(NA,length(index))

simu_null=simulation(n,p,m,beta=rep(0,p),RP.testornpt=TRUE,CSX.testornot=TRUE,TC.testornot=TRUE,Gamma,data.gen=data.bino,
                    V.assume=logi.deri,g.assume=logi,g.assume.deri=logi.deri,
                     rho=0.4,sign=0.95,index)
simu_sp=simulation(n,p,m,beta=beta_sp,RP.testornpt=TRUE,CSX.testornot=TRUE,TC.testornot=TRUE,Gamma,data.gen=data.bino,
                   V.assume=logi.deri,g.assume=logi,g.assume.deri=logi.deri,
                   rho=0.4,sign=0.95,index)
simu_ran=simulation(n,p,m,beta=beta_ran,RP.testornpt=TRUE,CSX.testornot=TRUE,TC.testornot=TRUE,Gamma,data.gen=data.bino,
                    V.assume=logi.deri,g.assume=logi,g.assume.deri=logi.deri,
                    rho=0.4,sign=0.95,index)
## The default model used is the logistic model.
##   For Poission model:
##    data.gen=data.pois; V.assume=exp; g.assume=exp; g.assume.deri=exp.
##   For Model 1-5:
##    only run RP test:
##       CSX.testornot = FALSE; TC.testornot = FALSE.
##    generate data from  Model 1 as an example :
##       data.gen = data.model1.
##     V.assume, g.assume, g.assume.deri makes no influence in the simulation.
## Change rho for different projection dimension

true_m[1,]=as.numeric(simu_null$count)
true_m[2,]=as.numeric(simu_sp$count)
true_m[3,]=as.numeric(simu_ran$count)
true_null=as.numeric(simu_null$record.RP)  
## this is used to generate the figure of kernel density estimation of the proposed test, see Figure 4 in Appendix D of the Supplementary Material.
 
## save data ------------------------------ 
filename <- paste0("true_m",core,".txt")
write.table(true_m,
            file=filename,
            row.names = FALSE,
            col.names = FALSE)

filename <- paste0("true_null",core,".txt")
write.table(true_null,
            file=filename,
            row.names = FALSE,
            col.names = FALSE)

##-----------Task 2: Comparison in Figure-------------------------
### If one wants to generate figure such as Figure 1 in Section 5.1, please comment the Task 1 and uncomment Task 2. Correpsonding revision in Ser.R is also required.
#
# sign_seq=seq(from=0.01,to=0.99,by=0.04)
# true_array=array(NA,dim=c(3,3,length(sign_seq)))
# 
# 
# for (i in 1:length(sign_seq)){
#   
#   sign=NULL
#   sign=(1-sign_seq[i])
#   
#   # for the logistic model:
#   simu_null=simulation(n,p,m,beta=rep(0,p),RP.testornpt=TRUE,CSX.testornot=TRUE,TC.testornot=TRUE,Gamma,data.gen=data.bino,
#                        V.assume=logi.deri,g.assume=logi,g.assume.deri=logi.deri,
#                        rho=0.4,sign=sign,index)
#   simu_sp=simulation(n,p,m,beta=beta_sp,RP.testornpt=TRUE,CSX.testornot=TRUE,TC.testornot=TRUE,Gamma,data.gen=data.bino,
#                      V.assume=logi.deri,g.assume=logi,g.assume.deri=logi.deri,
#                      rho=0.4,sign=sign,index)
#   simu_ran=simulation(n,p,m,beta=beta_ran,RP.testornpt=TRUE,CSX.testornot=TRUE,TC.testornot=TRUE,Gamma,data.gen=data.bino,
#                       V.assume=logi.deri,g.assume=logi,g.assume.deri=logi.deri,
#                       rho=0.4,sign=sign,index)
#   ## Simulation for Poission model:
#   ##    data.gen=data.pois; V.assume=exp;g.assume=exp;g.assume.deri=exp.
#   
#   ## Simulation for Model 1-5:
#   ##    only run RP test:
#   ##       CSX.testornot = FALSE; TC.testornot = FALSE.
#   ##    generate data from  Model 1 as an example :
#   ##       data.gen = data.model1.
#   ##     V.assume, g.assume, g.assume.deri makes no influence in the simulation.
#   
#   true_array[1,,i]=as.numeric(simu_null$count)
#   true_array[2,,i]=as.numeric(simu_sp$count)
#   true_array[3,,i]=as.numeric(simu_ran$count)
# } 
# 
# ## save data ---------------------------
# filename <- paste0("true_array",core,".txt")
# 
# saveRDS(true_array, 
#         file = filename)



 