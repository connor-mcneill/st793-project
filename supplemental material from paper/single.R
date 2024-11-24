## This is the code for the simulation in Appendix D.6 of the Supplementary material

## Parallel computing
#######################################################################
args <- commandArgs(trailingOnly=TRUE)

core <- as.numeric(args[1])
#######################################################################

#### Input:
##  n: sample size
##  p1, p2: dimensions of covariates x1 and x2 respectively. Here x^T=(x1^T,x2^T).
##  X:  n x p matrix of covariates.
##  Y:  n dimensional vector.
##  rho2: projection ratio.
##  D: random projection times.
##  S: parameter of sparse projection matrix
RP.test<-function(n,p1,p2,X,Y,rho2,D,S){
  p=p1+p2
  X1=X[,1:p1]
  X2=X[,(p1+1):p]
  k2=ceiling(rho2*n)
  rho1=p1/n
  P1=(1/n)*matrix(1,n,n)
  I.P1=diag(n)-P1
  P.array=array(NA,dim = c(k2,p2,D))
  for(i in 1:D){
    data.p=NULL
    data.p=rnorm((k2*p2),0,1)
    P.tempt=NULL
    P.tempt=matrix(data.p,nrow=k2,ncol=p2)
    P.array[,,i]=P.tempt
  }
  P=NULL
  P=as.matrix(apply(P.array,c(1,2),mean))
  P.new=NULL
  P.new=as.matrix(P.array[,,1])
  prob=1/sqrt(S)
  data.p=NULL
  data.p=sample(c(-1,0,1),(k2*p2),replace = T,prob =c((1/2)*prob,(1-prob),(1/2)*prob) )
  P_sparse=NULL
  P_sparse=matrix(data.p,nrow=k2,ncol=p2)

  X2.P.new=X2%*%t(P.new)
  Z.new=I.P1%*%cbind(X1,X2.P.new)
  PX.new=Z.new%*%solve(t(Z.new)%*%Z.new)%*%t(Z.new)
  UX2.P.new=I.P1%*%X2.P.new
  PX2.P.new=UX2.P.new%*%solve(t(UX2.P.new)%*%UX2.P.new)%*%t(UX2.P.new)
  T.new=((t(Y)%*%(PX.new-PX2.P.new)%*%Y)/p1)/((t(Y)%*%(I.P1-PX.new)%*%Y)/(n-k2-1-p1))
  
  X2.P.multi=X2%*%t(P)
  Z.multi=I.P1%*%cbind(X1,X2.P.multi)
  PX.multi=Z.multi%*%solve(t(Z.multi)%*%Z.multi)%*%t(Z.multi)
  UX2.P.multi=I.P1%*%X2.P.multi
  PX2.P.multi=UX2.P.multi%*%solve(t(UX2.P.multi)%*%UX2.P.multi)%*%t(UX2.P.multi)
  T.multi=((t(Y)%*%(PX.multi-PX2.P.multi)%*%Y)/p1)/((t(Y)%*%(I.P1-PX.multi)%*%Y)/(n-k2-1-p1))

  X2.P.sparse=X2%*%t(P_sparse)
  Z.sparse=I.P1%*%cbind(X1,X2.P.sparse)
  PX.sparse=Z.sparse%*%solve(t(Z.sparse)%*%Z.sparse)%*%t(Z.sparse)
  UX2.P.sparse=I.P1%*%X2.P.sparse
  PX2.P.sparse=UX2.P.sparse%*%solve(t(UX2.P.sparse)%*%UX2.P.sparse)%*%t(UX2.P.sparse)
  T.sparse=((t(Y)%*%(PX.sparse-PX2.P.sparse)%*%Y)/p1)/((t(Y)%*%(I.P1-PX.sparse)%*%Y)/(n-k2-1-p1))

  multi<-T.multi 
  new<-T.new
  sparse<-T.sparse
  return(list(multi=multi,new=new,sparse=sparse))
}

#### Input:
##  n: sample size
##  p1, p2: dimensions of covariates x1 and x2 respectively. Here x^T=(x1^T,x2^T).
##  m: dimension of Z, required by generating X.
##  rho2: projection ratio.
##  beta 
##  Gamma: required by generating X.
##  D: random projection times.
##  S: parameter of sparse projection matrix
simulation<-function(n,p1,p2,m,rho2,beta,Gamma,sign=0.95,data_gen,D,S,index){
    
  p=p1+p2
  k2=ceiling(rho2*n)
  rho1=p1/n
  
  count.new=0
  count.multi=0
  count.P=0
  
  L=length(index)
  new=rep(NA,L)
  multi=rep(NA,L)
  P=rep(NA,L)
  
  for(i in 1:L){
    seed=index[i]
    set.seed(seed)
    data=NULL
    data=data_gen(n,m,p,Gamma,beta)
    X=as.matrix(data$X)
    Y=as.vector(data$Y)
    
    pro<-RP.test(n,p1,p2,X,Y,rho2,D,S)
    
    modify<-as.numeric(pro$new)
    modify.multi<-as.numeric(pro$multi)
    modify.sparse<-as.numeric(pro$sparse)
    
   if(modify>qf(sign, df1=p1, df2=(n-1-p1-k2), lower.tail = TRUE)){
      count.new=count.new+1
    }
    
    if(modify.multi>qf(sign, df1=p1, df2=(n-1-p1-k2), lower.tail = TRUE)){
      count.multi=count.multi+1
    }
    
    if(modify.sparse>qf(sign, df1=p1, df2=(n-1-p1-k2), lower.tail = TRUE)){
      count.P=count.P+1
    }
    
    new[i]=modify
    multi[i]=modify.multi
    P[i]=modify.sparse
    
  }
  count=list(count.new=count.new,count.multi=count.multi,count.P=count.P)
  value=list(new=new,multi=multi,P=P)
  return(list(count=count,value=value))
}


## Data generation ---------------------------------

## Generate data following the logistic model --------

logi<-function(x){
    y = NULL
    y = exp(x)/(1+exp(x))
    return(y) 
}

data.bino<-function(n,m,p, Gamma,beta){
    Z1 = matrix(sqrt(3/5)*rt((n*(m/2)),5,ncp=0), n, (m/2))   
    Z2 = matrix(((rbinom((n*(m/2)), 1, prob=0.5)-0.5)*2), n, (m/2))
    Z = NULL;   Z=cbind(Z1,Z2)
    X = NULL;   X = Z %*% Gamma 
    pro = NULL; pro = logi(X %*% beta)
    Y = NULL;   Y = rbinom(n, 1, prob=pro)
    return(list(X=X,Y=Y)) }


## Generate data following Poisson model ------- 

data.pois<-function(n,m,p,Gamma,beta){
    Z1 = matrix(sqrt(3/5)*rt((n*(m/2)),5,ncp=0), n, (m/2))   
    Z2 = matrix(((rbinom((n*(m/2)), 1, prob=0.5)-0.5)*2), n, (m/2))
    Z = NULL;   Z=cbind(Z1,Z2)
    X = NULL;   X = Z %*% Gamma 
    pro = NULL; pro = as.vector(exp(X %*% beta))
    Y = NULL;   Y = rpois(n, lambda=pro)
    return(list(X=X,Y=Y)) }


## Generate data following Model 1 ------- 

data.model1<-function(n,m,p,Gamma,beta){
    Z1 = NULL;  Z1 = matrix(sqrt(3/5)*rt((n*(m/2)),5,ncp=0), n, (m/2)) 
    Z2 = NULL;  Z2 = matrix(((rbinom((n*(m/2)), 1, prob=0.5)-0.5)*2), n, (m/2))
    Z = NULL;   Z=cbind(Z1,Z2)
    X = NULL;   X = Z %*% Gamma 
    e = NULL;   e = sqrt(3/5)*rt(n,5,ncp=0)
    Y = NULL;   Y = sin(X %*% beta)+e  
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


##-----------------------------------------
n=400;p1=1;p2=999;p=p1+p2

sign_value=0.95; no_seq=c(1,100,5)
matrix_cor=array(NA,dim = c(3,3,length(no_seq)))
index=(20*core-19):(20*core) 

for (i in 1:length(no_seq)){
  set=setting(n,p,no_block=no_seq[i], sp_small=10, omega2=1, seed=15)
  Gamma=NULL;Gamma=as.matrix(set$Gamma)
  m=as.numeric(set$m)
  
  data_gen=data.bino 
  beta=NULL; beta=rep(0,p)

  simu_m1_1=simulation(n,p1,p2,m,rho2=0.2,beta,Gamma,sign=sign_value,data_gen,D=10,S=400,index)
  simu_m1_2=simulation(n,p1,p2,m,rho2=0.4,beta,Gamma,sign=sign_value,data_gen,D=10,S=400,index)
  simu_m1_3=simulation(n,p1,p2,m,rho2=0.6,beta,Gamma,sign=sign_value,data_gen,D=10,S=400,index)
  
  matrix_cor[1,1,i]=as.numeric(simu_m1_1$count$count.multi)
  matrix_cor[1,2,i]=as.numeric(simu_m1_2$count$count.multi)
  matrix_cor[1,3,i]=as.numeric(simu_m1_3$count$count.multi)
  
  data_gen=data.pois  
  beta=NULL; beta=rep(0,p)
  
  simu_m2_1=simulation(n,p1,p2,m,rho2=0.2,beta,Gamma,sign=sign_value,data_gen,D=10,S=400,index)
  simu_m2_2=simulation(n,p1,p2,m,rho2=0.4,beta,Gamma,sign=sign_value,data_gen,D=10,S=400,index)
  simu_m2_3=simulation(n,p1,p2,m,rho2=0.6,beta,Gamma,sign=sign_value,data_gen,D=10,S=400,index)
  
  matrix_cor[2,1,i]=as.numeric(simu_m2_1$count$count.multi)
  matrix_cor[2,2,i]=as.numeric(simu_m2_2$count$count.multi)
  matrix_cor[2,3,i]=as.numeric(simu_m2_3$count$count.multi)
  
  data_gen=data.model1  
  beta=NULL; beta=rep(0,p)
  
  simu_m3_1=simulation(n,p1,p2,m,rho2=0.2,beta,Gamma,sign=sign_value,data_gen,D=10,S=400,index)
  simu_m3_2=simulation(n,p1,p2,m,rho2=0.4,beta,Gamma,sign=sign_value,data_gen,D=10,S=400,index)
  simu_m3_3=simulation(n,p1,p2,m,rho2=0.6,beta,Gamma,sign=sign_value,data_gen,D=10,S=400,index)
  
  matrix_cor[3,1,i]=as.numeric(simu_m3_1$count$count.multi)
  matrix_cor[3,2,i]=as.numeric(simu_m3_2$count$count.multi)
  matrix_cor[3,3,i]=as.numeric(simu_m3_3$count$count.multi)
  
}

filename <- paste0("matrix_cor",core,".txt")
saveRDS(matrix_cor, file=filename)

