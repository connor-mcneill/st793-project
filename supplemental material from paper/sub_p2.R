## This is the code for Appendix D.5 of the Supplementary material

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

  UX1=I.P1%*%X1
  PX1=UX1%*%solve(t(UX1)%*%UX1)%*%t(UX1)
  
  X2.P.new=X2%*%t(P.new)
  Z.new=I.P1%*%cbind(X1,X2.P.new)
  PX.new=Z.new%*%solve(t(Z.new)%*%Z.new)%*%t(Z.new)
  T.new=((t(Y)%*%(PX.new-PX1)%*%Y)/k2)/((t(Y)%*%(I.P1-PX.new)%*%Y)/(n-k2-1-p1))
  
  X2.P.multi=X2%*%t(P)
  Z.multi=I.P1%*%cbind(X1,X2.P.multi)
  PX.multi=Z.multi%*%solve(t(Z.multi)%*%Z.multi)%*%t(Z.multi)
  T.multi=((t(Y)%*%(PX.multi-PX1)%*%Y)/k2)/((t(Y)%*%(I.P1-PX.multi)%*%Y)/(n-k2-1-p1))
  
  X2.P.sparse=X2%*%t(P_sparse)
  Z.sparse=I.P1%*%cbind(X1,X2.P.sparse)
  PX.sparse=Z.sparse%*%solve(t(Z.sparse)%*%Z.sparse)%*%t(Z.sparse)
  T.sparse=((t(Y)%*%(PX.sparse-PX1)%*%Y)/k2)/((t(Y)%*%(I.P1-PX.sparse)%*%Y)/(n-k2-1-p1))
  
  deno<-sqrt((2*(1-rho1))/(n*rho2*(1-rho1-rho2)))
  multi<-(T.multi-1)/deno
  new<-(T.new-1)/deno
  sparse<-(T.sparse-1)/deno
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
simulation<-function(n,p1,p2,m,rho2,beta,Gamma,sign=0.95,D,S,index){
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
    data=data.bino(n,m,p,Gamma,beta)
    X=as.matrix(data$X)
    Y=as.vector(data$Y)
    pro<-RP.test(n,p1,p2,X,Y,rho2,D,S)
    modify<-as.numeric(pro$new)
    modify.multi<-as.numeric(pro$multi)
    modify.sparse<-as.numeric(pro$sparse)
    
    if(modify>qnorm(sign,0,1)){
      count.new=count.new+1
    }
    
    if(modify.multi>qnorm(sign,0,1)){
      count.multi=count.multi+1
    }
    
    if(modify.sparse>qnorm(sign,0,1)){
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


##----------Data generation ---------------------------------
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

data.bino<-function(n,m,p,Gamma,beta){
  Z = NULL;   Z = matrix(rnorm(n*m,0,1), n, m)  
  ## Distribution of Z:
  ##   Uniform case:   Z = matrix(runif(n*m,-sqrt(3),sqrt(3)), n, m)   
  ##   Rademacher case:  Z = matrix(sample(c(-1,1),n*m,replace=TRUE,pro=c(0.5,0.5)),n,m) 
  X = NULL;   X = Z %*% t(Gamma)  
  pro = NULL; pro = logi(X %*% beta)
  Y = NULL;   Y = rbinom(n, 1, prob=pro)
  return(list(X=X,Y=Y)) 
}

setting<-function(n,p1,p2,seed,sp,omega2){
    set.seed(seed)
    save1=rnorm(p1^2,0,1)
    orignal1=matrix(save1,p1,p1)
    symtric1=orignal1+t(orignal1)            
    O1=as.matrix(eigen(symtric1)$vectors)    
    diag1=diag(sqrt(abs(rnorm(p1,0,1))))
    Gamma1=O1%*%diag1
    s=ceiling(n^0.8)
    M=ceiling(n^0.72)
    diag.t=rep(0,(p2-M))
    v=4
    for (i in 1:(p2-M)){
    diag.t[i]=i^(-v)
    }
    sum.v=sum(diag.t)
    dia.tempt=((s-M)/sum.v)*diag.t
    dia=sqrt(c(rep(1,M),dia.tempt))
    dia=sqrt(dia)
    p_sub=NULL; p_sub=p2/100
    orignal=matrix(0,p2,p2)
    for(i in 1:100){
      orignal[((i-1)*p_sub+1):(i*p_sub), ((i-1)*p_sub+1):(i*p_sub)]=matrix(rnorm(p_sub^2,0,1),p_sub,p_sub)
    }
    symtric=orignal+t(orignal)           
    O=as.matrix(eigen(symtric)$vectors)
    Gamma2=NULL
    Gamma2=O%*%diag(dia) 
    p=p1+p2
    m=p
    Gamma=matrix(0,p,m)
    Gamma[1:p1,1:p1]= Gamma1
    Gamma[(p1+1):p,(p1+1):m]= Gamma2
    Sigma2=Gamma2%*%t(Gamma2)
    beta_sp2_tempt = rep(0,p2)
    beta_sp2_tempt[sample(1:p2,sp)]=rep(1,sp)
    O2_eigen=as.matrix(eigen(Sigma2)$vectors)
    L=100
    pick=rnorm(L)
    beta_ran2_tempt=O2_eigen[,1:L]%*%pick
    t10=as.numeric(t(beta_sp2_tempt)%*%Sigma2%*%beta_sp2_tempt)
    t20=as.numeric(t(beta_ran2_tempt)%*%Sigma2%*%beta_ran2_tempt)
    beta2_sp=sqrt(omega2/t10)* beta_sp2_tempt
    beta2_ran = sqrt(omega2/t20)*beta_ran2_tempt
    beta1=NULL; beta1_tempt=rnorm(p1)
    beta1=beta1_tempt/sqrt(sum(beta1_tempt^2))
    beta_sp=c(beta1,beta2_sp)
    beta_ran=c(beta1,beta2_ran)
    beta_null=c(beta1,rep(0,p2))
  return(list(Gamma=Gamma,m=m,
              beta_sp=beta_sp,beta_ran=beta_ran,beta_null=beta_null))
}


n=400;p1=40;p2=1000
set=setting(n,p1,p2,sp=10,omega2=0.4,seed=20)  
## omega2=b2^2. take 0.4 or 0.8
Gamma=as.matrix(set$Gamma)
m=as.numeric(set$m)
beta_sp=as.vector(set$beta_sp)
beta_ran=as.vector(set$beta_ran)
beta_null=as.vector(set$beta_null)
index=(20*core-19):(20*core) 

matrix_cor=array(NA,dim = c(2,3,3))

beta=NULL; beta=beta_null
simu_null_cor.1=simulation(n,p1,p2,m,rho2=0.2,beta,Gamma,sign=0.95,D=10,S=400,index)
simu_null_cor.2=simulation(n,p1,p2,m,rho2=0.4,beta,Gamma,sign=0.95,D=10,S=400,index)

matrix_cor[1,,1]=as.numeric(simu_null_cor.1$count)
matrix_cor[2,,1]=as.numeric(simu_null_cor.2$count)

beta=NULL; beta=beta_sp
simu_sp_cor.1=simulation(n,p1,p2,m,rho2=0.2,beta,Gamma,sign=0.95,D=10,S=400,index)
simu_sp_cor.2=simulation(n,p1,p2,m,rho2=0.4,beta,Gamma,sign=0.95,D=10,S=400,index)

matrix_cor[1,,2]=as.numeric(simu_sp_cor.1$count)
matrix_cor[2,,2]=as.numeric(simu_sp_cor.2$count)

beta=NULL; beta=beta_ran
simu_ran_cor.1=simulation(n,p1,p2,m,rho2=0.2,beta,Gamma,sign=0.95,D=10,S=400,index)
simu_ran_cor.2=simulation(n,p1,p2,m,rho2=0.4,beta,Gamma,sign=0.95,D=10,S=400,index)

matrix_cor[1,,3]=as.numeric(simu_ran_cor.1$count)
matrix_cor[2,,3]=as.numeric(simu_ran_cor.2$count)

## save data-----------------------------------------

filename <- paste0("matrix_cor",core,".txt")
saveRDS(matrix_cor, file=filename)
